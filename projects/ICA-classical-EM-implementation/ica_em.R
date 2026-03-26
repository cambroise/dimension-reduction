# ============================================================================
# ICA via EM Algorithm - R Implementation
# ============================================================================
# This implementation uses a mixture of Gaussians model for each source
# (following Attias, 1999) which makes the E-step tractable.
# ============================================================================

library(MASS)  # For mvrnorm

# ============================================================================
# Helper Functions
# ============================================================================

#' Compute log-likelihood for a Gaussian mixture component
#'
#' @param x Vector of observations
#' @param mu Mean of component
#' @param sigma2 Variance of component
#' @return Log-likelihood for each observation
log_gaussian <- function(x, mu, sigma2) {
  -0.5 * log(2 * pi * sigma2) - 0.5 * (x - mu)^2 / sigma2
}

#' Soft-thresholding for sparsity (optional enhancement)
#'
#' @param x Input value
#' @param lambda Threshold parameter
#' @return Thresholded value
soft_threshold <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}

# ============================================================================
# Main ICA-EM Implementation
# ============================================================================

#' ICA via Expectation-Maximization Algorithm
#'
#' @param X Data matrix (d x n): d = dimensions, n = observations
#' @param n_components Number of independent components to extract
#' @param n_mixtures Number of Gaussian mixture components per source
#' @param max_iter Maximum number of EM iterations
#' @param tol Convergence tolerance on log-likelihood
#' @param noise_variance Observation noise variance (fixed)
#' @param verbose Print iteration details
#' @param init_method Initialization method: "random", "pca", or "fastica"
#'
#' @return List containing:
#'   - A: Estimated mixing matrix (d x n_components)
#'   - S: Estimated sources (n_components x n)
#'   - W: Unmixing matrix (pseudo-inverse of A)
#'   - log_likelihood: Log-likelihood at each iteration
#'   - converged: Whether algorithm converged
#'
#' @export
ica_em <- function(X, 
                   n_components = NULL,
                   n_mixtures = 3,
                   max_iter = 100,
                   tol = 1e-6,
                   noise_variance = 0.1,
                   verbose = TRUE,
                   init_method = "pca") {
  
  # ========================================================================
  # 1. Input validation and preprocessing
  # ========================================================================
  
  if (is.null(dim(X))) {
    X <- matrix(X, nrow = 1)
  }
  
  d <- nrow(X)  # Number of observed dimensions
  n <- ncol(X)  # Number of observations
  
  if (is.null(n_components)) {
    n_components <- d
  }
  
  if (n_components > d) {
    stop("n_components cannot exceed number of observed dimensions")
  }
  
  # Center the data
  X_mean <- rowMeans(X)
  X_centered <- X - X_mean
  
  if (verbose) {
    cat(sprintf("ICA-EM Algorithm\n"))
    cat(sprintf("================\n"))
    cat(sprintf("Data dimensions: %d x %d\n", d, n))
    cat(sprintf("Components: %d\n", n_components))
    cat(sprintf("Mixture components per source: %d\n", n_mixtures))
    cat(sprintf("Initialization: %s\n\n", init_method))
  }
  
  # ========================================================================
  # 2. Initialize mixing matrix A
  # ========================================================================
  
  A <- initialize_mixing_matrix(X_centered, n_components, method = init_method)
  
  # ========================================================================
  # 3. Initialize source model parameters
  # ========================================================================
  
  # Each source has n_mixtures Gaussian components
  # Mixture weights (alpha), means (mu), and variances (sigma2)
  source_params <- list()
  for (k in 1:n_components) {
    source_params[[k]] <- list(
      alpha = rep(1/n_mixtures, n_mixtures),  # Equal mixture weights initially
      mu = seq(-2, 2, length.out = n_mixtures),  # Spread means
      sigma2 = rep(1, n_mixtures)  # Unit variance initially
    )
  }
  
  # ========================================================================
  # 4. EM iterations
  # ========================================================================
  
  log_likelihood <- numeric(max_iter)
  converged <- FALSE
  
  for (iter in 1:max_iter) {
    
    # ----------------------------------------------------------------------
    # E-STEP: Infer posterior over sources
    # ----------------------------------------------------------------------
    
    # Compute current source estimates (MAP estimate)
    W <- MASS::ginv(A)  # Pseudo-inverse for unmixing
    S_hat <- W %*% X_centered
    
    # For each source, compute responsibilities for mixture components
    gamma <- array(0, dim = c(n_components, n_mixtures, n))
    
    for (k in 1:n_components) {
      params <- source_params[[k]]
      
      # Compute log-likelihood for each mixture component
      log_prob <- matrix(0, n_mixtures, n)
      for (m in 1:n_mixtures) {
        log_prob[m, ] <- log(params$alpha[m]) + 
          log_gaussian(S_hat[k, ], params$mu[m], params$sigma2[m])
      }
      
      # Normalize to get responsibilities (posterior probabilities)
      log_sum <- apply(log_prob, 2, function(x) {
        max_log <- max(x)
        max_log + log(sum(exp(x - max_log)))
      })
      
      for (m in 1:n_mixtures) {
        gamma[k, m, ] <- exp(log_prob[m, ] - log_sum)
      }
    }
    
    # Compute expected sufficient statistics
    E_S <- S_hat
    E_SS <- matrix(0, n_components, n_components)
    
    for (i in 1:n) {
      # For simplicity, use MAP estimate for second moment
      # (Full EM would integrate over posterior)
      E_SS <- E_SS + S_hat[, i, drop = FALSE] %*% t(S_hat[, i, drop = FALSE])
    }
    E_SS <- E_SS / n
    
    # ----------------------------------------------------------------------
    # M-STEP: Update parameters
    # ----------------------------------------------------------------------
    
    # (a) Update source model parameters (mixture components)
    for (k in 1:n_components) {
      for (m in 1:n_mixtures) {
        gamma_km <- gamma[k, m, ]
        N_km <- sum(gamma_km)
        
        if (N_km > 1e-10) {
          # Update mixture weight
          source_params[[k]]$alpha[m] <- N_km / n
          
          # Update mean
          source_params[[k]]$mu[m] <- sum(gamma_km * S_hat[k, ]) / N_km
          
          # Update variance
          source_params[[k]]$sigma2[m] <- 
            sum(gamma_km * (S_hat[k, ] - source_params[[k]]$mu[m])^2) / N_km
          
          # Prevent numerical issues
          source_params[[k]]$sigma2[m] <- 
            max(source_params[[k]]$sigma2[m], 1e-6)
        }
      }
      
      # Normalize mixture weights
      source_params[[k]]$alpha <- 
        source_params[[k]]$alpha / sum(source_params[[k]]$alpha)
    }
    
    # (b) Update mixing matrix A
    # Using weighted least squares: A = X * E[S]^T * (E[SS^T])^-1
    A_new <- X_centered %*% t(E_S) %*% solve(E_SS + diag(1e-6, n_components))
    
    # Normalize columns to prevent scaling ambiguity
    for (k in 1:n_components) {
      norm_k <- sqrt(sum(A_new[, k]^2))
      if (norm_k > 1e-10) {
        A_new[, k] <- A_new[, k] / norm_k
      }
    }
    
    A <- A_new
    
    # ----------------------------------------------------------------------
    # Compute log-likelihood
    # ----------------------------------------------------------------------
    
    # Reconstruction error term
    reconstruction <- X_centered - A %*% S_hat
    ll_obs <- -0.5 * n * d * log(2 * pi * noise_variance) - 
              0.5 * sum(reconstruction^2) / noise_variance
    
    # Source prior term (mixture model)
    ll_source <- 0
    for (k in 1:n_components) {
      params <- source_params[[k]]
      for (i in 1:n) {
        ll_comp <- numeric(n_mixtures)
        for (m in 1:n_mixtures) {
          ll_comp[m] <- log(params$alpha[m]) + 
            log_gaussian(S_hat[k, i], params$mu[m], params$sigma2[m])
        }
        # Log-sum-exp trick
        max_ll <- max(ll_comp)
        ll_source <- ll_source + max_ll + log(sum(exp(ll_comp - max_ll)))
      }
    }
    
    log_likelihood[iter] <- ll_obs + ll_source
    
    # ----------------------------------------------------------------------
    # Check convergence
    # ----------------------------------------------------------------------
    
    if (iter > 1) {
      delta_ll <- abs(log_likelihood[iter] - log_likelihood[iter - 1])
      rel_change <- delta_ll / (abs(log_likelihood[iter]) + 1e-10)
      
      if (verbose && (iter %% 10 == 0 || iter <= 5)) {
        cat(sprintf("Iter %3d: Log-Likelihood = %12.4f (Δ = %.2e)\n", 
                    iter, log_likelihood[iter], delta_ll))
      }
      
      if (rel_change < tol) {
        if (verbose) {
          cat(sprintf("\nConverged at iteration %d\n", iter))
        }
        converged <- TRUE
        log_likelihood <- log_likelihood[1:iter]
        break
      }
    } else {
      if (verbose) {
        cat(sprintf("Iter %3d: Log-Likelihood = %12.4f\n", 
                    iter, log_likelihood[iter]))
      }
    }
  }
  
  if (!converged && verbose) {
    cat(sprintf("\nMaximum iterations (%d) reached without convergence\n", max_iter))
  }
  
  # ========================================================================
  # 5. Final source estimation
  # ========================================================================
  
  W <- MASS::ginv(A)
  S <- W %*% X_centered
  
  # ========================================================================
  # 6. Return results
  # ========================================================================
  
  result <- list(
    A = A,                    # Mixing matrix
    S = S,                    # Estimated sources
    W = W,                    # Unmixing matrix
    source_params = source_params,  # Mixture model parameters
    log_likelihood = log_likelihood,
    converged = converged,
    n_iter = if (converged) iter else max_iter,
    X_mean = X_mean           # For reconstruction
  )
  
  class(result) <- "ica_em"
  return(result)
}

# ============================================================================
# Initialization Functions
# ============================================================================

#' Initialize mixing matrix
#'
#' @param X Data matrix (centered)
#' @param n_components Number of components
#' @param method Initialization method
#' @return Initial mixing matrix
initialize_mixing_matrix <- function(X, n_components, method = "pca") {
  
  d <- nrow(X)
  n <- ncol(X)
  
  if (method == "random") {
    # Random orthogonal initialization
    A <- matrix(rnorm(d * n_components), d, n_components)
    A <- qr.Q(qr(A))
    
  } else if (method == "pca") {
    # PCA-based initialization
    if (n_components == d) {
      A <- diag(d)
    } else {
      # Use top principal components
      svd_result <- svd(X)
      A <- svd_result$u[, 1:n_components, drop = FALSE]
      # Add small random rotation
      Q <- qr.Q(qr(matrix(rnorm(n_components * n_components), 
                          n_components, n_components)))
      A <- A %*% Q
    }
    
  } else if (method == "fastica") {
    # Simple FastICA-style initialization (requires additional implementation)
    warning("FastICA initialization not fully implemented, using PCA instead")
    A <- initialize_mixing_matrix(X, n_components, method = "pca")
    
  } else {
    stop("Unknown initialization method")
  }
  
  return(A)
}

# ============================================================================
# Utility Functions
# ============================================================================

#' Print method for ica_em objects
#' @export
print.ica_em <- function(x, ...) {
  cat("ICA-EM Results\n")
  cat("==============\n")
  cat(sprintf("Mixing matrix A: %d x %d\n", nrow(x$A), ncol(x$A)))
  cat(sprintf("Sources S: %d x %d\n", nrow(x$S), ncol(x$S)))
  cat(sprintf("Converged: %s\n", x$converged))
  cat(sprintf("Iterations: %d\n", x$n_iter))
  cat(sprintf("Final log-likelihood: %.4f\n", tail(x$log_likelihood, 1)))
}

#' Plot method for ica_em objects
#' @export
plot.ica_em <- function(x, type = "convergence", ...) {
  
  if (type == "convergence") {
    # Plot log-likelihood convergence
    plot(x$log_likelihood, type = "l", lwd = 2, col = "blue",
         xlab = "Iteration", ylab = "Log-Likelihood",
         main = "ICA-EM Convergence")
    grid()
    
  } else if (type == "sources") {
    # Plot estimated sources
    n_comp <- nrow(x$S)
    par(mfrow = c(ceiling(n_comp / 2), 2))
    for (k in 1:n_comp) {
      plot(x$S[k, ], type = "l", 
           main = sprintf("Source %d", k),
           ylab = "Amplitude", xlab = "Sample")
      grid()
    }
    par(mfrow = c(1, 1))
    
  } else if (type == "mixing") {
    # Visualize mixing matrix
    image(t(x$A[nrow(x$A):1, ]), 
          main = "Mixing Matrix A",
          xlab = "Component", ylab = "Dimension",
          col = heat.colors(20))
    
  } else {
    stop("Unknown plot type. Use 'convergence', 'sources', or 'mixing'")
  }
}

#' Reconstruct original data from ICA-EM results
#'
#' @param ica_result Object returned by ica_em()
#' @return Reconstructed data matrix
#' @export
reconstruct <- function(ica_result) {
  X_reconstructed <- ica_result$A %*% ica_result$S
  X_reconstructed <- X_reconstructed + ica_result$X_mean
  return(X_reconstructed)
}

# ============================================================================
# Example Usage
# ============================================================================

if (FALSE) {
  # This block demonstrates usage but won't run automatically
  
  # Generate synthetic data
  set.seed(42)
  n <- 1000
  t <- seq(0, 10, length.out = n)
  
  # True independent sources
  S_true <- rbind(
    sin(2 * pi * t),                    # Sinusoidal
    sign(sin(3 * pi * t)),              # Square wave
    runif(n, -1, 1)                     # Uniform noise
  )
  
  # Random mixing matrix
  A_true <- matrix(rnorm(3 * 3), 3, 3)
  
  # Mixed observations
  X <- A_true %*% S_true
  X <- X + matrix(rnorm(3 * n, sd = 0.1), 3, n)  # Add noise
  
  # Run ICA-EM
  result <- ica_em(X, n_components = 3, n_mixtures = 3, 
                   max_iter = 100, verbose = TRUE)
  
  # Print results
  print(result)
  
  # Plot convergence
  plot(result, type = "convergence")
  
  # Plot estimated sources
  plot(result, type = "sources")
  
  # Compare with true sources (need to handle permutation and scaling)
  par(mfrow = c(3, 2))
  for (i in 1:3) {
    plot(S_true[i, ], type = "l", main = sprintf("True Source %d", i))
    plot(result$S[i, ], type = "l", main = sprintf("Estimated Source %d", i))
  }
  par(mfrow = c(1, 1))
  
  # Check reconstruction error
  X_recon <- reconstruct(result)
  mse <- mean((X - X_recon)^2)
  cat(sprintf("\nReconstruction MSE: %.6f\n", mse))
}
