# ============================================================================
# ICA via EM with SGD-approximated E-step
# ============================================================================
# This implementation uses EM where:
# - E-step: Approximate posterior p(S|X,A) using SGD on the negative log-posterior
# - M-step: Update mixing matrix A using the posterior statistics from E-step
# ============================================================================

library(MASS)  # For matrix operations

# ============================================================================
# Helper Functions
# ============================================================================

#' Compute log-prior for a source component (mixture of Gaussians)
#'
#' @param s Source values
#' @param params List with alpha (weights), mu (means), sigma2 (variances)
#' @return Log-prior for each observation
log_prior_mog <- function(s, params) {
  n_mixtures <- length(params$alpha)
  n <- length(s)
  
  log_prob <- matrix(0, n_mixtures, n)
  for (m in 1:n_mixtures) {
    log_prob[m, ] <- log(params$alpha[m]) - 
      0.5 * log(2 * pi * params$sigma2[m]) - 
      0.5 * (s - params$mu[m])^2 / params$sigma2[m]
  }
  
  # Log-sum-exp for numerical stability
  max_log <- apply(log_prob, 2, max)
  log_sum <- max_log + log(colSums(exp(log_prob - rep(max_log, each = n_mixtures))))
  
  return(log_sum)
}

#' Compute negative log-posterior for sources (objective for SGD E-step)
#'
#' @param S Sources (k x n)
#' @param X Observations (d x n) 
#' @param A Mixing matrix (d x k)
#' @param source_params List of source prior parameters
#' @param noise_var Observation noise variance
#' @return Negative log-posterior
neg_log_posterior <- function(S, X, A, source_params, noise_var) {
  k <- nrow(S)
  n <- ncol(S)
  d <- nrow(X)
  
  # Log-likelihood: log p(X|S,A)
  # X = AS + noise, so log p(X|S,A) = -1/(2*sigma^2) * ||X - AS||^2 + const
  reconstruction <- X - A %*% S
  log_likelihood <- -0.5 * sum(reconstruction^2) / noise_var
  
  # Log-prior: log p(S) = sum_i log p(s_i)
  log_prior <- 0
  for (i in 1:k) {
    log_prior <- log_prior + sum(log_prior_mog(S[i, ], source_params[[i]]))
  }
  
  # Negative log-posterior (we minimize this in E-step)
  return(-(log_likelihood + log_prior))
}

#' Gradient of negative log-posterior w.r.t. S
#'
#' @param S Sources (k x n)
#' @param X Observations (d x n)
#' @param A Mixing matrix (d x k)
#' @param source_params List of source prior parameters
#' @param noise_var Observation noise variance
#' @return Gradient matrix (k x n)
grad_neg_log_posterior <- function(S, X, A, source_params, noise_var) {
  k <- nrow(S)
  n <- ncol(S)
  
  # Gradient from likelihood: d/dS[-1/(2*sigma^2) * ||X - AS||^2]
  #                         = 1/sigma^2 * A^T(AS - X)
  grad_likelihood <- (t(A) %*% (A %*% S - X)) / noise_var
  
  # Gradient from prior: d/dS[log p(S)]
  # For mixture of Gaussians, this requires computing responsibilities
  grad_prior <- matrix(0, k, n)
  
  for (i in 1:k) {
    s_i <- S[i, ]
    params <- source_params[[i]]
    n_mixtures <- length(params$alpha)
    
    # Compute responsibilities gamma[m,j] = p(component m | s_j)
    log_prob <- matrix(0, n_mixtures, n)
    for (m in 1:n_mixtures) {
      log_prob[m, ] <- log(params$alpha[m]) - 
        0.5 * log(2 * pi * params$sigma2[m]) - 
        0.5 * (s_i - params$mu[m])^2 / params$sigma2[m]
    }
    
    # Normalize to get responsibilities
    max_log <- apply(log_prob, 2, max)
    log_sum <- max_log + log(colSums(exp(log_prob - rep(max_log, each = n_mixtures))))
    gamma <- exp(log_prob - rep(log_sum, each = n_mixtures))
    
    # Gradient: d/ds_i[log p(s_i)] = sum_m gamma_m * (mu_m - s_i) / sigma2_m
    for (m in 1:n_mixtures) {
      grad_prior[i, ] <- grad_prior[i, ] + 
        gamma[m, ] * (params$mu[m] - s_i) / params$sigma2[m]
    }
  }
  
  # Combine (remember we're computing gradient of NEGATIVE log-posterior)
  return(grad_likelihood - grad_prior)
}

# ============================================================================
# Main ICA-EM with SGD E-step Algorithm
# ============================================================================

#' ICA via EM with SGD-approximated E-step
#'
#' @param X Data matrix (d x n): d = dimensions, n = observations
#' @param n_components Number of independent components to extract
#' @param n_mixtures Number of Gaussian mixture components per source
#' @param em_max_iter Maximum number of EM iterations
#' @param sgd_max_iter Maximum number of SGD iterations in E-step
#' @param sgd_batch_size Mini-batch size for SGD in E-step
#' @param sgd_lr Learning rate for SGD in E-step
#' @param noise_variance Observation noise variance
#' @param tol Convergence tolerance
#' @param verbose Print iteration details
#'
#' @return List containing:
#'   - A: Estimated mixing matrix (d x n_components)
#'   - S: Estimated sources (n_components x n)
#'   - W: Unmixing matrix
#'   - log_likelihood: Log-likelihood at each EM iteration
#'   - converged: Whether algorithm converged
#'
#' @export
ica_em_sgd_estep <- function(X, 
                              n_components = NULL,
                              n_mixtures = 3,
                              em_max_iter = 50,
                              sgd_max_iter = 100,
                              sgd_batch_size = 32,
                              sgd_lr = 0.01,
                              noise_variance = 0.1,
                              tol = 1e-4,
                              verbose = TRUE) {
  
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
    cat(sprintf("ICA-EM with SGD E-step\n"))
    cat(sprintf("======================\n"))
    cat(sprintf("Data dimensions: %d x %d\n", d, n))
    cat(sprintf("Components: %d\n", n_components))
    cat(sprintf("Mixture components per source: %d\n", n_mixtures))
    cat(sprintf("EM iterations: %d, SGD iterations per E-step: %d\n", 
                em_max_iter, sgd_max_iter))
    cat(sprintf("SGD batch size: %d, learning rate: %.4f\n\n", 
                sgd_batch_size, sgd_lr))
  }
  
  # ========================================================================
  # 2. Initialize mixing matrix A
  # ========================================================================
  
  # PCA-based initialization
  if (n_components == d) {
    A <- diag(d)
  } else {
    svd_result <- svd(X_centered)
    A <- svd_result$u[, 1:n_components, drop = FALSE]
    # Add small random rotation
    Q <- qr.Q(qr(matrix(rnorm(n_components * n_components), 
                        n_components, n_components)))
    A <- A %*% Q
  }
  
  # ========================================================================
  # 3. Initialize source model parameters (mixture of Gaussians)
  # ========================================================================
  
  source_params <- list()
  for (k in 1:n_components) {
    source_params[[k]] <- list(
      alpha = rep(1/n_mixtures, n_mixtures),
      mu = seq(-2, 2, length.out = n_mixtures),
      sigma2 = rep(1, n_mixtures)
    )
  }
  
  # Initialize sources (initial E-step)
  W <- MASS::ginv(A)
  S <- W %*% X_centered
  
  # ========================================================================
  # 4. EM iterations
  # ========================================================================
  
  log_likelihood <- numeric(em_max_iter)
  converged <- FALSE
  
  for (em_iter in 1:em_max_iter) {
    
    # ----------------------------------------------------------------------
    # E-STEP: Approximate MAP estimate of S using SGD
    # ----------------------------------------------------------------------
    
    # We want to find S that maximizes p(S|X,A) = p(X|S,A) * p(S) / p(X)
    # Equivalently, minimize negative log-posterior: -log p(X|S,A) - log p(S)
    
    S_current <- S  # Start from previous iteration's S
    
    # SGD to optimize S
    for (sgd_iter in 1:sgd_max_iter) {
      # Sample mini-batch
      if (sgd_batch_size >= n) {
        # Use full batch
        batch_idx <- 1:n
      } else {
        batch_idx <- sample(n, sgd_batch_size, replace = FALSE)
      }
      
      X_batch <- X_centered[, batch_idx, drop = FALSE]
      S_batch <- S_current[, batch_idx, drop = FALSE]
      
      # Compute gradient on mini-batch
      grad <- grad_neg_log_posterior(S_batch, X_batch, A, 
                                     source_params, noise_variance)
      
      # SGD update (only update the batch)
      S_current[, batch_idx] <- S_batch - sgd_lr * grad
      
      # Optional: gradient clipping for stability
      S_current[, batch_idx] <- pmax(pmin(S_current[, batch_idx], 10), -10)
    }
    
    S <- S_current
    
    # Compute posterior statistics for M-step
    # For simplicity, we use the MAP estimate S directly
    E_S <- S
    E_SS <- (S %*% t(S)) / n
    
    # Also compute responsibilities for mixture component updates
    gamma <- array(0, dim = c(n_components, n_mixtures, n))
    
    for (k in 1:n_components) {
      params <- source_params[[k]]
      
      log_prob <- matrix(0, n_mixtures, n)
      for (m in 1:n_mixtures) {
        log_prob[m, ] <- log(params$alpha[m]) - 
          0.5 * log(2 * pi * params$sigma2[m]) - 
          0.5 * (S[k, ] - params$mu[m])^2 / params$sigma2[m]
      }
      
      # Normalize
      max_log <- apply(log_prob, 2, max)
      log_sum <- max_log + log(colSums(exp(log_prob - rep(max_log, each = n_mixtures))))
      
      for (m in 1:n_mixtures) {
        gamma[k, m, ] <- exp(log_prob[m, ] - log_sum)
      }
    }
    
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
          source_params[[k]]$mu[m] <- sum(gamma_km * S[k, ]) / N_km
          
          # Update variance
          source_params[[k]]$sigma2[m] <- 
            sum(gamma_km * (S[k, ] - source_params[[k]]$mu[m])^2) / N_km
          
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
    
    # Normalize columns
    for (k in 1:n_components) {
      norm_k <- sqrt(sum(A_new[, k]^2))
      if (norm_k > 1e-10) {
        A_new[, k] <- A_new[, k] / norm_k
      }
    }
    
    A <- A_new
    
    # ----------------------------------------------------------------------
    # Compute log-likelihood (for monitoring)
    # ----------------------------------------------------------------------
    
    reconstruction <- X_centered - A %*% S
    ll_obs <- -0.5 * sum(reconstruction^2) / noise_variance
    
    ll_source <- 0
    for (k in 1:n_components) {
      ll_source <- ll_source + sum(log_prior_mog(S[k, ], source_params[[k]]))
    }
    
    log_likelihood[em_iter] <- ll_obs + ll_source
    
    # ----------------------------------------------------------------------
    # Check convergence
    # ----------------------------------------------------------------------
    
    if (em_iter > 1) {
      delta_ll <- abs(log_likelihood[em_iter] - log_likelihood[em_iter - 1])
      rel_change <- delta_ll / (abs(log_likelihood[em_iter]) + 1e-10)
      
      if (verbose && (em_iter %% 5 == 0 || em_iter <= 3)) {
        cat(sprintf("EM Iter %3d: Log-Likelihood = %12.4f (Δ = %.2e)\n", 
                    em_iter, log_likelihood[em_iter], delta_ll))
      }
      
      if (rel_change < tol) {
        if (verbose) {
          cat(sprintf("\nConverged at EM iteration %d\n", em_iter))
        }
        converged <- TRUE
        log_likelihood <- log_likelihood[1:em_iter]
        break
      }
    } else {
      if (verbose) {
        cat(sprintf("EM Iter %3d: Log-Likelihood = %12.4f\n", 
                    em_iter, log_likelihood[em_iter]))
      }
    }
  }
  
  if (!converged && verbose) {
    cat(sprintf("\nMaximum EM iterations (%d) reached without convergence\n", 
                em_max_iter))
  }
  
  # ========================================================================
  # 5. Final estimates
  # ========================================================================
  
  W <- MASS::ginv(A)
  
  # ========================================================================
  # 6. Return results
  # ========================================================================
  
  result <- list(
    A = A,
    S = S,
    W = W,
    source_params = source_params,
    log_likelihood = log_likelihood,
    converged = converged,
    n_em_iter = if (converged) em_iter else em_max_iter,
    X_mean = X_mean
  )
  
  class(result) <- "ica_em_sgd"
  return(result)
}

# ============================================================================
# Utility Functions
# ============================================================================

#' Print method for ica_em_sgd objects
#' @export
print.ica_em_sgd <- function(x, ...) {
  cat("ICA-EM with SGD E-step Results\n")
  cat("==============================\n")
  cat(sprintf("Mixing matrix A: %d x %d\n", nrow(x$A), ncol(x$A)))
  cat(sprintf("Sources S: %d x %d\n", nrow(x$S), ncol(x$S)))
  cat(sprintf("Converged: %s\n", x$converged))
  cat(sprintf("EM Iterations: %d\n", x$n_em_iter))
  cat(sprintf("Final log-likelihood: %.4f\n", tail(x$log_likelihood, 1)))
}

#' Plot method for ica_em_sgd objects
#' @export
plot.ica_em_sgd <- function(x, type = "convergence", ...) {
  
  if (type == "convergence") {
    plot(x$log_likelihood, type = "l", lwd = 2, col = "blue",
         xlab = "EM Iteration", ylab = "Log-Likelihood",
         main = "ICA-EM with SGD E-step Convergence")
    grid()
    
  } else if (type == "sources") {
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
    image(t(x$A[nrow(x$A):1, ]), 
          main = "Mixing Matrix A",
          xlab = "Component", ylab = "Dimension",
          col = heat.colors(20))
    
  } else {
    stop("Unknown plot type. Use 'convergence', 'sources', or 'mixing'")
  }
}

#' Reconstruct original data
#'
#' @param ica_result Object returned by ica_em_sgd_estep()
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
  # This demonstrates usage but won't run automatically
  
  # Generate synthetic data
  set.seed(42)
  n <- 1000
  t <- seq(0, 10, length.out = n)
  
  # True independent sources
  S_true <- rbind(
    sin(2 * pi * 0.5 * t),              # Sinusoidal
    sign(sin(2 * pi * 1.5 * t)),        # Square wave
    rt(n, df = 3) / 3                   # Heavy-tailed
  )
  
  # Random mixing matrix
  A_true <- matrix(rnorm(3 * 3), 3, 3)
  
  # Mixed observations
  X <- A_true %*% S_true
  X <- X + matrix(rnorm(3 * n, sd = 0.1), 3, n)
  
  # Run ICA-EM with SGD E-step
  result <- ica_em_sgd_estep(
    X, 
    n_components = 3,
    n_mixtures = 3,
    em_max_iter = 30,
    sgd_max_iter = 50,
    sgd_batch_size = 64,
    sgd_lr = 0.01,
    noise_variance = 0.1,
    verbose = TRUE
  )
  
  # Print results
  print(result)
  
  # Plot convergence
  plot(result, type = "convergence")
  
  # Plot estimated sources
  plot(result, type = "sources")
  
  # Compare with true sources
  par(mfrow = c(3, 2))
  for (i in 1:3) {
    plot(S_true[i, ], type = "l", main = sprintf("True Source %d", i),
         ylab = "Amplitude")
    plot(result$S[i, ], type = "l", main = sprintf("Estimated Source %d", i),
         ylab = "Amplitude")
  }
  par(mfrow = c(1, 1))
  
  # Reconstruction error
  X_recon <- reconstruct(result)
  mse <- mean((X - X_recon)^2)
  cat(sprintf("\nReconstruction MSE: %.6f\n", mse))
}

# ============================================================================
# Additional: Comparison Function
# ============================================================================

#' Compare ICA-EM-SGD with standard ICA-EM
#'
#' @param X Data matrix
#' @param n_components Number of components
#' @export
compare_ica_methods <- function(X, n_components = NULL) {
  
  cat("Running ICA-EM with SGD E-step...\n")
  result_sgd <- ica_em_sgd_estep(
    X,
    n_components = n_components,
    em_max_iter = 30,
    sgd_max_iter = 50,
    verbose = FALSE
  )
  
  cat("\nRunning standard ICA-EM (if available)...\n")
  # Note: Would need to source the standard ica_em.R here
  # For now, just show the SGD version results
  
  cat("\n=== Results ===\n")
  cat(sprintf("ICA-EM-SGD: %d EM iterations, final LL = %.4f\n",
              result_sgd$n_em_iter,
              tail(result_sgd$log_likelihood, 1)))
  
  # Plot comparison
  par(mfrow = c(1, 2))
  plot(result_sgd$log_likelihood, type = "l", lwd = 2, col = "blue",
       main = "ICA-EM with SGD E-step",
       xlab = "EM Iteration", ylab = "Log-Likelihood")
  grid()
  
  par(mfrow = c(1, 1))
  
  return(result_sgd)
}
