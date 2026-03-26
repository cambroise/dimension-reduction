# ============================================================================
# ICA-EM Demonstration Script
# ============================================================================
# This script demonstrates various use cases of the ICA-EM implementation
# ============================================================================

source("ica_em.R")

# ============================================================================
# Example 1: Simple Synthetic Data (3 sources)
# ============================================================================

cat("Example 1: Simple Synthetic Data\n")
cat("=================================\n\n")

set.seed(123)
n <- 1000
t <- seq(0, 10, length.out = n)

# Create three distinct source signals
S_true <- rbind(
  sin(2 * pi * 0.5 * t),              # Slow sinusoid
  sign(sin(2 * pi * 1.5 * t)),        # Square wave
  rt(n, df = 3) / 3                   # Heavy-tailed (Student-t)
)

# Random mixing
A_true <- matrix(c(
  0.8, 0.5, 0.3,
  0.2, 0.9, 0.4,
  0.5, 0.2, 0.8
), 3, 3, byrow = TRUE)

# Generate observations
X <- A_true %*% S_true
X <- X + matrix(rnorm(3 * n, sd = 0.05), 3, n)  # Small noise

# Run ICA-EM
result1 <- ica_em(X, n_components = 3, n_mixtures = 4, 
                  max_iter = 50, verbose = TRUE)

# Visualize results
png("example1_convergence.png", width = 800, height = 600)
plot(result1, type = "convergence")
dev.off()

png("example1_sources.png", width = 1000, height = 800)
plot(result1, type = "sources")
dev.off()

cat("\nMixing Matrix A (estimated):\n")
print(round(result1$A, 3))

cat("\nMixing Matrix A (true):\n")
print(round(A_true, 3))

# ============================================================================
# Example 2: Audio Signal Separation (Simulated)
# ============================================================================

cat("\n\nExample 2: Audio Signal Separation\n")
cat("===================================\n\n")

set.seed(456)
n <- 2000
t <- seq(0, 4, length.out = n)

# Simulate different "audio" sources
S_audio <- rbind(
  # Voice-like: modulated sinusoid
  sin(2 * pi * 200 * t) * (1 + 0.5 * sin(2 * pi * 5 * t)),
  
  # Music-like: sum of harmonics
  sin(2 * pi * 440 * t) + 0.5 * sin(2 * pi * 880 * t) + 0.25 * sin(2 * pi * 1320 * t),
  
  # Noise-like: filtered random
  filter(rnorm(n), rep(1/10, 10), sides = 1)
)

# Normalize sources
for (i in 1:3) {
  S_audio[i, ] <- S_audio[i, ] / sd(S_audio[i, ])
}

# Random mixing
A_audio <- matrix(rnorm(3 * 3), 3, 3)
A_audio <- A_audio / sqrt(rowSums(A_audio^2))

# Generate mixed signals
X_audio <- A_audio %*% S_audio
X_audio <- X_audio + matrix(rnorm(3 * n, sd = 0.1), 3, n)

# Run ICA-EM with more mixture components for complex distributions
result2 <- ica_em(X_audio, n_components = 3, n_mixtures = 5,
                  max_iter = 100, tol = 1e-5, verbose = TRUE)

# Plot results
png("example2_convergence.png", width = 800, height = 600)
plot(result2, type = "convergence")
dev.off()

cat("\nReconstruction error: ", 
    mean((X_audio - reconstruct(result2))^2), "\n")

# ============================================================================
# Example 3: Dimensionality Reduction (More observations than components)
# ============================================================================

cat("\n\nExample 3: Dimensionality Reduction\n")
cat("====================================\n\n")

set.seed(789)
n <- 500
d <- 5  # 5 observed dimensions
k <- 2  # 2 underlying sources

# Two true sources
S_latent <- rbind(
  rnorm(n),
  rexp(n, rate = 1) - 1  # Non-Gaussian
)

# Random projection to higher dimension
A_proj <- matrix(rnorm(d * k), d, k)

# Generate high-dimensional observations
X_hd <- A_proj %*% S_latent
X_hd <- X_hd + matrix(rnorm(d * n, sd = 0.2), d, n)

# Run ICA-EM to recover the 2 sources
result3 <- ica_em(X_hd, n_components = 2, n_mixtures = 3,
                  max_iter = 80, verbose = TRUE)

cat("\nEstimated mixing matrix (5x2):\n")
print(round(result3$A, 3))

png("example3_sources.png", width = 1000, height = 400)
par(mfrow = c(1, 2))
plot(result3$S[1, ], type = "l", main = "Estimated Source 1",
     ylab = "Amplitude", xlab = "Sample")
plot(result3$S[2, ], type = "l", main = "Estimated Source 2",
     ylab = "Amplitude", xlab = "Sample")
dev.off()

# ============================================================================
# Example 4: Comparison of Initialization Methods
# ============================================================================

cat("\n\nExample 4: Initialization Comparison\n")
cat("=====================================\n\n")

set.seed(999)
n <- 800
t <- seq(0, 8, length.out = n)

S_comp <- rbind(
  sin(2 * pi * t),
  (t %% 2) - 1,  # Sawtooth
  rnorm(n)
)

A_comp <- matrix(c(1, 0.5, 0.3, 0.4, 1, 0.2, 0.3, 0.4, 1), 3, 3)
X_comp <- A_comp %*% S_comp

# Test different initializations
init_methods <- c("random", "pca")
results_init <- list()

for (method in init_methods) {
  cat(sprintf("\nTesting %s initialization...\n", method))
  results_init[[method]] <- ica_em(X_comp, n_components = 3, 
                                   n_mixtures = 3, max_iter = 50,
                                   init_method = method, verbose = FALSE)
}

# Compare convergence
png("example4_init_comparison.png", width = 1000, height = 600)
plot(results_init[["pca"]]$log_likelihood, type = "l", lwd = 2, col = "blue",
     xlab = "Iteration", ylab = "Log-Likelihood",
     main = "Convergence Comparison: Initialization Methods",
     ylim = range(c(results_init[["pca"]]$log_likelihood,
                    results_init[["random"]]$log_likelihood)))
lines(results_init[["random"]]$log_likelihood, lwd = 2, col = "red")
legend("bottomright", legend = c("PCA", "Random"), 
       col = c("blue", "red"), lwd = 2)
grid()
dev.off()

cat("\nFinal log-likelihoods:\n")
cat(sprintf("  PCA:    %.4f\n", tail(results_init[["pca"]]$log_likelihood, 1)))
cat(sprintf("  Random: %.4f\n", tail(results_init[["random"]]$log_likelihood, 1)))

# ============================================================================
# Example 5: Robustness to Noise
# ============================================================================

cat("\n\nExample 5: Robustness to Noise\n")
cat("================================\n\n")

set.seed(111)
n <- 600
t <- seq(0, 6, length.out = n)

S_clean <- rbind(
  sin(2 * pi * 0.8 * t),
  sign(sin(2 * pi * 1.2 * t))
)

A_noise <- matrix(c(1, 0.5, 0.5, 1), 2, 2)

# Test different noise levels
noise_levels <- c(0.01, 0.1, 0.5)
results_noise <- list()

for (noise_sd in noise_levels) {
  X_noisy <- A_noise %*% S_clean + 
             matrix(rnorm(2 * n, sd = noise_sd), 2, n)
  
  cat(sprintf("\nNoise level: %.2f\n", noise_sd))
  results_noise[[as.character(noise_sd)]] <- 
    ica_em(X_noisy, n_components = 2, n_mixtures = 3,
           noise_variance = noise_sd^2, max_iter = 50, verbose = FALSE)
}

# Compare reconstruction errors
cat("\nReconstruction Errors:\n")
for (noise_sd in noise_levels) {
  result <- results_noise[[as.character(noise_sd)]]
  X_noisy <- A_noise %*% S_clean + 
             matrix(rnorm(2 * n, sd = noise_sd), 2, n)
  mse <- mean((X_noisy - reconstruct(result))^2)
  cat(sprintf("  Noise %.2f: MSE = %.6f\n", noise_sd, mse))
}

# ============================================================================
# Summary
# ============================================================================

cat("\n\n============================================\n")
cat("ICA-EM Demonstration Complete\n")
cat("============================================\n")
cat("Generated plots:\n")
cat("  - example1_convergence.png\n")
cat("  - example1_sources.png\n")
cat("  - example2_convergence.png\n")
cat("  - example3_sources.png\n")
cat("  - example4_init_comparison.png\n")
cat("\nAll examples completed successfully!\n")
