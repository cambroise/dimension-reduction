<img src="Silo-Images/logo.png" alt="Logo" width="200"/>

# Lecture Notes: Dimension Reduction and Visualization

> This repository is primarily a collection of lecture slides compiled into a Quarto book for convenience. It should not be considered a comprehensive textbook or self-contained learning resource.

A Quarto book covering the main methods for dimensionality reduction and visualization, from classical linear methods to modern deep learning approaches.

## Contents

| Chapter | Topic |
|---------|-------|
| SVD & PCA | Singular Value Decomposition, Principal Component Analysis |
| Factor Analysis | Latent factor models |
| ICA | Independent Component Analysis |
| Neural Unsupervised | Autoencoders, deep representations |
| VAE | Variational Autoencoders |
| MDS | Multidimensional Scaling |
| t-SNE | t-distributed Stochastic Neighbor Embedding |

## Build

```bash
# Preview (HTML, live reload)
quarto preview

# Full render (HTML + PDF)
quarto render

# Render with solutions
quarto render -P solution:true
```

Requires Python 3.11+ (conda env: `machine-learning`).

## Compiled Notes

The latest compiled PDF is available in [`_book/`](_book/).
