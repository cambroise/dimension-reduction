# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

Quarto book project for a lecture on **Dimension Reduction and Visualization**. The same material serves two purposes:
- **Book/notes** format: rendered as HTML (web) and PDF (print)
- **Slides** format: RevealJS presentations in `slides/`

## Build Commands

```bash
# Preview with live reload (HTML only, fast iteration)
quarto preview

# Full render (all formats: HTML + PDF)
quarto render

# Render a single file
quarto render chapters/svd.qmd

# Render with solutions visible
quarto render -P solution:true

# Render slides only
quarto render slides/01-svd.qmd
```

Output goes to `_book/`. Python 3.11.8 (conda env: `machine-learning`) is the default kernel.

## Content Architecture

Each topic follows a **chapter → exercises** pattern in `_quarto.yml`:

```
chapters/svd.qmd → exercises/exercices-svd.qmd
chapters/factor-analysis.qmd → exercises/exercices-fa.qmd
chapters/ICA.qmd → exercises/exercices-ica.qmd
chapters/neural-unsupervised.qmd → exercises/exercices-neural.qmd
chapters/vae.qmd → exercises/exercices-vae.qmd
chapters/mds.qmd → exercises/exercices-mds.qmd
chapters/t-SNE.qmd (no exercises yet)
```

Slides in `slides/` mirror chapter content but use RevealJS format.

## Key Conventions

**Solution visibility** — exercise solutions are hidden by default and toggled via a global parameter:
```yaml
# In _quarto.yml
params:
  solution: false
```
Solutions are wrapped with:
```markdown
::: {.content-hidden unless-meta="params.solution"}
::: {.callout-tip title="Solution" collapse="true"}
...
:::
:::
```

**Math** — uses standard LaTeX (`$...$` inline, `$$...$$` display).

**Citations** — `@key` style, bibliography in `references.bib`.

**Slide-specific content** — use `.allowframebreaks` class for long slides.

## Formats

- **HTML**: cosmo theme + custom brand styling (`brand.scss`)
- **PDF**: KOMA-Script `scrreprt` document class
- **Slides**: RevealJS (in `slides/`) + compiled to PDF

## Data and Assets

- `data/` — MNIST dataset, `Peppers.tiff`, `shakespeare-texts.txt`
- `Silo-Images/` — images used across chapters
- `documents/` — local copies of key academic papers
- `projects/` — student project descriptions + reference R implementations (ICA-EM, ICA-SGD)
