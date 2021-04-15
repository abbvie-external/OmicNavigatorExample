# How to create an OmicNavigator study

This repository contains an example of how to convert an omic analysis into a
study package to be explored with the OmicNavigator app.

**Files:**

* [`setup.R`](./setup.R) - Installs the required R packages for the example

* [`data/`](./data/) - The input RNA-seq counts and MSigDB annotations

* [`analyze.R`](./analyze.R) - Performs a differential expression and enrichment
analysis of the RNA-seq experiment in [`data/`](./data/). It uses limma+voom,
but OmicNavigator is agnostic to how you perform your analysis. If you prefer,
you could use Python or a GUI to perform the analysis, as long as you export the
results.

* [`results/`](./results/) - The analysis results exported by
[`analyze.R`](./analyze.R)

* [`build.R`](./build.R) - Builds the OmicNavigator study package from the
results files that [`analyze.R`](./analyze.R) exported to
[`results/`](./results/)
