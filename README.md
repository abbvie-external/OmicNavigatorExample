# How to create an OmicNavigator study

This repository contains an example of how to convert an omic analysis into a
study package to be explored with the OmicNavigator app.

**Files:**

* [`setup.R`](./setup.R) - Installs the required R packages for the example

* [`analyze.R`](./analyze.R) - Performs a differential expression and enrichment
analysis of an RNA-seq experiment. It uses limma+voom, but OmicNavigator is
agnostic to how you perform your analysis. If you prefer, you could use Python
or a GUI to perform the analysis, as long as you export the results.

* [`build.R`](./build.R) - Builds the OmicNavigator study package from the
results of `analyze.R`
