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
[`results/`](./results/). Installs the study package and starts the web app.

## Run the code

1. Install R package dependencies

    ```
    source("setup.R")
    ```

1. Perform the differential expression analysis. This reads the input files in
`data/` and exports the output files to `results/`

    ```
    library(rmarkdown)
    render("analyze.R")
    ```

1. Create and install the OmicNavigator study package. This reads the analysis
results files in `results/`, converts them to an OmicNavigator study package,
installs the package, and starts the app.

    ```
    source("build.R")
    ```

## Acknowledgements

The example limma+voom code was adapted from the Bioconductor workflow
[RNAseq123](https://bioconductor.org/packages/release/workflows/html/RNAseq123.html).

If you use the code, please cite:

> Law CW, Alhamdoosh M, Su S et al. RNA-seq analysis is easy as 1-2-3
> with limma, Glimma and edgeR [version 2; referees: 3 approved].
> F1000Research 2016, 5:1408 (doi: 10.12688/f1000research.9005.2)

If you use the data please cite:

> Sheridan JM, Ritchie ME, Best SA, et al.: A pooled shRNA screen for
> regulators of primary mammary stem and progenitor cells identifies
> roles for Asap1 and Prox1. BMC Cancer. 2015; 15(1): 221.
