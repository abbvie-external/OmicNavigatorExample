# Installs the required packages from CRAN, Bioconductor, and GitHub.

cran <- c("BiocManager", "ggplot2", "remotes")
install.packages(cran)

bioc <- c("edgeR", "limma", "Mus.musculus")
BiocManager::install(bioc, update = FALSE)

github <- c("abbvie-external/OmicNavigator@*release")
remotes::install_github(github, dependencies = TRUE, upgrade = FALSE)
