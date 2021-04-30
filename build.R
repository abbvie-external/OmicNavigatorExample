# Create OmicNavigator study from files in the directory results/

library(ggplot2)
library(OmicNavigator)

# Create a new study -----------------------------------------------------------

study <- createStudy("RNAseq123",
                     "Bioc workflow package converted to OmicNavigator",
                     version = "0.1.0")

# Models -----------------------------------------------------------------------

models <- list(
  main = "A standard group-means analysis of 3 mammary cell populations"
)
study <- addModels(study, models)

# Samples ----------------------------------------------------------------------

samples <- read.delim("results/samples.txt", stringsAsFactors = FALSE)
samples <- list(main = samples)
study <- addSamples(study, samples)

# Features ---------------------------------------------------------------------

features <- read.delim("results/features.txt", stringsAsFactors = FALSE)
# The featureIDs in the first column must be a character vector
features$entrez <- as.character(features$entrez)
features <- list(main = features)
study <- addFeatures(study, features)

# MetaFeatures -----------------------------------------------------------------

metaFeatures <- read.delim("results/metaFeatures.txt", stringsAsFactors = FALSE,
                           colClasses = c(entrez = "character"))
metaFeatures <- list(main = metaFeatures)
study <- addMetaFeatures(study, metaFeatures)

# Assays -----------------------------------------------------------------------

assays <- read.delim("results/assays.txt", stringsAsFactors = FALSE)
assays <- list(main = assays)
study <- addAssays(study, assays)

# Results (differential expression) --------------------------------------------

# The featureID in the first column must be a character vector
BasalvsLP <- read.delim("results/BasalvsLP.txt", stringsAsFactors = FALSE,
                        colClasses = c(entrez = "character"))
BasalvsML <- read.delim("results/BasalvsML.txt", stringsAsFactors = FALSE,
                        colClasses = c(entrez = "character"))
LPvsML <- read.delim("results/LPvsML.txt", stringsAsFactors = FALSE,
                     colClasses = c(entrez = "character"))
results <- list(
  main = list(
    BasalvsLP = LPvsML,
    BasalvsML = BasalvsML,
    LPvsML = LPvsML
  )
)
study <- addResults(study, results)

# Tests ------------------------------------------------------------------------

tests <- list(
  main = list(
    BasalvsLP = "Which genes are DE between Basal and LP cells?",
    BasalvsML = "Which genes are DE between Basal and ML cells?",
    LPvsML = "Which genes are DE between LP and ML cells?"
  )
)
study <- addTests(study, tests)

# Linkouts to external resources for the results table -------------------------

resultsLinkouts <- list(
  main = list(
    entrez = c("https://www.ncbi.nlm.nih.gov/gene/",
               "https://ensembl.org/Mus_musculus/Gene/Summary?g="),
    symbol = "http://www.informatics.jax.org/searchtool/Search.do?query="
  )
)
study <- addResultsLinkouts(study, resultsLinkouts)

# Annotations (used for enrichments) -------------------------------------------

load("data/mouse_H_v5p2.rdata")
annotations <- list(
  mouse_H_v5p2 = list(
    terms = Mm.H,
    description = "The gene signatures from the Broad Institute's MSigDB Hallmark collection",
    featureID = "entrez"
  )
)
study <- addAnnotations(study, annotations)

# Enrichments ------------------------------------------------------------------

enrichedBasalvsLP <- read.delim("results/enrichedBasalvsLP.txt",
                                stringsAsFactors = FALSE)
enrichedBasalvsML <- read.delim("results/enrichedBasalvsML.txt",
                                stringsAsFactors = FALSE)
enrichedLPvsML <- read.delim("results/enrichedLPvsML.txt",
                             stringsAsFactors = FALSE)

create_enrichments <- function(x) {
  data.frame(
    termID = row.names(x),
    description = gsub("_", " ", tolower(row.names(x))),
    nominal = x$PValue,
    adjusted = x$FDR,
    stringsAsFactors = FALSE
  )
}
enrichments <- list(
  main = list(
    mouse_H_v5p2 = list(
      BasalvsLP = create_enrichments(enrichedBasalvsLP),
      BasalvsML = create_enrichments(enrichedBasalvsML),
      LPvsML = create_enrichments(enrichedLPvsML)
    )
  )
)
study <- addEnrichments(study, enrichments)

# Linkouts to external resources for the enrichments table ---------------------

enrichmentsLinkouts <- list(
  mouse_H_v5p2 = "https://www.gsea-msigdb.org/gsea/msigdb/cards/"
)
study <- addEnrichmentsLinkouts(study, enrichmentsLinkouts)

# Barcodes (for barcode plot of features enriched in a given annotation) -------

barcodes <- list(
  default = list(
    statistic = "t", # Use the t-statistics from the results table
    absolute = TRUE, # Take the absolute value of the t-statistics
    logFoldChange = "logFC", # Use the column "logFC" from the results table for violin plot
    labelStat = "abs(t)",
    labelLow = "Small effect",
    labelHigh = "Large effect"
  )
)
study <- addBarcodes(study, barcodes = barcodes)

# Reports ----------------------------------------------------------------------

reports <- list(
  main = "results/report.html"
)
study <- addReports(study, reports)

# Custom plots -----------------------------------------------------------------

x <- getPlottingData(study, modelID = "main", featureID = "497097")

expression_by_cell_type <- function(x) {
  ggDataFrame <- cbind(x$samples,
                       feature = as.numeric(x$assays))
  ggplot(ggDataFrame, aes(x = .data$group, y = .data$feature, fill = .data$group)) +
    geom_boxplot() +
    labs(x = "Cell type", y = "Gene expression",
         title = sprintf("%s (Entrez %s)", x$features$symbol, x$features$entrez)) +
    scale_fill_manual("Cell type", values = c("pink", "purple", "gold")) +
    theme_classic(base_size = 24)
}
expression_by_cell_type(x)

plots <- list(
  main = list(
    expression_by_cell_type = list(
      displayName = "Expression by cell type",
      plotType = "singleFeature",
      packages = c("ggplot2")
    )
  )
)
study <- addPlots(study, plots = plots)

plotStudy(study, modelID = "main", featureID = "497097",
          plotID = "expression_by_cell_type")
plotStudy(study, modelID = "main", featureID = "27395",
          plotID = "expression_by_cell_type")

# Install study package and start app ------------------------------------------

installStudy(study)
message("Starting app. Should open in new browser tab.")
startApp()
