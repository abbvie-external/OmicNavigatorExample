# Create OmicNavigator study

library(ggplot2)
library(OmicNavigator)

study <- createStudy("RNAseq123",
                     "Bioc workflow package converted to OmicNavigator",
                     version = "0.1.0")

models <- list(
  main = "A standard group-means analysis of 3 mammary cell populations"
)
study <- addModels(study, models)

samples <- read.delim("results/samples.txt", stringsAsFactors = FALSE)
samples <- list(main = samples)
study <- addSamples(study, samples)

features <- read.delim("results/features.txt", stringsAsFactors = FALSE)
# The featureIDs in the first column must be a character vector
features$entrez <- as.character(features$entrez)
features <- list(main = features)
study <- addFeatures(study, features)

metaFeatures <- read.delim("results/metaFeatures.txt", stringsAsFactors = FALSE,
                           colClasses = c(entrez = "character"))
metaFeatures <- list(main = metaFeatures)
study <- addMetaFeatures(study, metaFeatures)

assays <- read.delim("results/assays.txt", stringsAsFactors = FALSE)
assays <- list(main = assays)
study <- addAssays(study, assays)

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

tests <- list(
  main = list(
    BasalvsLP = "Which genes are DE between Basal and LP cells?",
    BasalvsML = "Which genes are DE between Basal and ML cells?",
    LPvsML = "Which genes are DE between LP and ML cells?"
  )
)
study <- addTests(study, tests)

resultsLinkouts <- list(
  main = list(
    entrez = c("https://www.ncbi.nlm.nih.gov/gene/",
               "https://ensembl.org/Mus_musculus/Gene/Summary?g="),
    symbol = "http://www.informatics.jax.org/searchtool/Search.do?query="
  )
)
study <- addResultsLinkouts(study, resultsLinkouts)

load("data/mouse_H_v5p2.rdata")
annotations <- list(
  mouse_H_v5p2 = list(
    terms = Mm.H,
    description = "The gene signatures from the Broad Institute's MSigDB Hallmark collection",
    featureID = "entrez"
  )
)
study <- addAnnotations(study, annotations)

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

enrichmentsLinkouts <- list(
  mouse_H_v5p2 = "https://www.gsea-msigdb.org/gsea/msigdb/cards/"
)
study <- addEnrichmentsLinkouts(study, enrichmentsLinkouts)

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



reports <- list(
  main = "https://github.com/abbvie-external/OmicNavigatorExample/blob/main/analyze.R"
)
study <- addReports(study, reports)

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

exportStudy(study)
