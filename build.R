# Create OmicNavigator study from files in the directory results/

library(ggplot2)
library(gplots)
library(viridis)
library(OmicNavigator)
library(plotly)
library(heatmaply)
library(iheatmapr)

# Create a new study -----------------------------------------------------------

study <- createStudy("RNAseq123",
                     "Bioc workflow package converted to OmicNavigator",
                     version = "0.8.0")

# Models -----------------------------------------------------------------------

models <- list(
  Differential_Expression = "A standard group-means analysis of 3 mammary cell populations"
)
study <- addModels(study, models)

# Samples ----------------------------------------------------------------------

samples <- read.delim("results/samples.txt", stringsAsFactors = FALSE)
samples <- list(Differential_Expression = samples)
study <- addSamples(study, samples)

# Features ---------------------------------------------------------------------

features <- read.delim("results/features.txt", stringsAsFactors = FALSE)
# The featureIDs in the first column must be a character vector
features$entrez <- as.character(features$entrez)
features <- list(Differential_Expression = features)
study <- addFeatures(study, features)

# MetaFeatures -----------------------------------------------------------------

metaFeatures <- read.delim("results/metaFeatures.txt", stringsAsFactors = FALSE,
                           colClasses = c(entrez = "character"))
metaFeatures <- list(Differential_Expression = metaFeatures)
study <- addMetaFeatures(study, metaFeatures)

# Assays -----------------------------------------------------------------------

assays <- read.delim("results/assays.txt", stringsAsFactors = FALSE)
assays <- list(Differential_Expression = assays)
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
  Differential_Expression = list(
    BasalvsLP = BasalvsLP,
    BasalvsML = BasalvsML,
    LPvsML = LPvsML
  )
)
study <- addResults(study, results)

# Tests ------------------------------------------------------------------------

# here we can add a description of each test, along with columns (displayed as tt)
tests <- list(
  Differential_Expression = list(
    BasalvsLP = list(
      description = "Test of differential gene expression between Basal and LP cells",
      entrez = "NCBI entrez id",
      symbol = "HGNC symbol",
      chrom = "chromosomal location of gene",
      logFC = "the log2 fold change in expression between Basal and LP cells. log2(Basal/LP)",
      AveExpr = "the average log2 expression across all experiments",
      t = "limma moderated t statistic",
      P.Value = "raw p-value for the comparison",
      adj.P.Val = "adjusted p-value accounting for multiple tests. The method of Benjamini and Hochberg was used",
      B = "the log-odds that this gene is differentially expressed"
    ),
    BasalvsML = list(
      description = "Test of differential gene expression between Basal and ML cells",
      entrez = "NCBI entrez id",
      symbol = "HGNC symbol",
      chrom = "chromosomal location of gene",
      logFC = "the log2 fold change in expression between Basal and ML cells. log2(Basal/ML)",
      AveExpr = "the average log2 expression across all experiments",
      t = "limma moderated t statistic",
      P.Value = "raw p-value for the comparison",
      adj.P.Val = "adjusted p-value accounting for multiple tests. The method of Benjamini and Hochberg was used",
      B = "the log-odds that this gene is differentially expressed"
  ),
  LPvsML = list(
    description = "Test of differential gene expression between LP and ML cells",
    entrez = "NCBI entrez id",
    symbol = "HGNC symbol",
    chrom = "chromosomal location of gene",
    logFC = "the log2 fold change in expression between LP and ML cells. log2(LP/ML)",
    AveExpr = "the average log2 expression across all experiments",
    t = "limma moderated t statistic",
    P.Value = "raw p-value for the comparison",
    adj.P.Val = "adjusted p-value accounting for multiple tests. The method of Benjamini and Hochberg was used",
    B = "the log-odds that this gene is differentially expressed"
  )
  )
)

study <- addTests(study, tests)

# Linkouts to external resources for the results table -------------------------

resultsLinkouts <- list(
  Differential_Expression = list(
    entrez = c("https://www.ncbi.nlm.nih.gov/gene/",
               "https://ensembl.org/Mus_musculus/Gene/Summary?g="),
    symbol = "http://www.informatics.jax.org/searchtool/Search.do?query="
  )
)
study <- addResultsLinkouts(study, resultsLinkouts)

# Custom plots -----------------------------------------------------------------

x <- getPlottingData(
  study,
  modelID = "Differential_Expression",
  featureID = "497097",
  testID = "BasalvsLP"
)

#single feature
expression_by_cell_type <- function(x) {
  ggDataFrame <- cbind(x$samples,
                       feature = as.numeric(x$assays))
  pval <- x[["results"]][["P.Value"]]
  ggplot(ggDataFrame, aes(x = .data$group, y = .data$feature, fill = .data$group)) +
    geom_boxplot(alpha = .75) +
    labs(x = "Cell type", y = "Gene expression",
         title = sprintf("%s (Entrez %s)", x$features$symbol, x$features$entrez),
         subtitle = sprintf("p-value: %.2e", pval)) +
    scale_fill_viridis(discrete = TRUE, begin = .25) +
    theme_classic()
}
expression_by_cell_type(x)

#multi-feature
IntFeatures <- study$results$Differential_Expression$BasalvsLP[1:10,1]
plottingData <- getPlottingData(study, modelID = "Differential_Expression", featureID = IntFeatures)

heatmap.custom <- function(plottingData){
  if (nrow(plottingData[["assays"]]) < 2) {
    stop("This plotting function requires at least 2 features")
  }
  plotMatrix <- as.matrix(plottingData$assays)
  colnames(plotMatrix) <- paste(plottingData$samples[['group']], plottingData$samples[['lane']], sep = "_")
  row.names(plotMatrix) <- plottingData$features$symbol
  heatmap.2(x = plotMatrix,
            trace = "none",
            key = TRUE,
            key.title = NA,
            key.xlab = "Gene Expression",
            key.ylab = NA,
            col = viridis(75),
            cexRow = 1.25,
            cexCol = 1.25,
            srtCol= 60,
            margins = c(9,8)
  )
}
heatmap.custom(plottingData)

# multiTest (singleFeature)
dataMultiTest <- getPlottingData(
  study,
  modelID = "Differential_Expression",
  featureID = "497097",
  testID = names(getTests(study, modelID = "Differential_Expression"))
)

plotTstats <- function(x){
  if (length(x$results) < 2) {
    stop("This plotting function requires results from at least 2 testIDs")
  }

  tstats <- vapply(x$results, function(z) z$t, numeric(1))

  dotchart(
    x = tstats,
    xlab = "t-statistic",
    ylab = "testID",
    main = sprintf("%s (Entrez %s)", x$features$symbol, x$features$entrez)
  )
}
plotTstats(dataMultiTest)

# multiTest (multiFeature)
dataMultiTestMultiFeature <- getPlottingData(
  study,
  modelID = "Differential_Expression",
  featureID = IntFeatures,
  testID = names(getTests(study, modelID = "Differential_Expression"))
)

heatmapTstats <- function(x){
  if (length(x$results) < 2) {
    stop("This plotting function requires results from at least 2 testIDs")
  }

  tstats <- vapply(x$results, function(z) z$t, numeric(nrow(x$features)))
  rownames(tstats) <- x$features[[1]]

  heatmap.2(x = tstats,
            trace = "none",
            key = TRUE,
            key.title = NA,
            key.xlab = "t-statistics",
            key.ylab = NA,
            col = viridis(75),
            cexRow = 1,
            cexCol = 1.25,
            srtCol= 60,
            margins = c(9,8)
  )
}
heatmapTstats(dataMultiTestMultiFeature)

# multiTest (multiFeature) with iheatmapr

dataMultiTestMultiFeature <- getPlottingData(
  study,
  modelID = "Differential_Expression",
  featureID = IntFeatures,
  testID = names(getTests(study, modelID = "Differential_Expression"))
)

iheatmapr.custom <- function(dataMultiTestMultiFeature) {
  if (nrow(dataMultiTestMultiFeature[["assays"]]) < 2) {
    stop("This plotting function requires at least 2 features")
  }
  plotMatrix <- as.matrix(dataMultiTestMultiFeature$assays)
  colnames(plotMatrix) <- paste(dataMultiTestMultiFeature$samples[['group']], dataMultiTestMultiFeature$samples[['lane']], sep = "_")
  row.names(plotMatrix) <- dataMultiTestMultiFeature$features$symbol
  # assays
  cells       <- dataMultiTestMultiFeature$samples
  cells_array <- cells[order(cells$group), "group"]
  plotMatrix  <- plotMatrix[,order(cells$group)]

  # results
  results.ls <- list()
  res.ind    <- 1
  for (ii in seq_along(names(dataMultiTestMultiFeature$results))) {
    results_boolean <- dataMultiTestMultiFeature$results[[ii]]$P.Val >= 0.05
    results.ls[[res.ind]] <- ifelse(results_boolean, "p>=0.05", "p<0.05")
    res.ind <- res.ind + 1
  }
  names(results.ls) <- names(dataMultiTestMultiFeature$results)

  # iheatmapr
  p <- main_heatmap(plotMatrix,
                    name = "Assay",
                    colors = "Blues") %>%
    add_col_clustering(method = "groups",
                       groups = cells_array,
                       name = "Group",
                       colors = "Set1") %>%
    add_row_annotation(data.frame("BasalvsLP" = results.ls[[1]],
                                  "BasalvsML" = results.ls[[2]],
                                  "LPvsML" = results.ls[[3]]),
                       side = 'right',
                       size = 0.03,
                       colors = list("BasalvsLP" = c("#FF7F00", "white"),
                                     "BasalvsML" = c("purple", "white"),
                                     "LPvsML" = c("darkblue", "white"))) %>%
    add_row_summary(groups = TRUE,
                    type = "bar",
                    layout = list(title = "Average<br>per<br>group",
                                  font = list(size = 8)))
  p <- p %>%
    add_col_labels() %>%
    modify_layout(list(margin = list(b = 100))) %>%
    add_row_clustering() %>%
    add_row_title("Features")
  plotly::plotly_build(iheatmapr::to_plotly_list(p))
}
iheatmapr.custom(dataMultiTestMultiFeature)

# Interactive plotly plot (Single Feature)
single_feature_plotly <- function(x) {
  ggDataFrame <- cbind(x$samples,
                       feature = as.numeric(x$assays))
  pval <- x[["results"]][["P.Value"]]
  p <- ggplot(ggDataFrame, aes(x = .data$group, y = .data$feature, fill = .data$group)) +
    geom_boxplot(alpha = .75) +
    labs(x = "Cell type", y = "Gene expression",
         title = sprintf("%s (Entrez %s)", x$features$symbol, x$features$entrez),
         subtitle = sprintf("p-value: %.2e", pval)) +
    scale_fill_viridis(discrete = TRUE, begin = .25) +
    theme_classic()

  plotly::ggplotly(p)
}
single_feature_plotly(x)

# Interactive plotly plot (Multi Feature)

heatmap.plotly <- function(plottingData){
  if (nrow(plottingData[["assays"]]) < 2) {
    stop("This plotting function requires at least 2 features")
  }
  plotMatrix <- as.matrix(plottingData$assays)
  colnames(plotMatrix) <- paste(plottingData$samples[['group']], plottingData$samples[['lane']], sep = "_")
  row.names(plotMatrix) <- plottingData$features$symbol
  plot_ly(x = colnames(plotMatrix), y = rownames(plotMatrix), z = plotMatrix, type = "heatmap",
          colors = viridis(75))
}
heatmap.plotly(plottingData)

# Interactive heatmaply plot (Multi Feature)
heatmap.heatmaply <- function(plottingData){
  if (nrow(plottingData[["assays"]]) < 2) {
    stop("This plotting function requires at least 2 features")
  }
  plotMatrix <- round(plottingData$assays, 2)
  row.names(plotMatrix) <- ifelse(duplicated(plottingData$features$symbol) | is.na(plottingData$features$symbol),
                                  plottingData$features$entrez,
                                  plottingData$features$symbol) #avoid duplicated gene names or NAs
  heatmaply(
    x = plotMatrix,
    seriate = "OLO",
    label_names = c("Gene", "Sample", "log<sub>2</sub>(CPM)"),
    key.title = "Gene Expression
    (log<sub>2</sub>CPM)",
    col_side_colors = plottingData$samples[, c("group", "lane")],
    row_side_colors = plottingData$features[, c("chrom"), drop = F], #if you keep it a data.frame heatmaply will retain label
    plot_method = "ggplot"
    )
}
heatmap.heatmaply(plottingData)


plots <- list(
  Differential_Expression = list(
    expression_by_cell_type = list(
      displayName = "Expression by cell type",
      description = "boxplot of gene expression by cell type",
      plotType = "singleFeature",
      packages = c("ggplot2", "viridis")
    ),
    heatmap.custom = list(
      displayName = "Expression Heatmap",
      description = "heatmap of gene expression using gplots package",
      plotType = "multiFeature",
      packages = c("gplots", "viridis")
    ),
    plotTstats = list(
      displayName = "t-statistics",
      description = "plot of t-stats for each comparison",
      plotType = "multiTest"
    ),
    heatmapTstats = list(
      displayName = "t-statistics",
      description = "heatmap of t-stats for each comparison using gplots package",
      plotType = c("multiFeature", "multiTest"),
      packages = c("gplots", "viridis")
    ),
    single_feature_plotly = list(
      displayName = "Expression by cell type interactive with plotly",
      description = "interactive plotly boxplot of gene expression by cell type",
      plotType = c("singleFeature", "plotly"),
      packages = c("ggplot2", "viridis", "plotly")
    ),
    heatmap.plotly = list(
      displayName = "Expression Heatmap interactive with plotly",
      description = "interactive heatmap of gene expression created with the plotly package. Click on a row to highlight a gene and box select to zoom.",
      plotType = c("multiFeature", "plotly"),
      packages = c("viridis", "plotly")
    ),
    heatmap.heatmaply = list(
      displayName = "Expression Heatmap interactive with heatmaply",
      description = "interactive heatmap of gene expression created with the heatmaply package. Click on a row to highlight a gene and box select to zoom.",
      plotType = c("multiFeature", "plotly"),
      packages = c("viridis", "heatmaply")
    ),
    iheatmapr.custom = list(
      displayName = "Expression Heatmap interactive with iheatmapr",
      description = "interactive heatmap of gene expression created with iheatmapr. Click on a row to highlight a gene and box select to zoom. Vertical colorbars to the right of the heatmap indicate DE calls at p < 0.5. The average expression per group is displayed as a colored bargraph furthest to the right.",
      plotType = c("multiFeature", "multiTest", "plotly"),
      packages = c("plotly", "iheatmapr")
    )
  )
)
study <- addPlots(study, plots = plots)

plotStudy(study, modelID = "Differential_Expression", featureID = "497097",
          plotID = "expression_by_cell_type")
plotStudy(study, modelID = "Differential_Expression", featureID = "27395",
          plotID = "expression_by_cell_type")
plotStudy(study, modelID = "Differential_Expression", featureID = IntFeatures,
          plotID = "heatmap.custom")
plotStudy(study, modelID = "Differential_Expression", featureID = "27395",
          plotID = "plotTstats",
          testID = names(getTests(study, modelID = "Differential_Expression")))
plotStudy(study, modelID = "Differential_Expression", featureID = IntFeatures,
          plotID = "heatmapTstats",
          testID = names(getTests(study, modelID = "Differential_Expression")))
jsonBoxplot <- plotStudy(study, modelID = "Differential_Expression", featureID = "497097",
          plotID = "single_feature_plotly")
jsonHeatmap <- plotStudy(study, modelID = "Differential_Expression", featureID = IntFeatures,
          plotID = "heatmap.plotly")
jsonHeatmaply <- plotStudy(study, modelID = "Differential_Expression", featureID = IntFeatures,
                         plotID = "heatmap.heatmaply")
jsonHeatmapr <- plotStudy(study, modelID = "Differential_Expression", featureID = IntFeatures,
                          testID = names(getTests(study, modelID = "Differential_Expression")),
                          plotID = "iheatmapr.custom")

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
  Differential_Expression = list(
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
  Differential_Expression = "results/report.html"
)
study <- addReports(study, reports)

# Install study package and start app ------------------------------------------

installStudy(study)
if (interactive()) {
  message("Starting app. Should open in new browser tab.")
  startApp()
}
