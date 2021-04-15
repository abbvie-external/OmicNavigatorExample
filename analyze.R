# Analyze RNA-seq experiment with limma+voom

# Prepare data -----------------------------------------------------------------

# Load required packages
library("limma")
library("edgeR")
library("Mus.musculus")

# Import RNA-seq counts
files <- Sys.glob("data/*txt")
x <- readDGE(files, columns = c(1, 3))

# Organize sample data
samplenames <- basename(colnames(x))
colnames(x) <- samplenames
x$samples <- cbind(samplenames, x$samples)
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP",
                     "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004", "L006", "L008"), c(3, 4, 2)))
x$samples$lane <- lane

# Organize feature data
geneid <- rownames(x)
genes <- OrganismDbi::select(
  Mus.musculus,
  keys = geneid,
  columns = c("SYMBOL", "TXCHROM"),
  keytype = "ENTREZID"
)
genes <- genes[!duplicated(genes$ENTREZID), ]
x$genes <- genes
colnames(x$genes) <- c("entrez", "symbol", "chrom")

# Organize metaFeature data
metaFeatures <- OrganismDbi::select(
  Mus.musculus,
  keys = x$genes$entrez,
  columns = c("ENSEMBL", "ENSEMBLTRANS", "ENSEMBLPROT"),
  keytype = "ENTREZID"
)
colnames(metaFeatures)[1] <- "entrez"

# Filter and normalize counts
cpm <- cpm(x)
keep.exprs <- rowSums(cpm > 1) >= 3
x <- x[keep.exprs, , keep.lib.sizes = FALSE]
x <- calcNormFactors(x, method = "TMM")

# Differential expression analysis ---------------------------------------------

# Specify linear model and contrasts to test
design <- model.matrix(~0 + group + lane, data = x$samples)
colnames(design) <- sub("(group|lane)", "", colnames(design))
contrastsMatrix <- makeContrasts(BasalvsLP = Basal - LP,
                                 BasalvsML = Basal - ML,
                                 LPvsML = LP - ML,
                                 levels = colnames(design))

# Fit the model
v <- voom(x, design)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contrastsMatrix)
efit <- eBayes(vfit)

# Enrichment analysis ----------------------------------------------------------

# Load Mm.h, which contains the MSigDB Hallmark gene sets for the enrichment
# analysis
load("data/mouse_H_v5p2.rdata")
idx <- ids2indices(Mm.H, id = rownames(v))

# Calculate enrichment for each of the 3 tests
enrichedBasalvsLP <- camera(v, idx, design, contrast = contrastsMatrix[, 1])
enrichedBasalvsML <- camera(v, idx, design, contrast = contrastsMatrix[, 2])
enrichedLPvsML <- camera(v, idx, design, contrast = contrastsMatrix[, 3])

# Export results ---------------------------------------------------------------

dir.create("results", showWarnings = FALSE)
export <- function(x, file, row.names = FALSE, ...) {
  write.table(x, file = file, quote = FALSE, sep = "\t", row.names = row.names, ...)
}

export(x$samples, "results/samples.txt")
export(x$genes, "results/features.txt")
export(metaFeatures, "results/metaFeatures.txt")

# Convert counts to log-counts-per-million
assays <- as.data.frame(cpm(x, log = TRUE))
export(assays, "results/assays.txt", row.names = TRUE)

# Test results
create_results_table <- function(fit, coef) {
  topTable(fit, coef = coef, number = Inf, sort.by = "p")[, -2:-3]
}
export(create_results_table(efit, "BasalvsLP"), "results/BasalvsLP.txt")
export(create_results_table(efit, "BasalvsML"), "results/BasalvsML.txt")
export(create_results_table(efit, "LPvsML"), "results/LPvsML.txt")

# Enrichments
export(enrichedBasalvsLP, "results/enrichedBasalvsLP.txt", row.names = TRUE)
export(enrichedBasalvsML, "results/enrichedBasalvsML.txt", row.names = TRUE)
export(enrichedLPvsML, "results/enrichedLPvsML.txt", row.names = TRUE)
