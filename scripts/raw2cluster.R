
#################################################
## Libraries, Functions and Variables
#################################################

library(Matrix)
library(reshape2)
library(Seurat)

readConfig <- function(file = "seurat_config.txt") {
  rawConfig <- scan(file, what = "character", sep = "\n")
  tmp <- strsplit(rawConfig, split = "=")
  seuratParams <- list()
  x <- sapply(tmp, function(x) x[1])
  y <- sapply(tmp, function(x) x[2])
  for (i in 1:length(tmp)) {
    seuratParams[[x[i]]] <- y[i]
  }
  return(seuratParams)
}

readFG <- function(countFile = "expression_data.tsv", geneFile = "gene_metadata.tsv") {
  require(Matrix)
  countData <- read.table(countFile, header = TRUE)
  geneData <- scan(geneFile, "character")[-1]
  countData$index <- match(countData[, 2], geneData)
  cellIDs <- sort(unique(countData[, 1]))
  exprsMatrix <- Matrix(matrix(data = 0, 
                               ncol = length(cellIDs),
                               nrow = length(geneData),
                               dimnames = list(geneData, cellIDs)))
  for (i in cellIDs){
    tmp <- countData[countData[, 1]==i,]
    exprsMatrix[tmp$index, as.character(i)] <- tmp[, 3]
  }
  colnames(exprsMatrix) <- paste("cellID", colnames(exprsMatrix), sep  = ".")
  return(exprsMatrix)
}

writeFG <- function(countMatrix, countFile = "vargene_expression_data.tsv", geneFile = "vargene_metadata.tsv") {
  expression_data <- data.frame(gene = rownames(countMatrix),
                                as.matrix(countMatrix),
                                stringsAsFactors = F,
                                row.names = NULL)
  expression_data <- melt(expression_data, id.vars="gene", measure.vars = colnames(countMatrix))
  expression_data <- as.matrix(expression_data[, c(2, 1, 3)])
  colnames(expression_data) <- c("cellId*Ganzzahl", "entrezId*Ganzzahl", "expressionValue*Zahl")
  expression_data <- expression_data[as.numeric(expression_data[,3])!=0,]
  expression_data <- gsub(" ", "", expression_data)
  expression_data[, 1] <- sub("cellID.", "", expression_data[, 1])
  write.table(expression_data, countFile, row.names = F, col.names = T, quote = F, sep = "\t")
  write(c("entrezId*Ganzzahl", unique(expression_data[,2])), geneFile)
}


#################################################
## Data Input
#################################################

setwd("/config/")

seuratParams <- readConfig()


setwd("/data/")

countMatrix <- readFG(countFile = seuratParams$EXPRESSION_INPUT_NAME, geneFile = seuratParams$GENEMETA_INPUT_NAME)


#################################################
## Data Analysis
#################################################

setwd("/output/")
PROJECT_NAME
seuratCounts <- new("seurat", raw.data = countMatrix)
seuratCounts <- Setup(object = seuratCounts, 
                      project = seuratParams$PROJECT_NAME, 
                      min.cells = as.numeric(seuratParams$MIN_CELLS), 
                      min.genes = as.numeric(seuratParams$MIN_GENES), 
                      is.expr = as.numeric(seuratParams$MIN_EXPRS),
                      do.logNormalize = TRUE, 
                      total.expr = 10000, 
                      do.scale = TRUE,
                      do.center = TRUE, 
                      names.field = 1, 
                      names.delim = "_",
                      meta.data = NULL, 
                      save.raw = TRUE)

seuratCounts <- AddMetaData(seuratCounts, colSums(seuratCounts@data), "nGene")
seuratCounts <- RegressOut(seuratCounts, latent.vars = "nGene")

png(paste(seuratParams$PROJECT_NAME, "Mean Variance Plot.png"), width = 2000, height = 1200, res = 140)
seuratCounts <- MeanVarPlot(object = seuratCounts, 
                            fxn.x = expMean, 
                            fxn.y = logVarDivMean,
                            do.plot = TRUE, 
                            set.var.genes = TRUE, 
                            do.text = TRUE,
                            x.low.cutoff = as.numeric(seuratParams$X_LOW_CUTOFF), 
                            x.high.cutoff = as.numeric(seuratParams$X_HIGH_CUTOFF), 
                            y.cutoff = as.numeric(seuratParams$Y_LOW_CUTOFF),
                            y.high.cutoff = as.numeric(seuratParams$Y_HIGH_CUTOFF), 
                            cex.use = 0.6, 
                            cex.text.use = 0.5,
                            do.spike = FALSE, 
                            pch.use = 16, 
                            col.use = "black",
                            spike.col.use = "red", 
                            plot.both = FALSE, 
                            do.contour = TRUE,
                            contour.lwd = 3, 
                            contour.col = "white", 
                            contour.lty = 2, 
                            num.bin = 20,
                            do.recalc = TRUE)
dev.off()

writeFG(seuratCounts@raw.data[seuratCounts@var.genes, seuratCounts@cell.names],
        countFile = paste(seuratParams$PROJECT_NAME, "variable Gene Expression Data.tsv"), 
        geneFile = paste(seuratParams$PROJECT_NAME, "variable Gene Metadata.tsv"))

seuratCounts <- PCA(seuratCounts, 
                    pc.genes = seuratCounts@var.genes, 
                    do.print = FALSE)
seuratCounts <- ProjectPCA(seuratCounts, do.print = FALSE)


seuratCounts <- FindClusters(seuratCounts,
                             pc.use = 1:as.numeric(seuratParams$PCA_DIMS_FOR_TSNE), 
                             resolution = 0.6, 
                             print.output = 0, 
                             save.SNN = T)

seuratCounts <- RunTSNE(seuratCounts, 
                        k.seed = 42, 
                        dims.use = 1:as.numeric(seuratParams$PCA_DIMS_FOR_TSNE), 
                        do.fast = T)

png(paste(seuratParams$PROJECT_NAME, "Principal Component Analysis.png"), width = 1500, height = 1200, res = 140)
PCAPlot(seuratCounts, 
        dim.1 = 1, 
        dim.2 = 2)
dev.off()

png(paste(seuratParams$PROJECT_NAME, "tSNE.png"), width = 1500, height = 1200, res = 140)
TSNEPlot(seuratCounts)
dev.off()

clusterTab <- data.frame(cell = sub("cellID.", "", names(seuratCounts@ident)),
                         cluster = seuratCounts@ident)

write.table(clusterTab, 
            paste(seuratParams$PROJECT_NAME, "Cluster Assignments.tsv"),
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

q("no")


