library("edgeR")
library("DESeq2")
library("biomaRt")
library("ggplot2")
#library("pheatmap")

#--------------------

args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
output <- args[2]
#glengths <- args[3]
baseName <- args[3]
#geneAnnotation <- args[5]

#--------------------
# DESeq Counts Table Formation
#--------------------

# Counts
read <- list.files(path=input, pattern="*ReadsPerGene.out.tab$", full.names=TRUE)
readFiles <- lapply(read, read.table, skip=4)
counts <- as.data.frame(sapply(readFiles, function(x) x[ , 4])) # This '4' will have to later be for user input - column to be used for counts dependent on strandedness of data
read <- gsub("_ReadsPerGene.out.tab", "", read)
read <- gsub(input, "", read) # Will need to delete the "/" before the sample name here which will fix it in conditions table
read <- gsub(".*/", "", read)
colnames(counts) <- read
row.names(counts) <- readFiles[[1]]$V1
#head(counts)

# Condition
cond <- gsub(".*_","", read)
print(cond)
colData <- data.frame(row.names=colnames(counts), condition=factor(cond, levels=c('non-COVID-19', 'COVID-19')))
#head(colData)

# DESeq DataSet
data <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~condition)
data

# Filtration
#--------------------

keep <- rowSums(counts(data)) > 0 # Drops everything with 0 reads
data <- data[keep,]

#--------------------
# Normalization & DE Analysis
#--------------------

dds <- DESeq(data, quiet=TRUE)
results <- results(dds)
#head(results)

# Results with ENSEMBL IDs:
print("Results:")
mcols(results, use.names=TRUE)
print("Summary")
summary(results)
write.csv(as.data.frame(results), file=paste(output, "/", baseName, "_ENSEMBL_results.csv", sep=""))

print("Condition vs NonCondition Contrast Results")
contrastResults <- results(dds, contrast=c("condition", "COVID-19", "non-COVID-19"))
write.csv(as.data.frame(contrastResults), file=paste(output, "/", baseName, "_condition_treated_ENSEMBL_results.csv", sep=""))

# Mapping Gene Names to ENSEMBL IDs
#--------------------

results$ensembl <- sapply(strsplit(rownames(results), split="\\+"), '[',1 )

ensembl=useMart("ensembl", dataset="hsapiens_gene_ensembl")
geneMap <- getBM( attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'), 
		 filters = 'ensembl_gene_id', 
		 values = results$ensembl, 
		 mart = ensembl)
ids <- match(results$ensembl, geneMap$ensembl_gene_id)
results$entrez_id <- geneMap$entrezgene_id[ids]
results$hgnc_symbol <- geneMap$hgnc_symbol[ids]

write.csv(as.data.frame(results), file=paste(output, "/", baseName, "_geneID_results.csv", sep=""))

print("Results:")
dim(results)

#--------------------
# Visualization & Interpretation
#--------------------

# Top Up/Down Regulated Genes
#--------------------

# TO DO: Double check this is the appropriate padj threshold
# TO DO: Create a log2FoldChange filter to determine up and down regulated genes, this just creates a list
# TO DO: Do these sigResults need to be for the contrast results or for the general results? - No as there is only one condition COVID v NONCOVID

# Using a false positive threshold of 5% / padj = 0.05

# Upregulated genes:
sigResults <- results[ which(results$padj < 0.05), ]
head(sigResults[order(sigResults$log2FoldChange),])
sigResults$diffexpressed <- NA
sigResults$diffexpressed[sigResults$log2FoldChange > 1.0] <- "UP" #Play around with this threshold?

# Downregulated genes:
tail(sigResults[order(sigResults$log2FoldChange),])
sigResults$diffexpressed[sigResults$log2FoldChange < -1.0] <- "DOWN"

write.csv(as.data.frame(sigResults), file=paste(output, "/", baseName, "_geneID_condition_treated_sigResults.csv", sep=""))

print("Significant Results:")
dim(sigResults)


# PCA
#--------------------

d <- rlogTransformation(dds, blind=TRUE) # Could VST also be used here?
png(file=paste(output, "/", baseName, "_PCA.png", sep=""))
plotPCA(d)
dev.off()

# MA Plot
#--------------------

png(file=paste(output, "/", baseName, "_MA_Plot.png", sep=""), width=2048, height=2048)
plotMA(results) # , ylim=c(-1,1)
dev.off()

# Dispersion Plot
#--------------------

png(file=paste(output, "/", baseName, "_Dispersion_Plot.png", sep=""), width=2048, height=2048)
plotDispEsts(dds, ylim=c(1e-6,1e1)) #ylim based on default for DESeq
dev.off()

# HeatMap
#--------------------

# TO DO: See how this plot outputs as it uses dds but also the geneID mapped results table, need to make sure everything maps correctly

#d <- vst(dds, blind=FALSE)
#dim(d)
#df <- as.data.frame(colData(dds)[,"condition"])
#selectCounts <- head(order(rowVars(assay(d)), decreasing=TRUE), 5)
#topCounts <- assay(d)[selectCounts,]

#png(file=paste(output, "/", baseName, "_significantGenes_Heatmap.png", sep=""), width=2048, height=2048, pointsize=90)
#pheatmap(topCounts, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
#dev.off()

# Volcano Plot
#--------------------

results$diffexpressed <- NA
results$diffexpressed[results$log2FoldChange > 1.0 & results$padj < 0.05] <- "UP"
results$diffexpressed[results$log2FoldChange < -1.0 & results$padj < 0.05] <- "DOWN"

png(file=paste(output, "/", baseName,"_significantGenes_VolcanoPlot.png", sep=""))
ggplot(results) +
       geom_point(aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=hgnc_symbol)) +
       xlab("log2 Fold Change") +
       ylab("-log10 Adjusted P-Value") +
       geom_hline(yintercept=-log10(0.05), col="red")
dev.off()
