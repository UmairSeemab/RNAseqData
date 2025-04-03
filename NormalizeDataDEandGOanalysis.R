# Install the latest version of DEseq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# load the library
library(DESeq2)

setwd("/home/umair/Documents/RNA_Seq/humanData_20250401/below60Data/readcountFiles")


# Read in the raw read counts
rawCounts <- read.delim("ngs-data-table_countFile.tsv")
head(rawCounts)

# Read in the sample mappings
sampleData <- read.delim("phenodata.tsv")
head(sampleData)

# Also save a copy for later
sampleData_v2 <- sampleData

# Convert count data to a matrix of appropriate form that DESeq2 can read
geneName <- rawCounts$GeneName
sampleIndex <- grepl("ERR\\d+", colnames(rawCounts))
rawCounts <- as.matrix(rawCounts[,sampleIndex])
rownames(rawCounts) <- geneName
head(rawCounts)

write.csv(rawCounts,"rawCounts.csv",sep = "\t", row.names = TRUE)

# Convert sample variable mappings to an appropriate form that DESeq2 can read
head(sampleData)
rownames(sampleData) <- sampleData$sample
keep <- c("group_details", "group")
sampleData <- sampleData[,keep]
colnames(sampleData) <- c("group_details", "group")
sampleData$group <- factor(sampleData$group)
head(sampleData)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts <- rawCounts[,unique(rownames(sampleData))]
all(colnames(rawCounts) == rownames(sampleData))

# Order the tissue types so that it is sensible and make sure the control sample is first: normal sample -> primary tumor -> metastatic tumor
sampleData$group_details <- factor(sampleData$group_details, levels=c("Diseased", "DrugTreated"))

# Create the DESeq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ group)

dim(deseq2Data)
dim(deseq2Data[rowSums(counts(deseq2Data)) > 5, ])

# Perform pre-filtering of the data
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 5, ]

# Run pipeline for differential expression
deseq2Data <- DESeq(deseq2Data)

# Extract differential expression results
deseq2Results <- results(deseq2Data, contrast=c("group", "1", "2"))

# View summary of results
summary(deseq2Results)

# Using DEseq2 built in method
plotMA(deseq2Results)

# Load libraries for ggplot
library(ggplot2)
library(scales) 
library(viridis)

# Coerce to a data frame
deseq2ResDF <- as.data.frame(deseq2Results)

# Examine this data frame
head(deseq2ResDF)

# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .1, "Significant", NA)

# Plot the results similar to DEseq2
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) + 
  geom_point(size=1) + 
  scale_y_continuous(limits=c(-3, 3), oob=squish) + 
  scale_x_log10() + 
  geom_hline(yintercept = 0, colour="tomato1", size=2) + 
  labs(x="mean of normalized counts", y="log fold change") + 
  scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + 
  theme_bw()

# Variance stabilizing transformation
deseq2VST <- vst(deseq2Data)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

write.csv(deseq2VST,"/home/umair/Documents/RNA_Seq/humanData_20250401/below60Data/readcountFiles/deseq2VST.csv", row.names = TRUE) # Save after variance stabilizing transform

# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 3,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

head(deseq2VST)

# Convert the VST counts to long format for ggplot2
library(reshape2)

deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

# Plot the heatmap
heatmap <- ggplot(deseq2VST_long, aes(x=variable, y=Gene, fill=value)) + 
  geom_raster() + 
  scale_fill_viridis(trans="sqrt") + 
  theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

# Save the results
write.csv(deseq2ResDF,"/home/umair/Documents/RNA_Seq/humanData_20250401/below60Data/readcountFiles/deseq2ResDF.csv", row.names = TRUE)


library(reshape2)

# First compare wide vs long version
deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))

head(deseq2VST_wide)
head(deseq2VST_long)

# Now overwrite our original data frame with the long format
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

# Make a heatmap
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap


# Convert the significant genes back to a matrix for clustering
deseq2VSTMatrix <- dcast(deseq2VST, Gene ~ variable)
rownames(deseq2VSTMatrix) <- deseq2VSTMatrix$Gene
deseq2VSTMatrix$Gene <- NULL

# Compute a distance calculation on both dimensions of the matrix
distanceGene <- dist(deseq2VSTMatrix)
distanceSample <- dist(t(deseq2VSTMatrix))

# Cluster based on the distance calculations
clusterGene <- hclust(distanceGene, method="average")
clusterSample <- hclust(distanceSample, method="average")

# Construct a dendogram for samples
#install.packages("ggdendro")
library(ggdendro)
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

# Re-factor samples for ggplot2
deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])

# Construct the heatmap. note that at this point we have only clustered the samples NOT the genes
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

# Combine the dendrogram and the heatmap
#install.packages("gridExtra")
library(gridExtra)
grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))


# Load in libraries necessary for modifying plots
#install.packages("gtable")
library(gtable)
library(grid)

# Modify the ggplot objects
sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0))
heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))

# Convert both grid based objects to grobs
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)

# Check the widths of each grob
sampleDendrogramGrob$widths
heatmapGrob$widths

# Add in the missing columns
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

# Make sure every width between the two grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)

# Arrange the grobs into a plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))

# Draw the plot
grid.draw(finalGrob)

# Re-order the sample data to match the clustering we did
sampleData_v2$Run <- factor(sampleData_v2$sample, levels=clusterSample$labels[clusterSample$order])

# Construct a plot to show the clinical data
colours <- c("#743B8B", "#8B743B")
sampleClinical <- ggplot(sampleData_v2, aes(x=Run, y=1, fill=group_details)) + geom_tile() + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)) + scale_fill_manual(name="Tissue", values=colours) + theme_void()

# Convert the clinical plot to a grob
sampleClinicalGrob <- ggplotGrob(sampleClinical)

# Make sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleClinicalGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleClinicalGrob$widths <- as.list(maxWidth)

# Arrange and output the final plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, sampleClinicalGrob, heatmapGrob, ncol=1, heights=c(2,1,5))
grid.draw(finalGrob)

write.csv(deseq2ResDF,"deseq2ResDF.csv", row.names = TRUE)
















############################### GSEA ##########################

# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("gage", "GO.db", "AnnotationDbi", "org.Hs.eg.db"))
library(gage)

# Conduct DE analysis for female vs male comparison
female_v_male_DE <- results(deseq2Data, contrast=c("group", "1", "2"))

# Set up KEGG database for human
kg.hs <- kegg.gsets(species="hsa")
kegg.sigmet.gs <- kg.hs$kg.sets[kg.hs$sigmet.idx]
kegg.dise.gs <- kg.hs$kg.sets[kg.hs$dise.idx]

# Set up GO database for human
go.hs <- go.gsets(species="Human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]

# Load required annotation packages
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)

# Annotate DESeq2 results with additional gene identifiers (ENSEMBL, ENTREZ, GeneName)
female_v_male_DE$ENSEMBL <- mapIds(org.Hs.eg.db, keys=row.names(female_v_male_DE), column="ENSEMBL", keytype="SYMBOL", multiVals="first")
female_v_male_DE$entrez <- mapIds(org.Hs.eg.db, keys=row.names(female_v_male_DE), column="ENTREZID", keytype="SYMBOL", multiVals="first")
female_v_male_DE$name <- mapIds(org.Hs.eg.db, keys=row.names(female_v_male_DE), column="GENENAME", keytype="SYMBOL", multiVals="first")

# Grab the log2 fold changes for all genes
female_v_male_DE.fc <- female_v_male_DE$log2FoldChange
names(female_v_male_DE.fc) <- female_v_male_DE$entrez

# Run GAGE enrichment analysis for all log2 fold changes
fc.kegg.sigmet.p <- gage(female_v_male_DE.fc, gsets = kegg.sigmet.gs)
fc.kegg.dise.p <- gage(female_v_male_DE.fc, gsets = kegg.dise.gs)
fc.go.bp.p <- gage(female_v_male_DE.fc, gsets = go.bp.gs)
fc.go.mf.p <- gage(female_v_male_DE.fc, gsets = go.mf.gs)
fc.go.cc.p <- gage(female_v_male_DE.fc, gsets = go.cc.gs)

# Convert the KEGG results to data frames for up-regulated pathways
fc.kegg.sigmet.p.up <- as.data.frame(fc.kegg.sigmet.p$greater)
fc.kegg.dise.p.up <- as.data.frame(fc.kegg.dise.p$greater)

# Write the results to files
write.table(fc.kegg.sigmet.p.up, "/home/umair/Documents/RNA_Seq/humanData_20250401/below60Data/readcountFiles/GSEA/fc.kegg.sigmet.p.up.txt", sep = ";", row.names = TRUE)
write.table(fc.kegg.dise.p.up, "/home/umair/Documents/RNA_Seq/humanData_20250401/below60Data/readcountFiles/GSEA/fc.kegg.dise.p.up.txt", sep = ";", row.names = TRUE)

# Convert the KEGG results for down-regulated pathways
fc.kegg.sigmet.p.down <- as.data.frame(fc.kegg.sigmet.p$less)
fc.kegg.dise.p.down <- as.data.frame(fc.kegg.dise.p$less)

# Write the results for down-regulated pathways
write.table(fc.kegg.sigmet.p.down, "/home/umair/Documents/RNA_Seq/humanData_20250401/below60Data/readcountFiles/GSEA/fc.kegg.sigmet.p.down.txt", sep = ";", row.names = TRUE)
write.table(fc.kegg.dise.p.down, "/home/umair/Documents/RNA_Seq/humanData_20250401/below60Data/readcountFiles/GSEA/fc.kegg.dise.p.down.txt", sep = ";", row.names = TRUE)

# Convert the GO results to data frames for up-regulated terms
fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
fc.go.mf.p.up <- as.data.frame(fc.go.mf.p$greater)
fc.go.cc.p.up <- as.data.frame(fc.go.cc.p$greater)

# Write the results for GO up-regulated terms
write.table(fc.go.bp.p.up, "/home/umair/Documents/RNA_Seq/humanData_20250401/below60Data/readcountFiles/GSEA/fc.go.bp.p.up.txt", sep = ";", row.names = TRUE)
write.table(fc.go.mf.p.up, "/home/umair/Documents/RNA_Seq/humanData_20250401/below60Data/readcountFiles/GSEA/fc.go.mf.p.up.txt", sep = ";", row.names = TRUE)
write.table(fc.go.cc.p.up, "/home/umair/Documents/RNA_Seq/humanData_20250401/below60Data/readcountFiles/GSEA/fc.go.cc.p.up.txt", sep = ";", row.names = TRUE)

# Convert the GO results for down-regulated terms
fc.go.bp.p.down <- as.data.frame(fc.go.bp.p$less)
fc.go.mf.p.down <- as.data.frame(fc.go.mf.p$less)
fc.go.cc.p.down <- as.data.frame(fc.go.cc.p$less)

# Write the results for GO down-regulated terms
write.table(fc.go.bp.p.down, "/home/umair/Documents/RNA_Seq/humanData_20250401/below60Data/readcountFiles/GSEA/fc.go.bp.p.down.txt", sep = ";", row.names = TRUE)
write.table(fc.go.mf.p.down, "/home/umair/Documents/RNA_Seq/humanData_20250401/below60Data/readcountFiles/GSEA/fc.go.mf.p.down.txt", sep = ";", row.names = TRUE)
write.table(fc.go.cc.p.down, "/home/umair/Documents/RNA_Seq/humanData_20250401/below60Data/readcountFiles/GSEA/fc.go.cc.p.down.txt", sep = ";", row.names = TRUE)

# Install pathview for pathway visualization
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")
setwd("/home/umair/Documents/RNA_Seq/humanData_20250401/below60Data/readcountFiles/GSEA/")

library(pathview)

# Display KEGG pathway example (e.g., mmu04744)
pathview(gene.data=female_v_male_DE.fc, species="hsa", pathway.id="hsa00071")

# Display another KEGG pathway example (e.g., mmu04080)
pathview(gene.data=female_v_male_DE.fc, species="hsa", pathway.id="hsa00071", kegg.native=FALSE)
