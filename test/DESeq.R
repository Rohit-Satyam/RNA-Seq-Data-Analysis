## Some Good reads on deseq2: https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
# Others: https://lashlock.github.io/compbio/R_presentation.html
# Others: http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html


library(DESeq2)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(ExpressionNormalizationWorkflow)
library(sva)
library(dplyr)
library(hexbin)
library(PoiClaClu)
library(RColorBrewer)
library(RUVSeq)
library(ggbeeswarm)
library(genefilter)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ExpressionNormalizationWorkflow)
library(Biobase)
library(limma)
library(pvca)
## Reading samplesheet and saving it into sampleTable Object

sampleTable<-read.csv("samplesheet.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")
#saving row names of all samples
rownames(sampleTable)<-sampleTable[,1]
rownames(sampleTable)

##Saving in deseq object
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = ".", design = ~ Sampletype)

## Prefiltering: Here via filtering, we try to reduce the memory size of the dds data object, and we increase the speed of the transformation and testing functions within DESeq2. Here we perform a minimal pre-filtering to keep only rows that have at least 1 read.
## This was done since I had six columns for each row. Setting it to > 1 will results in some samples with zero read counts.
## more stringent # at least 3 samples with a count of 10 or higher
#keep <- rowSums(counts(dds) >= 10) >= 3

keep <- rowSums(counts(ddsHTSeq)) > 10
ddsHTSeq <- ddsHTSeq[keep,]

#######################################################
######################### Results Before filtering
#class: DESeqDataSet 
#dim: 60662 6 
#metadata(1): version
#assays(1): counts
#rownames(60662): ENSG00000000003.15 ENSG00000000005.6 ... ENSG00000288587.1 ENSG00000288588.1
#rowData names(0):
#colnames(6): Control1 Control2 ... TRF2-shRNA2 TRF2-shRNA3
#colData names(2): Sampletype time
#########################Results after filtering
#class: DESeqDataSet 
#dim: 27946 6 
#metadata(1): version
#assays(1): counts
#rownames(27946): ENSG00000000003.15 ENSG00000000419.12 ... ENSG00000288582.1 ENSG00000288586.1
#rowData names(0):
#colnames(6): Control1 Control2 ... TRF2-shRNA2 TRF2-shRNA3
#colData names(2): Sampletype time
#################################################################
ddsHTSeq
## You can check number of rows by 
nrow(ddsHTSeq)
##To plot log3 notmalized read counts for both control and test
rawcounts <-counts(ddsHTSeq) ##store roaw counts from ddsHTseq
rawcounts <- data.frame(rawcounts) ##saving it in a dataframe format
log2_rawcounts<-data.frame(log2(rawcounts+1)) ## log2 normalization of raw reads
df_boxplot <-melt(log2_rawcounts, variable.names("Samples")) ## melting function convert multivariate data into single column, Read here:http://www.datasciencemadesimple.com/melting-casting-r/
df_boxplot <- data.frame(df_boxplot, Condition = substr(df_boxplot$variable, 6, 6))
ggplot(df_boxplot, aes(x= variable, y = value, fill = Condition)) + geom_boxplot(alpha = 0.4) + scale_fill_manual(values = c("#619CFF", "#F564E3")) + theme_classic() + xlab("samples") + ylab("log2 converted raw counts")

#VST looks at the trend between variance and mean in the data, and then tries to find a strictly monotonous transformation of the data so that this trend is removed. In practice, the transformation will approach the logarithm function for high values and the square root function for small values (incl. 0), and smoothly interpolate inbetween.
vsd <- vst(ddsHTSeq, blind = FALSE)
# view first three line in vsd assay (chage ,3 to 5 or more to view more lines). Clearly the VST estimates are nearly equal to log2 of rawcounts
# Note: For genes with high counts, both the VST and the rlog will give similar result to the ordinary log2 transformation of normalized counts. For genes with lower counts, however, the values are shrunken towards a middle value.
###The VST is much faster to compute and is less sensitive to high count outliers than the rlog. The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST when there is a wide range of sequencing depth across samples
####The two transformations offered by DESeq2 are provided for applications other than differential testing.
head(assay(vsd), 3)
## To view size factor
colData(vsd)

rld <- rlog(ddsHTSeq, blind = FALSE)
head(assay(rld), 3)
colData(rld)

df <- bind_rows(
  as_data_frame(log2(counts(ddsHTSeq, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
colnames(df)[1:2] <- c("x", "y")

##Shown are scatterplots using the log2 transform of normalized counts (left), using the VST (middle), and using the rlog (right). While the rlog is on roughly the same scale as the log2 counts, the VST has a upward shift for the smaller values. It is the differences between samples (deviation from y=x in these scatterplots) which will contribute to the distance calculations and the PCA plot.
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

###Sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
# colnames(sampleDistMatrix) <- NULL  ##not really required
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists, col = colors)

##Another option for calculating sample distances is to use the Poisson Distance 
poisd <- PoissonDistance(t(counts(ddsHTSeq)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- colnames(assay(ddsHTSeq))
colnames(samplePoisDistMatrix) <- colnames(assay(ddsHTSeq))
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

###PCA plot
plotPCA(vsd, intgroup = c("Sampletype")) ## Using function plotPCA that comes with DESeq2. The two terms specified by intgroup are the interesting groups for labeling the samples; they tell the function to use them to choose colors.
##build the PCA plot from scratch using the ggplot2 package
pcaData <- plotPCA(vsd, intgroup = c( "Sampletype", "time"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = time, shape=Sampletype)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

##Another plot, very similar to the PCA plot, can be made using the multidimensional scaling (MDS) function in base R. This is useful when we don’t have a matrix of data, but only a matrix of distances.
mds <- as.data.frame(colData(vsd))  %>% cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = time, shape = Sampletype)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")
mdsPois <- as.data.frame(colData(ddsHTSeq)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = time, shape = Sampletype)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")

##Running the differential expression pipeline
dds <- DESeq(ddsHTSeq)

################################################################
#Log of DESeq calculation above
#using pre-existing size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing

#######estimateSizeFactors
#This calculates the relative library depth of each sample 

########estimateDispersions
#estimates the dispersion of counts for each gene 

##########nbinomWaldTest
#calculates the significance of coefficients in a Negative Binomial GLM using the size and dispersion outputs
##################################################
