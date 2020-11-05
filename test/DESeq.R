## Some Good reads on deseq2: https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
# Others: https://lashlock.github.io/compbio/R_presentation.html
# Others: http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
# Harvard training: https://hbctraining.github.io/Intro-to-R-with-DGE/schedule/
# DataCamp
# NOTE: DESeq uses Unnormalized counts meaning that the counts were not adjusted for the library size, which is the total number of reads counted per sample.

library("DESeq2")
library("ggplot2")
library("pheatmap")
library("reshape2")
library("ExpressionNormalizationWorkflow")
library("sva")
library("dplyr")
library("hexbin")
library("PoiClaClu")
library("RColorBrewer")
library("RUVSeq")
library("ggbeeswarm")
library("genefilter")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("Biobase")
library("limma")
library("pvca")
library("ReportingTools")
library("ashr")
library(EnhancedVolcano)
## Reading samplesheet and saving it into sampleTable Object

sampleTable <- read.csv("samplesheet.csv",  header = T, stringsAsFactors = FALSE, sep = ",")
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

tail(counts(ddsHTSeq))
counts(ddsHTSeq)

ddsHTSeq
## You can check number of rows by 
nrow(ddsHTSeq)

#######################################################
######################### Results After filtering
#class: DESeqDataSet 
#dim: 23728 6 
#metadata(1): version
#assays(1): counts
#rownames(23728): A1BG A1BG-AS1 ... ZYX ZZEF1
#rowData names(0):
#colnames(6): control1 test1 ... control3 test3
#colData names(2): Sampletype Batch
#################################################################
######################################################################## QCying the data before DESeq normalization (Optional) ###########################################################################


##To plot log2 normalized read counts for both control and test
rawcounts <-counts(ddsHTSeq) ##store roaw counts from ddsHTseq
rawcounts <- data.frame(rawcounts) ##saving it in a dataframe format
log2_rawcounts<-data.frame(log2(rawcounts+1)) ## log2 normalization of raw reads
df_boxplot <-melt(log2_rawcounts, variable.names("Samplename")) #variable.names:name of variable used to store values# melting function convert multivariate data into single column, Read here:http://www.datasciencemadesimple.com/melting-casting-r/
df_boxplot <- data.frame(df_boxplot, Condition = substr(df_boxplot$variable, 1, 4)) ## Substring will add a new column and will write string of 1st four characters (cont $ test)
ggplot(df_boxplot, aes(x= variable, y = value, fill = Condition)) + geom_boxplot(alpha = 0.4) + scale_fill_manual(values = c("#619CFF", "#F564E3")) + theme_classic() + xlab("samples") + ylab("log2 converted raw counts")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(df_boxplot, aes(x = value, colour = variable, fill = variable)) + ylim(c(0, 0.25)) +  geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ Condition) +  theme(legend.position = "top") + xlab(expression(log[2](count + 1))) # density plots

##making and understanding MA Plots for control1 & control 2 to check if gene expression among biological replicates is uniform
y = log2_rawcounts[, 2]
x = log2_rawcounts[, 1]
## M-values
M = x - y
A = (x + y)/2

df_ma = data.frame(A, M)
ggplot(df_ma, aes(x = A, y = M)) + geom_point(size = 1.5, alpha = 1/5) + geom_hline(yintercept= 0, color = "blue3") + stat_smooth(se = FALSE, method = "loess", color = "red3")



################################################################## QC Finished ###############################################################

################################################################## Normalizing using DESeq2 for various factors (Mandatory) ###############################################
#To calculate the normalized counts with DESeq2, we use estimatesizefactors() function. By assigning the results back to the dds object we are filling in the slots of the DESeqDataSet object with the appropriate information. 
ddsHTSeq <- estimateSizeFactors(ddsHTSeq)
#We can take a look at the normalization factor applied to each sample using:
sizeFactors(ddsHTSeq)
#Now these size factors are used by DESeq to normalize Raw counts. The raw counts of each sample are divided by their respective size factor for normalization
#Now, to retrieve the normalized counts matrix from dds, we use the counts() function and add the argument normalized=TRUE.
normalized_counts <- counts(ddsHTSeq, normalized=TRUE)
#We can save this normalized data matrix to file for later use:
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)

#NOTE: DESeq2 doesn’t actually use normalized counts, rather it uses the raw counts and models the normalization inside the Generalized Linear Model (GLM). 
#DESEq uses a "median of ratios" methof for normalization. This method adjusts the raw counts for library size and is resistant to large numbers of differentially expressed genes.
#These normalized counts will be useful for downstream visualization of results, but cannot be used as input to DESeq2 or any other tools that peform differential expression analysis which use the negative binomial model.
#######################################################################Normalization finished ########################################################

############################################## Unsupervised clustering analysis 

vsd_qc <- vst(ddsHTSeq, blind = FALSE)
vsd_qc_mat <- assay(vsd_qc)
vsd_qc_cor <- cor(vsd_qc_mat)
vsd_qc_cor
pheatmap(vsd_qc_cor)




#VST (Variance Stabilising Transformation) looks at the trend between variance and mean in the data, and then tries to find a strictly monotonous transformation of the data so that this trend is removed. In practice, the transformation will approach the logarithm function for high values and the square root function for small values (incl. 0), and smoothly interpolate inbetween.
vsd <- vst(ddsHTSeq, blind = FALSE)
# view first three line in vsd assay (chage ,3 to 5 or more to view more lines). Clearly the VST estimates are nearly equal to log2 of rawcounts
# Note: For genes with high counts, both the VST and the rlog will give similar result to the ordinary log2 transformation of normalized counts. For genes with lower counts, however, the values are shrunken towards a middle value.
###The VST is much faster to compute and is less sensitive to high count outliers than the rlog. The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST when there is a wide range of sequencing depth across samples
####The two transformations offered by DESeq2 are provided for applications other than differential testing.
head(assay(vsd), 3)
## To view size factor
colData(vsd)
#Reason for blind-FALSE here:https://www.biostars.org/p/333597/#440342
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


##Another option for calculating sample distances is to use the Poisson Distance. The PoissonDistance function takes the original count matrix (not normalized) 
poisd <- PoissonDistance(t(counts(ddsHTSeq)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- colnames(assay(ddsHTSeq))
colnames(samplePoisDistMatrix) <- colnames(assay(ddsHTSeq))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
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

##Another plot, very similar to the PCA plot, can be made using the multidimensional scaling (MDS) function in base R. This is useful when we don't have a matrix of data, but only a matrix of distances.
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
## The fold change is always in reference to something. Here the fold change of Silenced TRF2 samples over the controlled sample
#results of WaltTesting can be extracted using results(dds, alpha =0.5), where alpha value specify a significance value . You can choose alpha values on how stringent you want to be with your analysis
#Lower alpha value indicates less probability of ifentifying a gene as differentialy expressed when it is actually not. 0.05 is standard. A log2 fold change threshold of 0.32 can also be used.
res <- results(dds, contrast=c("Sampletype","test","control")) ## test will be in numerator and control in denominator
res

#to understand the meaning of the columns:

mcols(res, use.names = TRUE)

summary(res)
#############################################################
#out of 27946 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 401, 1.4%
#LFC < 0 (down)     : 127, 0.45%
#outliers [1]       : 30, 0.11%
#low counts [2]     : 7044, 25%
#(mean count < 7)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
######################################################################

#More Strict Criteria: There are two ways to be more strict about which set of genes are considered significant:
#1. lower the false discovery rate threshold (the threshold on padj in the results table) ie. lower the value of alpha; eg: res.05 <- results(dds, alpha = 0.05), summary(res), table(res.05$padj < 0.05)
#2.raise the log2 fold change threshold from 0 using the lfcThreshold argument of results; eg: resLFC1 <- results(dds, lfcThreshold=1), table(resLFC1$padj < 0.1)
#If we lower the false discovery rate threshold, we should also inform the results() function about it, so that the function can use this threshold for the optimal independent filtering that it performs:


###########Plotting results##########################
## use the plotCounts function that takes as arguments the DESeqDataSet, a gene name, and the group over which to plot the counts
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("Sampletype"))
##Normalized counts for a single gene over treatment group.
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("Sampletype","time"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Sampletype, y = count, color = time)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
#Normalized counts with lines connecting the groups according to time period
ggplot(geneCounts, aes(x = Sampletype, y = count, color = time, group = time)) +
  scale_y_log10() +  geom_beeswarm(size = 3) + geom_line()
##Histogram of p values for genes with mean normalized count larger than 1.
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
##Gene clustering ## Clustering top 100 highly expressed genes (,100)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 100)
topVarGenes_20 <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat_20  <- assay(vsd)[ topVarGenes_20, ]
mat <- mat - rowMeans(mat)
mat_20  <- mat_20 - rowMeans(mat_20)
anno <- as.data.frame(colData(vsd)[, c("Sampletype","time")])
pheatmap(mat, annotation_col = anno, border_color = NA, show_rownames = FALSE)
pheatmap(mat_20, annotation_col = anno, show_rownames = TRUE)

###############################Making MA Plots################################################
resultsNames(dds)
##Output: [1] "Intercept"                  "Sampletype_test_vs_control"
res.MA <- lfcShrink(dds, coef="Sampletype_test_vs_control", type="apeglm")
summary(res.MA)
###############################
#We found similar results
#out of 27946 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 401, 1.4%
#LFC < 0 (down)     : 127, 0.45%
#outliers [1]       : 30, 0.11%
#low counts [2]     : 7044, 25%
#(mean count < 7)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
#########################
DESeq2::plotMA(res.MA, ylim = c(-5,5))
## Without applying lfcshrink that remove biases and noise
res.noshr <- results(dds, name="Sampletype_test_vs_control")
plotMA(res.noshr, ylim = c(-5, 5))
#####################################################Volcano Plot ####################################

EnhancedVolcano(res.MA,
                lab = rownames(res.MA),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 8), pCutoff = 10e-14, FCcutoff = 1.0, labSize = 4.0, selectLab = 'TERF2')
				
				
##Labeling the top genes in MA Plot
topGene.MA <- rownames(res.MA)[which.min(res.MA$padj)]
with(res.MA[topGene.MA, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})



#################################Annotating and exporting results#######################
####Adding gene symbols
columns(org.Hs.eg.db) ##To get a list of all available ENtrez Gene IDs
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
resOrdered <- res[order(res$pvalue),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)
##resOrderedDF <- as.data.frame(resOrdered)[1:100, ] to get only top 100 genes
write.csv(resOrderedDF, file = "results_without_batch.csv")

## To generate an HTML version of Differentially expressed genes
htmlRep <- HTMLReport(shortName="report", title="My report",
                      reportDirectory="./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)

####ENSG to gene symbol
gene_annotate<-cbind(rownames(res), res$symbol)
colnames(gene_annotate)<-c("ENSG_ID", "Gene_Symbol")
###################################To plot heat map of log2 of normalised read counts generated by deseq2 ########################
###Extract the normalized counts
normalized_counts<-counts(ddsHTSeq, normalized =  TRUE)
log2_normalized_counts<-log2(normalized_counts+1)
############???????????????????????????????######################
#idx_gene_name<-gene_annotate[sapply(row.names(normalized_counts), function(i){return(which(i==gene_annotate[,1]))}),2]
#row.names(normalized_counts)<-idx_gene_name
#write.csv(log2_normalized_counts, file = "log2_normalized_counts_without_age.csv")
###################################################################

#Note:log2FCgiven by results are usually shrunken and should not be directly cob=nverted to FC. To do that use results(deseq, addMLE = T), This will add a column named lfcMLE in output, which should be closer to the log2FC calculated from normalized data.
#refer to this blog: https://www.biostars.org/p/248486/#248495 , https://www.biostars.org/p/331166/, https://www.biostars.org/p/101727/, http://rstudio-pubs-static.s3.amazonaws.com/13988_bb11d85b79b2436280de434988558140.html
df_normalized <-melt(log2_normalized_counts, variable.names("Samples"))

df_normalized <- data.frame(df_normalized, Condition = substr(df_normalized$Var2, 6, 6))
ggplot(df_normalized, aes(x= Var2, y = value, fill = Condition)) + geom_boxplot(alpha = 0.4) + scale_fill_manual(values = c("#619CFF", "#F564E3")) + theme_classic() + xlab("samples") + ylab("log2 converted raw counts")

#log2normalized counts of filtered reads with p value <=0.05
gene_filtered<-row.names(res)[which(res$pvalue<=0.05)]
idx_filtered<-which(row.names(log2_normalized_counts) %in% gene_filtered)
log2_normalized_counts_filtered<-log2_normalized_counts[idx_filtered,]
pheatmap(log2_normalized_counts_filtered, fontsize = 5, cellwidth = 20, show_rownames = FALSE, border_color = NA)

###Distance plot of significant genes filtered with p value <=0.05
sampleDists_res <- dist(t(log2_normalized_counts_filtered))
sampleDists_res
sampleDistMatrix_res <- as.matrix( sampleDists_res )
colnames(sampleDistMatrix_res) <- colnames(log2_normalized_counts_filtered)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix_res,
         clustering_distance_rows = sampleDists_res,
         clustering_distance_cols = sampleDists_res,
         col = colors)



#########################################################################################
##########Batch correction######################

###SVA analysis using Expression Normalization workflow package
inpData <- expSetobj(log2_normalized_counts, sampleTable)
## Set the covariates whose effect size on the data needs to be calculated
cvrts_eff_var <- c("Sampletype", "Batch")
## Set a PVCA Threshold Value between 0 & 1
## PVCA Threshold Value is the percentile value of the minimum amount of the variabilities that the selected principal components need to explain, here requiring 75% of the expression variance to be captured by the PCs
pct_thrsh <- 0.75 
## Perform the PVCA analysis
pvcAnaly(inpData, pct_thrsh, cvrts_eff_var)

##Surrogate variable analysis
## Choose a biological variable that is to be used to calculate the surrogate variables
biol_var_sva <- "Sampletype" 
##If the sva function is called without the n.sv argument specified, the number of factors will be estimated for you. The number of factors can also be estimated using the num.sv.
#sur_var_obj <- surVarAnaly(inpData, biol_var_sva, n.sv = 5)
##Eg explicitely defining n.sv =5 threw error: Error in solve.default(t(mod) %*% mod) : system is computationally singular: reciprocal condition number = 5.98609e-18
sur_var_obj <- surVarAnaly(inpData, biol_var_sva)
########Output###
#Number of significant surrogate variables is:  3 
#Iteration (out of 5 ):1  2  3  4  5  

## The newly generated surrogate variables sv1 through sv4 are appended to the ExpressionSet object
inpData_sv <- sur_var_obj$expSetobject

#PVCA of the raw data with the surrogate variables included as covariates. Keep the number of surrogate variables as described in output (Here 3)
var_names <- c("sv1","sv2","sv3") 
## View them appended to the covariate matrix as additional covariate columns
View(pData(inpData_sv))
## Include the surrogate variables as covariates in addition to BMI3, Rin3, CAD and Study (be sure to use categorical measures of BMI and Rin rather than the continuous measure)
pData(inpData_sv)<-conTocat(pData(inpData_sv), var_names) 
## Again set the PVCA Threshold to explain 75% of the expression variation
cvrts_eff_var <- c("time", "Sampletype", "sv1_cat","sv2_cat","sv3_cat")
pct_thrsh <- 0.75
pvcAnaly(inpData_sv, pct_thrsh, cvrts_eff_var) 

##check without time as variable
cvrts_eff_var <- c("Sampletype", "sv1_cat","sv2_cat","sv3_cat")
pct_thrsh <- 0.75 
pvcAnaly(inpData_sv, pct_thrsh, cvrts_eff_var)

#This implies that together the 3 SVs (majorly the SV1) capture as much variance as time alone.Both analyses have slightly reduced the CAD contribution relative to the analysis without surrogate variables.
pdata_final<-pData(inpData_sv)
pdata_final<-pdata_final[,c(-2,-6,-7:-11)]





