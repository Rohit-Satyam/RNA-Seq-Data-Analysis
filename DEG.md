# Sources Referred
1. Harvard DEG tutorial: [here](https://hbctraining.github.io/DGE_workshop/lessons/01_DGE_setup_and_overview.html) 

The data we currently use in this tutorial can be found [here](https://github.com/Rohit-Satyam/RNA-Seq-Data-Analysis/tree/master/test)

## Differential gene expression analysis overview

So what does this count data actually represent? The count data used for differential expression analysis represents the number of sequence reads that originated from a particular gene. The higher the number of counts, the more reads associated with that gene, and the assumption that there was a higher level of expression of that gene in the sample.
![enter image description here](https://hbctraining.github.io/DGE_workshop/img/deseq_counts_overview.png)
With differential expression analysis, we are looking for genes that change in expression between two or more groups (defined in the metadata)

-   case vs. control
-   correlation of expression with some variable or clinical outcome
![enter image description here](https://hbctraining.github.io/DGE_workshop/img/de_variation.png)
The “uninteresting” presents as sources of variation in your data, and so even though the mean expression levels between sample groups may appear to be quite different, it is possible that the difference is not actually significant. **We need to take into account the variation in the data (and where it might be coming from) when determining whether genes are differentially expressed.**
The goal of differential expression analysis is to determine, for each gene, whether the differences in expression (counts) **between groups** (control & case) is significant while also considering the amount of variation observed **within groups** (replicates). To test for significance, we need an appropriate statistical model that accurately performs normalization (to account for differences in sequencing depth, etc.) and variance modeling (to account for few numbers of replicates and large dynamic expression range).
### RNA-seq count distribution
To determine the appropriate statistical model, we need information about the distribution of counts. To get an idea about how RNA-seq counts are distributed, let’s plot the counts for a single sample *Contr_S14_L004.STAR.ReadCounts:*
```R
##Reading the readcount file
raw_read_counts <- read.csv("Contr_S14_L004.STAR.ReadCounts", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
raw_read_counts.df <- data.frame(raw_read_counts)

##Removing the las 5 lines containing the feature counts etc..
raw_read_counts2 <- (raw_read_counts.df[1:(dim(raw_read_counts.df)[1]-5),])
ggplot(raw_read_counts2)+ geom_histogram(aes(x = V2), stat = "bin", bins = 200) +
    xlab("Raw expression counts") +
    ylab("Number of genes")

##If we zoom in close to zero, we can see a large number of genes with counts of zero:
ggplot(raw_read_counts2)+ geom_histogram(aes(x = V2), stat = "bin", bins = 200) +
+     xlab("Raw expression counts") +
+     ylab("Number of genes")+xlim(-5, 500)


```
