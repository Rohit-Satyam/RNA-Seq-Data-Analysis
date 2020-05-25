# Sources Referred
1. Harvard DEG tutorial: [here](https://hbctraining.github.io/DGE_workshop/lessons/01_DGE_setup_and_overview.html) 

2. The data we currently use in this tutorial can be found [here](https://github.com/Rohit-Satyam/RNA-Seq-Data-Analysis/tree/master/test)

3. A wonderful RNASeq tutorial is [here](https://github.com/griffithlab/rnaseq_tutorial/wiki) 
 
4. RNASeq Nature Protocol [here](https://www.nature.com/articles/nprot.2013.099)
## Feature Count using HTSeq

    while read p; 
    do 
    name=$(basename ${p} .genome.sorted.bam); 
    bsub -o ./readcounts/${name}_readcount.o -e ./readcounts/${name}_readcount.e -n 8 -q trisutra "htseq-count -f bam -r pos -a 10 -t exon -s reverse -i gene_name -m intersection-strict $p ~/archive/rnaseq/gencode.v33.annotation.gtf > ./readcounts/${p}_ReadCounts.csv 2>./readcounts/${p}.stderr"; 
    done < list




## Differential gene expression analysis overview

So what does this count data actually represent? The count data used for differential expression analysis represents the number of sequence reads that originated from a particular gene. The higher the number of counts, the more reads associated with that gene, and the assumption that there was a higher level of expression of that gene in the sample.
![enter image description here](https://hbctraining.github.io/DGE_workshop/img/deseq_counts_overview.png)
With differential expression analysis, we are looking for genes that change in expression between two or more groups (defined in the metadata)

-   case vs. control
-   correlation of expression with some variable or clinical outcome
![enter image description here](https://hbctraining.github.io/DGE_workshop/img/de_variation.png)
The “uninteresting” presents as sources of variation in your data, and so even though the mean expression levels between sample groups may appear to be quite different, it is possible that the difference is not actually significant. **We need to take into account the variation in the data (and where it might be coming from) when determining whether genes are differentially expressed.**
The goal of differential expression analysis is to determine, for each gene, whether the differences in expression (counts) **between groups** (control & case) is significant while also considering the amount of variation observed **within groups** (replicates). To test for significance, we need an appropriate statistical model that accurately performs normalization (to account for differences in sequencing depth, etc.) and variance modeling (to account for few numbers of replicates and large dynamic expression range).

##  Data exploration and quality assessment

### Data Transformation
For data exploration and visualisation, it is useful to work with transformed versions of the count data. As the count values distribution **is highly skewed**, the log2 transformation helps to approximately **normalize the distributions**. Log base 2 is typically used as it facilitates the conversion back to the original scale: a difference of 1 on the log base 2 scale corresponds to a fold change of 2 on the original count scale. Since count values for a gene can be zero in some conditions (and non-zero in others), we advocates the use of pseudocounts, i.e. transformations of the form

`y = log2 (K + 1) or more generally, y = log2 (K + k0)`);

where K represents the count values and k0 is a positive constant
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
![enter image description here](https://i.ibb.co/PtHZf9J/raw-count-plot.png)

```R
    #Saving the pseudocounts in form of dataframe
    pseudoCount = data.frame(log2(raw_read_counts2$V2 + 1))
    ggplot(pseudoCount, aes(x = log2.raw_read_counts2.V2...1.)) + ylab(expression(log[2](count + 1))) + geom_histogram(colour = "white", fill = "#525252", binwidth = 0.6)+xlab("Control1")
```
![enter image description here](https://i.imgur.com/VcNZtnq.png)

Even more, common statistical methods for exploratory analysis of multidimensional data, especially methods for clustering and ordination (e. g., **principal-component analysis** ), work best for (at least approximately) **homoskedastic data**; this means that the variance of an observable quantity (i.e., here, **the expression strength of a gene**) does not depend on the mean (mean expression value). **In RNA-Seq data, however, variance grows with the mean, with larger variances for larger counts.** For example, if one performs PCA directly on a matrix of “normalized” read counts, the result typically depends only on the few most strongly expressed genes because they show the largest absolute differences between samples. A simple strategy to avoid this is to take the logarithm of the “normalized” count values plus a small constant; however, now the genes with low counts tend to dominate the results because, due to the strong “Poisson” noise inherent to small count values, they show the strongest relative differences between samples. 

In order to make counts approximately homoskedastic, the packages DESeq, DESeq2 and edgeR offers functions to transform the data:

#### How do I know if my data should be modeled using the Poisson distribution or Negative Binomial distribution?

If it’s count data, it should fit the negative binomial, as discussed previously. However, it can be helpful to plot the  _mean versus the variance_  of your data.  _Remember for the Poisson model, mean = variance, but for NB, mean < variance._

Run the following code to plot the  _mean versus variance_  for the control set of samples
```R
## To merge the readcounts files in form of a matrix
DF = do.call(cbind,
             lapply( list.files(pattern=".*.ReadCounts"),
                     FUN=function(x) { 
                         aColumn = read.table(x,header=F);
                         colnames(aColumn)[2] = x;
                         aColumn;
                     }
             )
)
#Removing repeated gene ID columns
DF = DF[,!duplicated(colnames(DF))]

##plot the _mean versus variance_ for the ‘Control and silenced’ replicates
mean_counts <- apply(DF[, 2:4], 1, mean)
variance_counts <- apply(DF[, 3:5], 1, var)
df <- data.frame(mean_counts, variance_counts)
ggplot(df) +
    geom_point(aes(x=mean_counts, y=variance_counts)) + 
    geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
    scale_y_log10() +
    scale_x_log10()
```
![enter image description here](https://i.ibb.co/0CNCS8H/mean-Vs-variance.png)
Note that in the above figure, the variance across replicates tends to be greater than the mean (red line), especially for genes with large mean expression levels. _This is a good indication that our data do not fit the Poisson distribution and we need to account for this increase in variance using the Negative Binomial model (i.e. Poisson will underestimate variability leading to an increase in false positive DE genes).

### Between-sample distribution
It iss useful to contrast the distribution of gene-level expression values on different samples. It can for example be used to display the effects of between-samples before and after filtering and/or normalization.

#### Boxplots
Boxplot method provides an easy way to visualize the distribution of pseoudocounts in each sample. Below is the source code for the same.Note the dots depicts the outliers
```R
##To plot log2 normalized read counts for both control and test
rawcounts <-counts(ddsHTSeq) ##store roaw counts from ddsHTseq
rawcounts <- data.frame(rawcounts) ##saving it in a dataframe format
log2_rawcounts<-data.frame(log2(rawcounts+1)) ## log2 normalization of raw reads
df_boxplot <-melt(log2_rawcounts, variable.names("Samplename")) #variable.names:name of variable used to store values# melting function convert multivariate data into single column, Read here:http://www.datasciencemadesimple.com/melting-casting-r/
df_boxplot <- data.frame(df_boxplot, Condition = substr(df_boxplot$variable, 1, 4)) ## Substring will add a new column and will write string of 1st four characters (cont $ test)
ggplot(df_boxplot, aes(x= variable, y = value, fill = Condition)) + geom_boxplot(alpha = 0.4) + scale_fill_manual(values = c("#619CFF", "#F564E3")) + theme_classic() + xlab("samples") + ylab("log2 converted raw counts")
```
