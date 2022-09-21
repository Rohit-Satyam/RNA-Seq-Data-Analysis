# RNA-Seq-Data-Analysis
Blogpost by Crazyhottomy on htseq: http://crazyhottommy.blogspot.com/2013/10/

> I now use the following commands to get SRR IDs

```
pysradb gse-to-srp GSE149118
pysradb srp-to-srr SRP257905 > samples.txt
awk -F'\t' '{print $2}' samples.txt | tail +2 | parallel -j 3 "fasterq-dump  --outdir . -e 10 --skip-technical --split-3 {}"
```

## RNA Seq Best Practices Highlights
![enter image description here](https://i.imgur.com/HZJyWvL.png)
![enter image description here](https://media.springernature.com/full/springer-static/image/art:10.1186/s13059-016-0881-8/MediaObjects/13059_2016_881_Fig1_HTML.gif?as=webp)

## Experimental design
1. Remove the highly abundant **ribosomal RNA (rRNA)**, which typically constitutes over 90 % of total RNA in the cell, leaving the 1–2 % comprising messenger RNA (mRNA) that we are normally interested in. For eukaryotes, this involves choosing whether to enrich for mRNA using poly(A) selection or to deplete rRNA. Poly(A) selection typically requires a relatively **high proportion of mRNA with minimal degradation** as **measured by RNA integrity number (RIN**), which normally yields a higher overall fraction of reads falling onto known exons. Many biologically relevant samples (such as tissue biopsies) cannot, however, be obtained in great enough quantity or good enough mRNA integrity to produce good poly(A) RNA-seq libraries and therefore require ribosomal depletion. For bacterial samples, in which mRNA is not polyadenylated, the only viable alternative is ribosomal depletion.
2. Sequencing can involve single-end (SE) or paired-end (PE) reads, although the latter is preferable for **de novo transcript discovery or isoform expression analysis** [[4](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8#ref-CR4 "Katz Y, Wang ET, Airoldi EM, Burge CB. Analysis and design of RNA sequencing experiments for identifying isoform regulation. Nat Methods. 2010;7:1009–15."), [5](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8#ref-CR5 "Garber M, Grabherr MG, Guttman M, Trapnell C. Computational methods for transcriptome annotation and quantification using RNA-seq. Nat Methods. 2011;8:469–77.")]. Similarly, longer reads improve mappability and transcript identification [[5](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8#ref-CR5 "Garber M, Grabherr MG, Guttman M, Trapnell C. Computational methods for transcriptome annotation and quantification using RNA-seq. Nat Methods. 2011;8:469–77."), [6](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8#ref-CR6 "Łabaj PP, Leparc GG, Linggi BE, Markillie LM, Wiley HS, Kreil DP. Characterization and improvement of RNA-Seq precision in quantitative transcript expression profiling. Bioinformatics. 2011;27:i383–91.")]. The best sequencing option depends on the analysis goals. The cheaper, short SE reads are normally sufficient for studies of gene expression levels in well-annotated organisms, whereas longer and PE reads are preferable to characterize poorly annotated transcriptomes.
3. **Depth:** While some authors will argue that as few as **5 million mapped reads** are sufficient to quantify accurately medium to highly expressed genes in most eukaryotic transcriptomes, others will sequence up to **100 million reads to quantify precisely genes and transcripts that have low expression levels** [[7](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8#ref-CR7 "Sims D, Sudbery I, Ilott NE, Heger A, Ponting CP. Sequencing depth and coverage: key considerations in genomic analyses. Nat Rev Genet. 2014;15:121–32.")]. When studying **single cells**, which have limited sample complexity, quantification is often carried out with **just 1 million reads** but may be done reliably for highly expressed genes with as few **as 50,000 reads** [[8](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8#ref-CR8 "Pollen AA, Nowakowski TJ, Shuga J, Wang X, Leyrat AA, Lui JH, et al. Low-coverage single-cell mRNA sequencing reveals cellular heterogeneity and activated signaling pathways in developing cerebral cortex. Nat Biotechnol. 2014;32:1053–8.")]; even 20,000 reads have been used to differentiate cell types in splenic tissue [[9](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8#ref-CR9 "Jaitin DA, Kenigsberg E, Keren-Shaul H, Elefant N, Paul F, Zaretsky I, et al. Massively parallel single-cell RNA-seq for marker-free decomposition of tissues into cell types. Science. 2014;343:776–9.")]. Moreover, optimal library size depends on the complexity of the targeted transcriptome. Experimental results suggest that deep sequencing improves quantification and identification but might also result in the detection of transcriptional noise and off-target transcripts [[10](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8#ref-CR10 "Tarazona S, Garcia-Alcalde F, Dopazo J, Ferrer A, Conesa A. Differential expression in RNA-seq: a matter of depth. Genome Res. 2011;21:2213–23.")]. Saturation curves can be used to assess the improvement in transcriptome coverage to be expected at a given sequencing depth [[10](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8#ref-CR10 "Tarazona S, Garcia-Alcalde F, Dopazo J, Ferrer A, Conesa A. Differential expression in RNA-seq: a matter of depth. Genome Res. 2011;21:2213–23.")].
4. **Number of replicates:** Depends on three factors:  1)  **Technical Variability:** variability in the measurements, which is influenced by the technical noise and the biological variation. Occurs at the level of RNA extraction and library preparation and may introduce biases. **2) Biological Variability**: biological replication is required if inference on the population is to be made, with **three replicates being the minimum for any inferential analysis.** For a proper statistical power analysis, estimates of the within-group variance and gene expression levels are required. This information is typically not available beforehand but can be obtained from similar experiments. The **exact power will depend on the method used for differential expression analysis**, and software packages exist that provide a theoretical estimate of power over a range of variables, given the within-group variance of the samples, which is intrinsic to the experiment.

> Filtering out genes that are expressed at low levels prior to differential expression analysis reduces the severity of the correction and may improve the power of detection. **Increasing sequencing depth** also can **improve statistical power** for lowly expressed genes, and for any given sample there exists a level of sequencing at which power improvement is best achieved **by increasing the number of replicates**. Tools such as **Scotty** are available to calculate the best trade-off between sequencing depth and replicate number given some budgetary constraints.
5. If reads primarily accumulate at the 3’ end of transcripts in poly(A)-selected samples, this might indicate low RNA quality in the starting material. The GC content of mapped reads may reveal PCR biases. Use [RSeQC](http://rseqc.sourceforge.net/) 
6. we expect between 70 and 90 % of regular RNA-seq reads to map onto the human genome (depending on the read mapper used), with a significant fraction of reads mapping to a limited number of identical regions equally well (‘multi-mapping reads’). **When reads are mapped against the transcriptome, we expect slightly lower total mapping percentages** because reads coming from unannotated transcripts will be lost, and significantly more multi-mapping reads because of reads falling onto exons that are shared by different transcript isoforms of the same gene. Genomic multireads are primarily due to repetitive sequences or shared domains of paralogous genes. They normally account for a significant fraction of the mapping output when mapped onto the genome and should not be discarded



![enter image description here](https://upload.wikimedia.org/wikipedia/commons/thumb/5/54/Gene_structure_eukaryote_2_annotated.svg/800px-Gene_structure_eukaryote_2_annotated.svg.png)

## RSeqQC Tutorial
Create conda environment using python 2.7
Install RSeqQC using `pip install RSeQC`
Check installation success via typing `bam2fq.py --version`

To use RseQC you need a GFF file in bed format
Use the gtf2bed.pl script as follows: `./gtf2bed.pl file.gff > file_gff.bed`

RNASeq Nature Protocol [here](https://www.nature.com/articles/nprot.2013.099)
