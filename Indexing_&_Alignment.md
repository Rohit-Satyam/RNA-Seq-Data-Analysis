## Using RSem

    rsem-prepare-reference --gtf /path to/gencode.v33.annotation.gtf --star --star-path ~/anaconda3/bin/ /path to/GRCh38.primary_assembly.genome.fa /home/parashar/scratch/name_of_output_file 2> rsem.stderr

## Alignments
bsub -e rsem.e -o rsem.o -n 32 -q regularq rsem-calculate-expression --output-genome-bam --strandedness reverse -p 32 --calc-pme --calc-ci --keep-intermediate-files --append-names --sort-bam-by-coordinate --paired-end --estimate-rspd --star --star-path /home/parashar/anaconda3/bin/ /home/parashar/archive/ananda_rnaseq/trimmomatic/paired/Contr_S14_L004_R1_p.fastq /home/parashar/archive/ananda_rnaseq/trimmomatic/paired/Contr_S14_L004_R2_p.fastq /home/parashar/archive/rnaseq/rsem /home/parashar/scratch/last_rsem_test/Contr_S14_L004 2> /home/parashar/scratch/last_rsem_test/Contr_S14_L004.stderr
