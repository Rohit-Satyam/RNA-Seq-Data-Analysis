#!/bin/sh

## The following script uses conda environment and install entrez utilities if absent. If you don't work with conda comment out the conda command ##

##................................................................###
## To check if entrez-direct is installed ##

if ! type "efetch" > /dev/null; then            ## check if 
  echo "E-utilities Not Found."
  echo "installing entrez-direct from conda"
  conda install -c bioconda entrez-direct
fi
##................................................................##

##To input File

input=$1

##................................................................##

## To convert GSM number to SRR number ##

for GSM in `cat $input`
do
echo "Retriving $GSM from NCBI GEO";
SRR=`esearch -db sra -query $GSM | efetch -format docsum | xtract -pattern DocumentSummary -element Run@acc`;
fastq-dump --outdir fastq --gzip --skip-technical --origfmt --readids --read-filter pass --dumpbase --split-3 --clip -A $SRR
echo 'retrival_successful_FASTQ DUmped';
done
