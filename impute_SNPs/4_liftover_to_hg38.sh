#!/bin/bash
#BSUB -J liftover_chr22
#BSUB -R "rusage[mem=24000]"
#BSUB -o liftover_chr22.out
#BSUB -e liftover_chr22.err
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 1

module load picard/3.1.1
# module load samtools/1.3

## Get the chain file
#wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
#gunzip hg19ToHg38.over.chain.gz

## Get the reference genome
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#gunzip hg38.fa.gz

## Remove chr prefix from hg38 reference genome
#sed 's/^>chr/>/' hg38.fa > hg38_nochr.fa

## Make the dictionary and index for the hg38 reference genome
#java -jar $PICARDJAR CreateSequenceDictionary R=hg38_nochr.fa O=hg38_nochr.dict
# samtools faidx hg38_nochr.fa

## Run Picard LiftoverVcf
java -jar $PICARDJAR LiftoverVcf \
    I=/home/genevieve.roberts-umw/imputed_genotypes/1KG/test_chr22.dose.vcf \
    O=/home/genevieve.roberts-umw/imputed_genotypes/1KG/test_hg38_chr22.dose.vcf \
    CHAIN=/home/genevieve.roberts-umw/liftover_chain/hg19ToHg38.over.chain \
    REJECT=/home/genevieve.roberts-umw/imputed_genotypes/1KG/test_hg38_rejected_variants.vcf \
    R=/home/genevieve.roberts-umw/liftover_chain/hg38_nochr.fa