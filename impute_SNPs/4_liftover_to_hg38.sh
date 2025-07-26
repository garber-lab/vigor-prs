#!/bin/bash
#BSUB -J liftover_chr22
#BSUB -R "rusage[mem=24000]"
#BSUB -o liftover_chr22.out
#BSUB -e liftover_chr22.err
#BSUB -q short
#BSUB -W 2:00
#BSUB -n 1

module load picard/3.1.1

# Run Picard LiftoverVcf
java -jar $PICARDJAR LiftoverVcf \
    I=/home/genevieve.roberts-umw/imputed_genotypes/1KG/test_chr22.dose.vcf \
    O=/home/genevieve.roberts-umw/imputed_genotypes/1KG/test_hg38_chr22.dose.vcf \
    CHAIN=/home/genevieve.roberts-umw/liftover_chain/hg19ToHg38.over.chain \
    REJECT=/home/genevieve.roberts-umw/imputed_genotypes/1KG/test_hg38_rejected_variants.vcf \
    R=/home/genevieve.roberts-umw/liftover_chain/hg38.fa