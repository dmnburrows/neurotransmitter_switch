#/!/bin/bash

cellranger count --id=cellranger  --transcriptome=/cndd/Public_Datasets/references/mm10_cellranger_H2BmCherry/mm10_H2B-mCherry --fastqs=/cndd3/dburrows/DATA/neurotransmitter_switch/raw/pilot/combined  --sample=combined-pilot --chemistry=SC3Pv3 --include-introns=true --localcores=10

cellranger mkref --genome=mm10_all --fasta=../GRCm38.p3.genome_all.fa --genes=../gencode.vM3.chr_patch_hapl_scaff.annotation.gtf


