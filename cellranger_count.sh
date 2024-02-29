#/!/bin/bash

cellranger count --id=cellranger  --transcriptome=/cndd3/dburrows/DATA/annotations/genome/grcm38.p3/cellranger/mm10_all --fastqs=/cndd3/dburrows/DATA/neurotransmitter_switch/raw/pilot/combined  --sample=combined-pilot --chemistry=SC3Pv3 --include-introns=true --localcores=10



