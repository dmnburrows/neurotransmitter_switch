library(DESeq2)
library(tidyverse)

setwd("/cndd3/dburrows/DATA/neurotransmitter_switch/analysis/deseq/")

coldata <- read.csv("design.csv", row.names=1)
cnts <- read.csv("cnts.csv", header=TRUE, check.names=FALSE, row.names=1)

dds <-DESeqDataSetFromMatrix(countData=cnts, 
                       colData=coldata, 
                       design=~condition) 



dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "PCP", "Sal"), alpha = 0.1)
write.csv(as.data.frame(res), 
          file="outs.csv")


