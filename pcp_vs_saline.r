library(DESeq2)
library(tidyverse)

setwd("/cndd3/dburrows/DATA/neurotransmitter_switch/analysis/deseq/")

# Get list of all files that start with 'cnts'
cnts_files <- list.files(pattern = "^cnts.*\\.csv$")

# Read the coldata file once (assuming it's the same for all count files)
coldata <- read.csv("design.csv", row.names = 1)

# Loop over each count file
for (cnts_file in cnts_files) {
  # Read the count data
  cnts <- read.csv(cnts_file, header = TRUE, check.names = FALSE, row.names = 1)
  
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = cnts, 
                                colData = coldata, 
                                design = ~condition)
  
  # Run DESeq2
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "PCP", "Sal"), alpha = 0.1)
  
  # Define the output filename
  output_file <- paste0("outs_", sub("\\.csv$", "", cnts_file), ".csv")
  
  # Save the results
  write.csv(as.data.frame(res), file = output_file)
}


