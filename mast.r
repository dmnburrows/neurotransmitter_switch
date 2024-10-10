# # library(MAST)
# # library(tidyverse)
# # library(lme4)
# # options(mc.cores=10)
# # print('running')

# # setwd("/cndd3/dburrows/DATA/neurotransmitter_switch/analysis/mast/")
# # cnts <- read.csv("cpm_01_IT-ET_Glut.csv", header=TRUE, check.names=FALSE, row.names=1)
# # cnt_mat <- as.matrix(cnts)
# # dimnames(cnt_mat) <- NULL

# # design <- read.csv("design_01_IT-ET_Glut.csv", header=TRUE)
# # genes <- read.csv("genes.csv", header=FALSE)
# # colnames(design) <- c("wellKey", "condition", 'Sample')
# # colnames(genes) <- c("primerid")

# # mast <- FromMatrix(cnt_mat, design, genes)
# # filtered<-filterLowExpressedGenes(mast,threshold=0.1)

# # #random intercept for sample, fixed for dist
# # zlmint <- zlm(~condition + (1|Sample), sca = filtered, method = 'glmer', ebayes = F, strictConvergence=F, fitArgsD = list(nAGQ = 0))
# # sum_zlmint <- summary(zlmint, logFC=TRUE, doLRT='condition') 
# # write.csv(sum_zlmint$datatable, file="MAST-LRT_sample-intercept.csv")

# # #random intercept and slope for sample, fixed for dist
# # # zlmint_slope <- zlm(~dist_nearest_plaq + (1+dist_nearest_plaq|sample), sca = plq_filtered, method = 'glmer', ebayes = F, strictConvergence=F, fitArgsD = list(nAGQ = 0))
# # # sum_zlmint_slope <- summary(zlmint_slope, logFC=TRUE, doLRT='dist_nearest_plaq') 
# # # write.csv(sum_zlmint_slope$datatable, file="old-APP-cortex_plqdist_MAST-LRT_sample-slopeintercept_reduced.csv")

# library(MAST)
# library(tidyverse)
# library(lme4)
# options(mc.cores=10)

# print('running')

# # Set working directory
# setwd("/cndd3/dburrows/DATA/neurotransmitter_switch/analysis/mast/")

# # Get list of all cpm_ files
# cpm_files <- list.files(pattern = "^cpm_.*\\.csv$")
# # Get list of all design_ files
# design_files <- list.files(pattern = "^design_.*\\.csv$")

# # Loop over all cpm_ files
# for (cpm_file in cpm_files) {
  
#   # Extract the base name (without the prefix and extension)
#   base_name <- sub("^cpm_", "", sub("\\.csv$", "", cpm_file))
  
#   # Match with the corresponding design file
#   design_file <- paste0("design_", base_name, ".csv")
  
#   if (design_file %in% design_files) {
#     print(paste("Processing", cpm_file, "and", design_file))
    
#     # Read the cpm and design files
#     cnts <- read.csv(cpm_file, header = TRUE, check.names = FALSE, row.names = 1)
#     cnt_mat <- as.matrix(cnts)
#     dimnames(cnt_mat) <- NULL
    
#     design <- read.csv(design_file, header = TRUE)
#     colnames(design) <- c("wellKey", "condition", "Sample")
    
#     # Read the genes file and set the column name to "primerid"
#     genes <- read.csv("genes.csv", header = FALSE)
#     colnames(genes) <- c("primerid")
    
#     # Create the SingleCellAssay object for MAST
#     mast <- FromMatrix(cnt_mat, design, genes)
#     filtered <- filterLowExpressedGenes(mast, threshold = 0.1)
    
#     # Run the MAST model with a random intercept for sample
#     zlmint <- zlm(~condition + (1|Sample), sca = filtered, method = 'glmer', ebayes = FALSE, strictConvergence = FALSE, fitArgsD = list(nAGQ = 0))
#     sum_zlmint <- summary(zlmint, logFC = TRUE, doLRT = 'conditionSal') 
    
#     # Save the output uniquely based on the base_name
#     output_file <- paste0("MAST-LRT_", base_name, "_sample-intercept.csv")
#     write.csv(sum_zlmint$datatable, file = output_file)
    
#     print(paste("Completed processing for", cpm_file))
#   } else {
#     print(paste("No matching design file found for", cpm_file))
#   }
# }


library(MAST)
library(tidyverse)
library(lme4)
options(mc.cores=10)
print('running')

setwd("/cndd3/dburrows/DATA/neurotransmitter_switch/analysis/mast/")
plq_cnts <- read.csv("cpm_02_NP-CT-L6b_Glut.csv", header=TRUE, check.names=FALSE, row.names=1)
plq_mat <- as.matrix(plq_cnts)
dimnames(plq_mat) <- NULL

plq_cData <- read.csv("design_02_NP-CT-L6b_Glut.csv", header=TRUE)
colnames(plq_cData) <- c("wellKey", "condition", "Sample")
plq_fData <- read.csv("genes.csv", header=FALSE)
colnames(plq_fData) <- c("primerid")

plq_mast <- FromMatrix(plq_mat, plq_cData, plq_fData)
plq_filtered<-filterLowExpressedGenes(plq_mast,threshold=0.1)

#random intercept for sample, fixed for dist
zlmint <- zlm(~condition + (1|Sample), sca = plq_filtered, method = 'glmer', ebayes = F, strictConvergence=F, fitArgsD = list(nAGQ = 0))
sum_zlmint <- summary(zlmint, logFC=TRUE, doLRT='conditionSal') 
write.csv(sum_zlmint$datatable, file="prac-intercept.csv")

