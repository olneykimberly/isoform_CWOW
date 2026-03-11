.libPaths(c("/tgen_labs/jfryer/kolney/R/rstudio-with_modules-4.4.0-3.sif", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
.libPaths()

#----------------- Libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(Matrix) #, lib.loc = "/usr/local/lib/R/site-library")
library(DESeq2) 
require(openxlsx)
library(edgeR)
library(ggrepel) #, lib.loc = "/usr/local/lib/R/site-library")
library(glmGamPoi)
library(devtools)
library(reshape2)
library(edgeR)  
library(limma)  
library(tximport)
library(tidyverse)
library(GenomicFeatures)
library(data.table)
library(gplots)
library(variancePartition)
library(NatParksPalettes) # colors
library(DRIMSeq)
library(ggpubr)
library(dittoSeq)
library(gridExtra)

#----------------- Define variables
path_ref <- c("/tgen_labs/jfryer/projects/references/human/GRCh38")

#----------------- Data
metadata_n580 <- read.delim("/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/metadata/metadata_n580.tsv", header = TRUE, sep = "\t")
metadata_kbase <- read.delim("/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/metadata/kbase_and_metadata.tsv", header = FALSE, sep = "\t")
colnames(metadata_kbase) <- c("NPID", "CWOW_NPID", "SMRT_cell", "SMRT_batch", "filename", "group", "sex_chr")
metadata <- merge(metadata_kbase, metadata_n580, by = "NPID")
rm(metadata_n580, metadata_kbase)

saveToPDF <- function(...) {
  d = dev.copy(pdf, ...)
  dev.off(d)
}


LBD <- "LBD"
AD <- "AD"
PA <- "PA"
CONTROL <- "CONTROL"
control_color <- "#4682B4" 
AD_color <- "#B4464B" 
PA_color <- "#B4AF46" 
LBD_color <- "gray35" 
control_shape <- c(15) # square
AD_shape <- c(16) # circle
PA_shape <- c(17) # triangle
LBD_shape <- c(18) # diamond

TypeColors <- c("#4682B4", "#B4AF46","#B4464B", "gray35")
ATSColors <- c("#4682B4", "#B4AF46","#B4464B", "gray35", "gray65", "gray", "gray85")
colorbindColors <- dittoColors()

SexColors <- c("#490092", "#D55E00")
correlationColors <-
  colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))
