library(DAseq)
library(Matrix)
library(Signac)
library(Seurat)
library(scales)
library(ggplot2)
library(tidyr)
library(dplyr)

X<-readRDS("/.../cna_milo_test_final_seurat.rds")
coldata <- read.csv("/.../CRClabel.csv")
X@meta.data$condition <- ifelse(data.frame(coldata)$synth_labels=="Condition1",0,1)
X@meta.data$sample_col <- data.frame(coldata)$synth_samples

Xlabel<-data.frame(label=c("Condition1_R1","Condition1_R2","Condition1_R3","Condition1_R4","Condition1_R5","Condition1_R6","Condition1_R7","Condition1_R8","Condition2_R1","Condition2_R2","Condition2_R3","Condition2_R4","Condition2_R5","Condition2_R6","Condition2_R7","Condition2_R8"),condition=c('NR','NR','NR','NR','NR','NR','NR','NR',"R","R","R","R","R","R","R","R"))

labels_res <- Xlabel[Xlabel$condition == "R", "label"]   #small score from R
labels_nonres <- Xlabel[Xlabel$condition == "NR", "label"] #big score from NR

X.harmony<-Embeddings(X@reductions[["harmony"]])[,2:30]
X2d<-Embeddings(X@reductions[["umap"]])
da_cells <- getDAcells(
  X = X.harmony,
  cell.labels = X@meta.data$sample_col,  
  labels.1 = labels_res,
  labels.2 = labels_nonres,
  k.vector = seq(50, 500, 50),
  plot.embedding = X2d
)

