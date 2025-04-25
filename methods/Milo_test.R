library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Seurat)

# main
data <- readRDS(".../ATACcydar_test.rds")

coldata <- read.csv("/.../CRClabel.csv")
colData(data)$condition <- ifelse(data.frame(coldata)$synth_labels=="Condition1",0,1)
colData(data)$sample_col <- data.frame(coldata)$synid

reducedDim(data, "HARMONY_sm") <- reducedDim(data, "HARMONY")[,2:30]

traj_milo<-Milo(data)
traj_milo <- buildGraph(traj_milo, k = 10, d = 29,reduced.dim="HARMONY_sm")
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 10, d=29, refined = TRUE,reduced_dims="HARMONY_sm")
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), samples="sample_col")
traj_design <- data.frame(colData(traj_milo))[,c("sample_col", "condition")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$sample_col
## Reorder rownames to match columns of nhoodCounts(milo)
traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop=FALSE]

traj_milo <- calcNhoodDistance(traj_milo, d=29,reduced.dim="HARMONY_sm")
rownames(traj_design) <- traj_design$sample_col
da_results <- testNhoods(traj_milo, design = ~ condition, design.df = traj_design,reduced.dim="HARMONY_sm")

da_results %>%
  arrange(- SpatialFDR) %>%
  head()


traj_milo <- buildNhoodGraph(traj_milo)

#test precision
milo2output <- function(milo, da_res, out_type="continuous", alpha=0.1){
  if (out_type=="continuous") {
    da.cell.mat <- milo@nhoods %*% da_res$logFC
    da.cell <- da.cell.mat[,1]
  } else {
    da.nhoods <- ifelse(da_res$SpatialFDR < alpha, ifelse(da_res$logFC > 0, "PosLFC", 'NegLFC'), "NotDA")
    da.nhoods.mat <- sapply(unique(da.nhoods), function(x) as.numeric(da.nhoods==x))
    da.cell.mat <- milo@nhoods %*% da.nhoods.mat
    da.cell <- apply(da.cell.mat, 1, function(x) colnames(da.cell.mat)[which.max(x)])
  }
  da.cell
}

out <- milo2output(milo, dares, out_type = "label", alpha=0.05)
