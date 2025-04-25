library(SingleCellExperiment)
library(tibble)
library(tidyr)
library(dplyr)
library(igraph)
library(cydar)
library(pdist)
library(reshape2)
library(Signac)
library(edgeR)
library(magrittr)



run_cydar <- function(sce, condition_col="synth_labels",
                      sample_col="synth_samples",
                      reduced.dim="pca.corrected",
                      d=20,
                      batch_col=NULL,
                      alpha=0.1,
                      tol=1.0,
                      downsample=10,
                      returnCd=TRUE){
  ## Make design matrix
  design_df <- as_tibble(colData(sce)[c(sample_col, condition_col, batch_col)]) %>%
    distinct() %>%
    column_to_rownames(sample_col)
  if (is.null(batch_col)) {
    design <- formula(paste('~', condition_col, collapse = ' '))  
  } else {
    design <- formula(paste('~', batch_col, "+", condition_col, collapse = ' '))
  }
  
  ## Make list for each sample
  sample_ls <- split(1:ncol(sce), sce[[sample_col]])
  processed.exprs <- lapply(sample_ls, function(s) reducedDim(sce[,s], reduced.dim)[,1:d])
  cd <- prepareCellData(processed.exprs)
  ## Count cells in hyperspheres
  cd <- cydar::countCells(cd, tol=tol, filter=1, downsample=downsample)
  print("1")  
  # do DA testing with edgeR
  cd.dge <- DGEList(assay(cd), lib.size=cd$totals)
  print("2")
  sim.design <- model.matrix(design, data=design_df)[colnames(cd.dge),]
  sim.dge <- estimateDisp(cd.dge, sim.design)
  sim.fit <- glmQLFit(sim.dge, sim.design)
  sim.res <- glmQLFTest(sim.fit, coef=2)
  print("3")
  # control the spatial FDR
  cydar.res <- sim.res$table
  cydar.res$SpatialFDR <- spatialFDR(intensities(cd), sim.res$table$PValue)
  is.sig <- cydar.res$SpatialFDR <= alpha
  if (returnCd) {
    list(Cd=cd, DAres=cydar.res)
  } else {
    cydar.res
  }
}

cydar2output <- function(sce, cd, da_res, out_type="continuous", alpha=0.1, sample_col="synth_samples", reduced.dim="pca.corrected", d=30){
  nhs <- lapply(cellAssignments(cd), function(hs) as.vector(hs))
  # hs_mat <- sapply(nhs, function(nh) ifelse(1:sum(cd@colData$totals) %in% nh, 1, 0))
  ## Recover cell ids
  ordered.cells <- colnames(cellIntensities(cd))
  hs_mat <- sapply(nhs, function(nh) ifelse(seq_along(ordered.cells) %in% nh, 1, 0))
  rownames(hs_mat) <- ordered.cells
  colnames(hs_mat) <- rownames(da_res)
  if (out_type=="continuous") { 
    da.cell.mat <- hs_mat %*% da_res$logFC
    da.cell <- da.cell.mat[,1]
  } else {
    da.hs <- ifelse(da_res$SpatialFDR < alpha, ifelse(da_res$logFC > 0,"PosLFC" , 'NegLFC'), "NotDA") #选细胞
    da.hs.mat <- sapply(unique(da.hs), function(x) as.numeric(da.hs==x))
    da.cell.mat <- hs_mat %*% da.hs.mat
    da.cell <- apply(da.cell.mat, 1, function(x) colnames(da.cell.mat)[which.max(x)])
  }
  da.cell[colnames(sce)]
}


# main
data <- readRDS("/.../ATACcydar_test.rds")

coldata <- read.csv("/.../CRClabel.csv")
colData(data)$condition <- ifelse(DataFrame(coldata)$synth_labels=="Condition1",0,1)
colData(data)$sample_col <- DataFrame(coldata)$synid

reducedDim(data, "HARMONY_sm") <- reducedDim(data, "HARMONY")[,2:30]

# run the method
cydar_res <- run_cydar(data, condition_col="condition", sample_col="sample_col",
                       reduced.dim = "HARMONY_sm", d=29, tol=2, downsample=5)
out <- cydar2output(data, cydar_res$Cd, cydar_res$DAres, out_type = "label", alpha=0.05)