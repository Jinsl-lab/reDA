library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)

library(future)
plan(strategy = "multisession")
options(future.globals.maxSize = 100000 * 1024^5)

data<-readRDS("/.../cna_milo_test_final_seurat.rds")

pfm <- getMatrixSet(
x = JASPAR2020,
opts = list(collection = "CORE", species = "Homo sapiens")
)

data <- AddMotifs(
object = data,
genome = BSgenome.Hsapiens.UCSC.hg38,
pfm = pfm
)

Idents(data) <- 'CellType'

da_peaks <- FindMarkers(
  object = data,
  ident.1 = 'Cancer Associated Fibroblasts',
  only.pos = TRUE,
  test.use = "LR",
  min.pct = 0.3,
  latent.vars = 'nCount_ATAC'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])

enriched.motifs <- FindMotifs(
  object = data,
  features = top.da.peak
)

