prc<- function(inpath, dimsta, dimend, group, assay="ATAC", hadpath, issave="NO", rdspath=NULL){
  library(Signac)
  library(Seurat)
  library(harmony)
  library(readr)
  library(GenomicRanges)
  library(SeuratDisk)

  x <- readRDS(inpath)
  x <- RunTFIDF(x)
  x <- FindTopFeatures(x, min.cutoff = 20)
  x <- RunSVD(x)
  #without harmony
  x <- RunUMAP(x, dims = dimsta:dimend, reduction = 'lsi' ,reduction.name = "umap_naive")
  x <- RunHarmony(object=x, group.by.vars = group, reduction = 'lsi', max.iter.harmony = 100, assay.use = assay, project.dim = FALSE)
  x <- RunUMAP(x, dims = dimsta:dimend, reduction = 'harmony',reduction.name = "umap")
  x <- FindNeighbors(object = x, reduction = 'harmony', dims = dimsta:dimend)
  SeuratDisk::SaveH5Seurat(x, filename = hadpath)
  SeuratDisk::Convert(hadpath,dest = "h5ad")       

  if(issave=="Yes"){
    saveRDS(x,file=rdspath)
  }
}



