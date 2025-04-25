library(simATAC)
library(rhdf5)
library(Matrix)
library(SeuratDisk)
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(Signac)
library(tidyr)

#initialization
object <- newsimATACCount()

count <- getCountFromh5("/.../GSM5766464_h3_act2_sort.snap")



object <- simATACEstimate(t(count))

object1 <- setParameters(object, noise.mean = -0.3, noise.sd = 0.3, seed=55205)
#simulation
##\item{\code{[species]}}{An string indicating the species of the input cells. 
##simATAC supports "hg38", "hg19", "mm9", and "mm10" in the current version.}

sim1 <- simATACSimulate(object1, nCells = 920,species='hg19')   #choose cell number
sim1.seurat<-CreateSeuratObject(counts = CreateChromatinAssay(counts=counts(sim1),ranges=StringToGRanges(regions=rownames(counts(sim1)),sep=c(":","-"))),assay = "ATAC")

