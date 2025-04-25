library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(Matrix)

peaks.1 <- read.table(
  file = "/.../A002-C-010-D_20200310_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

peaks.2 <- read.table(
  file = "/.../A002-C-010-D_20200702_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

peaks.3 <- read.table(
  file = "/.../A002-C-010-S2_20200310_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

peaks.4 <- read.table(
  file = "/.../A002-C-016-D_20200715_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

peaks.5 <- read.table(
  file = "/.../A002-C-021-D_20200811_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]


peaks.6 <- read.table(
  file = "/.../A002-C-106-D_20200702_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

peaks.7 <- read.table(
  file = "/.../A002-C-114-D_20200811_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

peaks.8 <- read.table(
  file = "/.../A002-C-116-D_20200310_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

peaks.9 <- read.table(
  file = "/.../A002-C-116-S2_20200310_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]


peaks.10 <- read.table(
  file = "/.../A002-C-201-D_08042020_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

peaks.11 <- read.table(
  file = "/.../A002-C-202-D-OCT_20200214_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

peaks.12 <- read.table(
  file = "/.../A002-C-203-D_08042020_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

peaks.13 <- read.table(
  file = "/.../A002-C-204-D_20200702_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

peaks.14 <- read.table(
  file = "/.../A002-C-205-D_20200811_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]


peaks.15 <- read.table(
  file = "/.../CRC-1-8810-D_20200917_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

peaks.16 <- read.table(
  file = "/.../CRC-2-15564-D_20200917_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

peaks.17 <- read.table(
  file = "/.../CRC-3-11773-D_20200917_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

peaks.18 <- read.table(
  file = "/.../CRC-4-8456-D_20200917_peaks.bed",
  col.names = c("chr", "start", "end","width","strand"),
  header=TRUE
)[,c(1,2,3)]

# convert to genomic ranges
gr.1 <- makeGRangesFromDataFrame(peaks.1)
gr.2 <- makeGRangesFromDataFrame(peaks.2)
gr.3 <- makeGRangesFromDataFrame(peaks.3)
gr.4 <- makeGRangesFromDataFrame(peaks.4)
gr.5 <- makeGRangesFromDataFrame(peaks.5)
gr.6 <- makeGRangesFromDataFrame(peaks.6)
gr.7 <- makeGRangesFromDataFrame(peaks.7)
gr.8 <- makeGRangesFromDataFrame(peaks.8)
gr.9 <- makeGRangesFromDataFrame(peaks.9)
gr.10 <- makeGRangesFromDataFrame(peaks.10)
gr.11 <- makeGRangesFromDataFrame(peaks.11)
gr.12 <- makeGRangesFromDataFrame(peaks.12)
gr.13 <- makeGRangesFromDataFrame(peaks.13)
gr.14 <- makeGRangesFromDataFrame(peaks.14)
gr.15 <- makeGRangesFromDataFrame(peaks.15)
gr.16 <- makeGRangesFromDataFrame(peaks.16)
gr.17 <- makeGRangesFromDataFrame(peaks.17)
gr.18 <- makeGRangesFromDataFrame(peaks.18)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.1,gr.2,gr.3,gr.4,gr.5,gr.6,gr.7,gr.8,gr.9,gr.10,gr.11,gr.12,gr.13,gr.14,gr.15,gr.16,gr.17,gr.18))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

# load cell.use
md.1 <- read.table(
  file = "/.../A002-C-010-D_20200310_cells.txt"
) [,1]

md.2 <- read.table(
  file = "/.../A002-C-010-D_20200702_cells.txt"
) [,1]

md.3 <- read.table(
  file = "/.../A002-C-010-S2_20200310_cells.txt"
) [,1]

md.4 <- read.table(
  file = "/.../A002-C-016-D_20200715_cells.txt"
) [,1]

md.5 <- read.table(
  file = "/.../A002-C-021-D_20200811_cells.txt"
) [,1]


md.6 <- read.table(
  file = "/.../A002-C-106-D_20200702_cells.txt"
) [,1]

md.7 <- read.table(
  file = "/.../A002-C-114-D_20200811_cells.txt"
) [,1]

md.8 <- read.table(
  file = "/.../A002-C-116-D_20200310_cells.txt"
) [,1]

md.9 <- read.table(
  file = "/.../A002-C-116-S2_20200310_cells.txt"
) [,1]


md.10 <- read.table(
  file = "/.../A002-C-201-D_08042020_cells.txt"
) [,1]

md.11 <- read.table(
  file = "/.../A002-C-202-D-OCT_20200214_cells.txt"
) [,1]

md.12 <- read.table(
  file = "/.../A002-C-203-D_08042020_cells.txt"
) [,1]

md.13 <- read.table(
  file = "/.../A002-C-204-D_20200702_cells.txt"
) [,1]

md.14 <- read.table(
  file = "/.../A002-C-205-D_20200811_cells.txt"
) [,1]

md.15 <- read.table(
  file = "/.../CRC-1-8810-D_20200917_cells.txt"
) [,1]

md.16 <- read.table(
  file = "/.../CRC-2-15564-D_20200917_cells.txt"
) [,1]

md.17 <- read.table(
  file = "/.../CRC-3-11773-D_20200917_cells.txt"
) [,1]

md.18 <- read.table(
  file = "/.../CRC-4-8456-D_20200917_cells.txt"
) [,1]

# create fragment objects
frags.1 <- CreateFragmentObject(
  path = "/.../A002-C-010-D_20200310_fragments.tsv.gz",
  cells = md.1
)

frags.2 <- CreateFragmentObject(
  path = "/.../A002-C-010-D_20200702_fragments.tsv.gz",
  cells = md.2
)

frags.3 <- CreateFragmentObject(
  path = "/.../A002-C-010-S2_20200310_fragments.tsv.gz",
  cells = md.3
)

frags.4 <- CreateFragmentObject(
  path = "/.../A002-C-016-D_20200715_fragments.tsv.gz",
  cells = md.4
)

frags.5 <- CreateFragmentObject(
  path = "/.../A002-C-021-D_20200811_fragments.tsv.gz",
  cells = md.5
)

frags.6 <- CreateFragmentObject(
  path = "/.../A002-C-106-D_20200702_fragments.tsv.gz",
  cells = md.6
)

frags.7 <- CreateFragmentObject(
  path = "/.../A002-C-114-D_20200811_fragments.tsv.gz",
  cells = md.7
)

frags.8 <- CreateFragmentObject(
  path = "/.../A002-C-116-D_20200310_fragments.tsv.gz",
  cells = md.8
)

frags.9 <- CreateFragmentObject(
  path = "/.../A002-C-116-S2_20200310_fragments.tsv.gz",
  cells = md.9
)


frags.10 <- CreateFragmentObject(
  path = "/.../A002-C-201-D_08042020_fragments.tsv.gz",
  cells = md.10
)

frags.11 <- CreateFragmentObject(
  path = "/.../A002-C-202-D-OCT_20200214_fragments.tsv.gz",
  cells = md.11
)

frags.12 <- CreateFragmentObject(
  path = "/.../A002-C-203-D_08042020_fragments.tsv.gz",
  cells = md.12
)

frags.13 <- CreateFragmentObject(
  path = "/.../A002-C-204-D_20200702_fragments.tsv.gz",
  cells = md.13
)

frags.14 <- CreateFragmentObject(
  path = "/.../A002-C-205-D_20200811_fragments.tsv.gz",
  cells = md.14
)


frags.15 <- CreateFragmentObject(
  path = "/.../CRC-1-8810-D_20200917_fragments.tsv.gz",
  cells = md.15
)

frags.16 <- CreateFragmentObject(
  path = "/.../CRC-2-15564-D_20200917_fragments.tsv.gz",
  cells = md.16
)

frags.17 <- CreateFragmentObject(
  path = "/.../CRC-3-11773-D_20200917_fragments.tsv.gz",
  cells = md.17
)

frags.18 <- CreateFragmentObject(
  path = "/.../CRC-4-8456-D_20200917_fragments.tsv.gz",
  cells = md.18
)

#feature matrix
frags1.counts <- FeatureMatrix(
  fragments = frags.1,
  features = combined.peaks,
  cells = md.1
)

frags2.counts <- FeatureMatrix(
  fragments = frags.2,
  features = combined.peaks,
  cells = md.2
)

frags3.counts <- FeatureMatrix(
  fragments = frags.3,
  features = combined.peaks,
  cells = md.3
)

frags4.counts <- FeatureMatrix(
  fragments = frags.4,
  features = combined.peaks,
  cells = md.4
)

frags5.counts <- FeatureMatrix(
  fragments = frags.5,
  features = combined.peaks,
  cells = md.5
)

frags6.counts <- FeatureMatrix(
  fragments = frags.6,
  features = combined.peaks,
  cells = md.6
)


frags7.counts <- FeatureMatrix(
  fragments = frags.7,
  features = combined.peaks,
  cells = md.7
)

frags8.counts <- FeatureMatrix(
  fragments = frags.8,
  features = combined.peaks,
  cells = md.8
)

frags9.counts <- FeatureMatrix(
  fragments = frags.9,
  features = combined.peaks,
  cells = md.9
)

frags10.counts <- FeatureMatrix(
  fragments = frags.10,
  features = combined.peaks,
  cells = md.10
)

frags11.counts <- FeatureMatrix(
  fragments = frags.11,
  features = combined.peaks,
  cells = md.11
)

frags12.counts <- FeatureMatrix(
  fragments = frags.12,
  features = combined.peaks,
  cells = md.12
)

frags13.counts <- FeatureMatrix(
  fragments = frags.13,
  features = combined.peaks,
  cells = md.13
)

frags14.counts <- FeatureMatrix(
  fragments = frags.14,
  features = combined.peaks,
  cells = md.14
)

frags15.counts <- FeatureMatrix(
  fragments = frags.15,
  features = combined.peaks,
  cells = md.15
)

frags16.counts <- FeatureMatrix(
  fragments = frags.16,
  features = combined.peaks,
  cells = md.16
)

frags17.counts <- FeatureMatrix(
  fragments = frags.17,
  features = combined.peaks,
  cells = md.17
)

frags18.counts <- FeatureMatrix(
  fragments = frags.18,
  features = combined.peaks,
  cells = md.18
)

#Create the objects
obj1 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags1.counts,
    fragments = frags.1
  ), assay = "ATAC")

obj2 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags2.counts,
    fragments = frags.2
  ), assay = "ATAC")

obj3 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags3.counts,
    fragments = frags.3
  ), assay = "ATAC")

obj4 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags4.counts,
    fragments = frags.4
  ), assay = "ATAC")

obj5 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags5.counts,
    fragments = frags.5
  ), assay = "ATAC")

obj6 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags6.counts,
    fragments = frags.6
  ), assay = "ATAC")

obj7 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags7.counts,
    fragments = frags.7
  ), assay = "ATAC")

obj8 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags8.counts,
    fragments = frags.8
  ), assay = "ATAC")

obj9 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags9.counts,
    fragments = frags.9
  ), assay = "ATAC")

obj10 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags10.counts,
    fragments = frags.10
  ), assay = "ATAC")

obj11 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags11.counts,
    fragments = frags.11
  ), assay = "ATAC")

obj12 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags12.counts,
    fragments = frags.12
  ), assay = "ATAC")

obj13 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags13.counts,
    fragments = frags.13
  ), assay = "ATAC")

obj14 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags14.counts,
    fragments = frags.14
  ), assay = "ATAC")

obj15 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags15.counts,
    fragments = frags.15
  ), assay = "ATAC")

obj16 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags16.counts,
    fragments = frags.16
  ), assay = "ATAC")

obj17 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags17.counts,
    fragments = frags.17
  ), assay = "ATAC")

obj18 <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = frags18.counts,
    fragments = frags.18
  ), assay = "ATAC")

# merge all datasets, adding a cell ID to make sure cell names are unique
data <- merge(
  x = obj1,
  y = list(obj2, obj3, obj4, obj5, obj6, obj7, obj8, obj9, obj10, obj11, obj12, obj13, obj14, obj15, obj16, obj17, obj18),
  add.cell.ids = c('A002-C-010-D-R1','A002-C-010-D','A002-C-010-S2','A002-C-016-D','A002-C-021-D','A002-C-106-D','A002-C-114-D','A002-C-116-D','A002-C-116-S2','A002-C-201-D','A002-C-202-D-OCT','A002-C-203-D','A002-C-204-D','A002-C-205-D','CRC-1-8810-D','CRC-2-15564-D','CRC-3-11773-D','CRC-4-8456-D')
)

saveRDS(data,file="/.../combined_seurat.rds")
