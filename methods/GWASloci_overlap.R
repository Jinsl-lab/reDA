########download loci
library(gwasrapidd)

studies <- get_studies(efo_trait = 'colorectal cancer')

my_associations <- get_associations(study_id = studies@studies$study_id)

# Get association ids for which pvalue is less than 1e-6.
dplyr::filter(my_associations@associations, pvalue < 1e-6) %>% # Filter by p-value
  tidyr::drop_na(pvalue) %>%
  dplyr::pull(association_id) -> association_ids # Extract column association_id

my_associations2 <- my_associations[association_ids]

var<-get_variants(association_id=my_associations2@associations$association_id)


########overlap test
library(GenomicRanges)
library(tidyverse)

library(tidyr)


pos.peak<-readRDS("/.../CAFpeaks.rds")

top.pos.peak <- pos.peak[pos.peak$p_val < 0.005 & pos.peak$pct.1 > 0.2, ]
dars<-top.pos.peak
dars<-data.frame(text = rownames(dars))

dars <- dars %>%
  separate(text, into = c("chr", "start","end"), sep = "-")

dars$start<-as.numeric(dars$start)-1000
dars$end<-as.numeric(dars$end)+1000
dars$info<-as.numeric(dars$end)-as.numeric(dars$start)+1

colnames(dars) <- c("chr", "start", "end", "info") 

# load GWAS data

gwas <- read.table("/.../CRCloci.txt",skip = 1)
colnames(gwas) <- c("row","variant_id", "merged", "functional_class", "chromosome_name", "chromosome_position", "chromosome_region", "last_update_date")

gwas <- gwas[!is.na(gwas$chromosome_name), ]

gwas$start<-gwas$chromosome_position
gwas$end<-gwas$chromosome_position
gwas$chr<-paste("chr",gwas$chromosome_name, sep = "")
gwas$info<-1

gwas<-gwas[,c("chr", "start", "end", "info","variant_id")]





# convert to GenomicRanges
dars_gr <- makeGRangesFromDataFrame(dars, keep.extra.columns = TRUE)
gwas_gr <- makeGRangesFromDataFrame(gwas, keep.extra.columns = TRUE)

# analysis
overlap <- findOverlaps(dars_gr, gwas_gr)


overlap_dars <- dars_gr[queryHits(overlap)]
overlap_gwas <- gwas_gr[subjectHits(overlap)]

# combine
overlap_df <- data.frame(
  dars_chr = seqnames(overlap_dars),
  dars_start = start(overlap_dars),
  dars_end = end(overlap_dars),
  dars_info = overlap_dars$info,
  gwas_chr = seqnames(overlap_gwas),
  gwas_start = start(overlap_gwas),
  gwas_end = end(overlap_gwas),
  gwas_info = overlap_gwas$info,
  gwas_variant = overlap_gwas$variant_id 
)

