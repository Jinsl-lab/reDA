
library(SingleCellExperiment)

library(miloR)
library(tibble)
library(dplyr)
library(igraph)

library(pdist)
library(reshape2)
library(Signac)
library(Seurat)



### SYNTHETIC LABELS ###

.find_centroid <- function(X_emb, cluster_membership){
  cl.ixs <- split(1:nrow(X_emb), cluster_membership)  
  centroid_emb <- sapply(cl.ixs, function(x) colMeans(X_emb[x, , drop=FALSE]))
  centroid_emb
}

.member_weight <- function(x, centroid_dist, m=2){
  # centroid_dist <- pdist(t(x), t(centroid_emb))@dist
  w_memb <- sapply(centroid_dist, function(x) 1/sum(x/centroid_dist)^(2/(m-1)))
}
                   
.scale_to_range <- function(x, min=1, max=10){
  ((x - min(x))/(max(x)-min(x)))*(max-min) + min
}


.logit <- function(x, a=1){
  1/(1+exp(- a * x))
}

add_synthetic_labels_pop <- function(sce, # SingleCellExperiment obj
                                     pop, pop_column="celltype",
                                     pop_enr = 0.7,
                                     redDim="HARMONY", # embedding to use to simulate differential abundance
                                     n_conditions=2, # number of conditions to simulate
                                     n_replicates=3, # number of replicates per condition
                                     n_batches = 2, # number of technical batches per condition (at least 2 replicates per batch)
                                     condition_balance = 1, # the distribution of cells across conditions (1 = equal)
                                     m=2, # Fuzziness parameter (higher m, more fuzziness)
                                     a_logit=0.5, # logit parameter
                                     cap_enr=NULL,
                                     seed=42){
  
  # pop_sce = sce[,sce[[pop_column]]==pop]
  
  set.seed(seed)
  conditions = paste0("Condition", 1:n_conditions)
  
  X_emb = reducedDim(sce, redDim)
  
  ## Find cluster center
  cluster_membership = sce[[pop_column]]
  centroid_emb <- .find_centroid(X_emb, cluster_membership)

  ## Assign weight to each cell for each cluster center
  centroid_dist <- pdist(X_emb, t(centroid_emb))
  centroid_dist <- as.matrix(centroid_dist)
  
  w <- sapply(1:ncol(centroid_dist),  function(j) 
    sapply(1:nrow(centroid_dist), function(i) 
      1/sum(centroid_dist[i,j]/centroid_dist[i,])^(2/(m-1))
    ) 
  )
  colnames(w) <- colnames(centroid_emb)
  rownames(w) <- rownames(X_emb)
  w <- apply(scale(w), 2, .logit, a=a_logit)
  ## Normalize weights from enr_score to 0.5
  enr_scores <- rep(0.5, ncol(w)) ## Generate enrichment prob for each cluster
  # enr_scores <- runif(ncol(w)) ## Generate enrichment prob for each cluster
  names(enr_scores) <- colnames(w)
  if(length(pop_enr) == length(pop)){
    enr_scores[pop] <- pop_enr
  } else{
    # assume all pops have the same enrichment
    pop_enr <- rep(pop_enr, length(pop))
    enr_scores[pop] <- pop_enr
  }
  
  
  # altering the baseline probability can induce a skew towards a condition across _all_ cells
  enr_prob <- sapply(1:ncol(w), function(i) .scale_to_range(w[,i], min=0.5*condition_balance,
                                                            max=enr_scores[i]))
  colnames(enr_prob) <- colnames(centroid_emb)
  
  # need to integrate over these to get the condition probabilities
  # need to set relevant pops only, force the others to ~0.5
  prob_matrix <- enr_prob[,pop]
  if(is(prob_matrix, "matrix")){
    cond_probability <- rowMeans(prob_matrix)
    for(x in seq_along(pop)){
      cond_probability[sce[[pop_column]] == pop[x]] <- prob_matrix[sce[[pop_column]] == pop[x], pop[x]]
    }
  } else{
    cond_probability <- prob_matrix
  }
  
  ## Cap probabilities (to have the same number of DA cells w different maximum Fold Change)
  if (!is.null(cap_enr)) {
    cond_probability <- ifelse(cond_probability > cap_enr, cap_enr, cond_probability)
  }
  # sim3$Condition1_prob <- ifelse(sim3$Condition1_prob > 0.8, 0.8, sim3$Condition1_prob)
  
  cond_probability = cbind(cond_probability, 1 - cond_probability)
  colnames(cond_probability) = conditions
  
  colnames(cond_probability) <- paste0(colnames(cond_probability), "_prob")
  colData(sce)[colnames(cond_probability)] <- cond_probability                   
                     
  # Generate labels for condition and replicates
  middata <- sort(sce$Condition2_prob)
  synth_labels<-ifelse(sce$Condition2_prob < quantile(middata, 0.5*condition_balance), "Condition1", "Condition2")
  #synth_labels <- ifelse(data.frame(cond_probability)$Condition1<0.5,"Condition2","Condition1")
  replicates <- paste0("R", 1:n_replicates)
  batches <- sample(paste0("B", rep(1:n_batches, each=n_replicates)))
  
  synth_samples <- paste0(synth_labels, "_", replicates)
  if(n_batches > 1){
   names(batches) <- sort(unique(synth_samples))
  } else{
    batches <- rep("B1", length(unique(synth_samples)))
    names(batches) <- unique(synth_samples)
  }
  synth_batches <- batches[synth_samples]
  
  # Add synthetic labels and probabilities to colData
  colData(sce)[["synth_labels"]] <- synth_labels
  # colData(sce)[["synth_replicates"]] <- synth_replicates
  colData(sce)[["synth_samples"]] <- synth_samples
  colData(sce)[["synth_batches"]] <- synth_batches
  return(sce)
}

seed <- 42
data_path <- "/.../ATACcydar_test.rds"

pop<-c('1')
pop_enr <- 0.8
pop_col<-"Diagnosis"

sce <- readRDS(data_path)


sce <- add_synthetic_labels_pop(sce, pop=pop, pop_column = pop_col, redDim="LSI", seed=seed, pop_enr=pop_enr,condition_balance = 1,n_replicates=8,n_batches = 1)


data <- sort(sce$Condition2_prob)
pop_tbl <- table(sce[[pop_col]])
por <- sum(pop_tbl[pop]) / (sum(pop_tbl))
quartile1 <- quantile(data, por)
quartile2 <- quantile(data, 1-por)

true_labels <- ifelse(sce$Condition2_prob < quartile1, "NegLFC", ifelse(sce$Condition2_prob > quartile2, "PosLFC", "NotDA"))
colData(sce)[["true_labels"]] <- true_labels

coldata <- data.frame(colData(sce)) %>% rownames_to_column()