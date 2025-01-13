
suppressPackageStartupMessages({
  library(splatter)
  library(scater)
  library(reticulate)
})

#reticulate::use_python(python = '/data/gpfs/projects/punim0614/sandeep_sk/Python_HPC/anaconda3/envs/r-env/bin/python',required = T)

args <- commandArgs(trailingOnly = TRUE)

num_cells <- args[1]
num_genes <- args[2]
num_simulations <- args[3] 
seed <- args[4] 

num_cells<-as.integer(num_cells)
num_genes <- as.integer(num_genes)
num_simulations<- as.integer(num_simulations)
seed <- as.integer(args[4]) 

set.seed(seed)




de_facLoc <- 0.4
path_sigma <- 0.1


total_path_genes_prob <- 0.4
de_prob <- 0.4

path_nonlinearprob <- (total_path_genes_prob - de_prob)/(1 - de_prob)
print(path_nonlinearprob)

num_genes_per_fac <-  floor(num_genes/num_simulations)
sim_paths_list <- list()
steps <- list()
unsorted_counts <- list()
sorted_ind <- list()
sorted_counts <- list()

for (i in 1:num_simulations) {
  sim_paths_list[[i]] <- splatSimulate(
    de.prob = de_prob,
    de.facLoc = de_facLoc,
    batchCells = num_cells,
    nGenes = num_genes_per_fac,
    method = "paths",
    verbose = FALSE,
    path.nonlinearProb = path_nonlinearprob,
    path.sigmaFac = path_sigma
  )
  unsorted_counts[[i]] <- counts(sim_paths_list[[i]])
  steps[[i]] <- colData(sim_paths_list[[i]])$Step
  
 row_sums <- colSums(unsorted_counts[[i]])
 # Sort the columns based on the increasing order of row sums
 sorted_ind[[i]] <- order(row_sums) 
 sorted_counts[[i]] <- unsorted_counts[[i]][, sorted_ind[[i]]]

}

sorted_counts_rbinded <- do.call(rbind, sorted_counts)
sce <- SingleCellExperiment(list(counts = sorted_counts_rbinded))


for (i in 1:num_simulations) {
  colData(sce)[paste0('steps_', i)] <- steps[[i]][sorted_ind[[i]]]
}

for (i in 1:num_simulations) {
  col_name <- paste0('DEFac_', i)
  prefix_length <- sum((i - 1) * num_genes_per_fac)  # Adjust based on the pattern
  # Set values for rowData
  rowData(sce)[col_name] <- c(
    rep(1, prefix_length),
    rowData(sim_paths_list[[i]])$DEFacPath1 + rowData(sim_paths_list[[i]])$SigmaFacPath1,
    rep(1, num_genes_per_fac *(num_simulations-1) - prefix_length))
}

rownames(sce) <- paste0("Gene", 1:(num_genes))

sce <- logNormCounts(sce)

sce <- runPCA(sce)

adata<-zellkonverter::SCE2AnnData(sce)
anndata::write_h5ad(adata,paste0('data/sim_anndata_', num_cells,'_cells_',num_genes,'_genes_',num_simulations,'_factors_',seed,'_seed.h5ad'))


