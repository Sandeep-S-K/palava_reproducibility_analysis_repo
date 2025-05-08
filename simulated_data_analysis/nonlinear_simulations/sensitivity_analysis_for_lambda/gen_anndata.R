
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
path_sigma <- 0.1 # 0.1 in mod splatter 


total_path_genes_prob <- 0.4
de_prob <- 0.4 # all de genes are nonlinear modified splatter, 0.4 in mod splatter

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
    path.sigmaFac = path_sigma,
    de.downProb = 0
  )
    # save unsorted counts in list
  unsorted_counts[[i]] <- counts(sim_paths_list[[i]]) 

  # deteremine library size of the counts
  row_sums <- colSums(unsorted_counts[[i]])

  # normalises based on library size, library size will be re introduced later
  # this is to prevent artifacts afte the concatentation of the each factor
  unsorted_counts[[i]] <- sweep(unsorted_counts[[i]], 2, row_sums, FUN = "/") #counts(sim_paths_list[[i]])/row_sums


  # saves steps of the data
  steps[[i]] <- colData(sim_paths_list[[i]])$Step
    

 # random subsets a fraction of the indecies 1 to 5000
 shuffle_idx <- sample(num_cells, floor(num_cells * 0.2)) 

 # determins the index of ordering based on steps
 order_ind <- order(colData(sim_paths_list[[i]])$Step)

 # shuffle the the randomly selected indecies    
 shuffled_values <- sample(order_ind[shuffle_idx])
 order_ind[shuffle_idx] <- shuffled_values

 sorted_ind[[i]] <- order_ind 
 sorted_counts[[i]] <-  unsorted_counts[[i]][, sorted_ind[[i]]] 
 steps[[i]] <- steps[[i]][sorted_ind[[i]]]

}

# merging counts
sorted_counts_rbinded <- do.call(rbind, sorted_counts)

# introducing library size effects
sorted_counts_rbinded<-sweep(sorted_counts_rbinded, 2, rlnorm(num_cells, 11, 0.2), FUN = "*") 

sce <- SingleCellExperiment(list(counts = sorted_counts_rbinded))


for (i in 1:num_simulations) {
  colData(sce)[paste0('steps_', i)] <- steps[[i]]
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

