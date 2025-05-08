library(anndata)
library(slalom)
library(SingleCellExperiment)
library(GSEABase)
library(Dict)
library(stringr)

# function to introduce zero cols or rows to matrix

insert_zeros <- function(mat, indices, axis = "row") {
  if (axis == "row") {
    # Insert zero rows
    for (idx in indices) {
      mat <- rbind(mat[1:(idx - 1),, drop = FALSE], 
                   matrix(0, nrow = 1, ncol = ncol(mat)), 
                   mat[idx:nrow(mat),, drop = FALSE])
    }
  } else if (axis == "col") {
    # Insert zero columns
    for (idx in indices) {
      mat <- cbind(mat[,1:(idx - 1), drop = FALSE], 
                   matrix(0, nrow = nrow(mat), ncol = 1), 
                   mat[,idx:ncol(mat), drop = FALSE])
    }
  } else {
    stop("Invalid axis. Please specify 'row' or 'col'.")
  }
  return(mat)
}




args <- commandArgs(trailingOnly = TRUE)

data_dr  <- args[1]

adata <- read_h5ad(data_dr)

sce <-zellkonverter::AnnData2SCE(adata)

num_cells<-length(colnames(sce))

# Indtroduces a non zero element (1) in first cell  when the logcount values for a gene are zero
# if not slalom gives the error **Some genes have zero variance in expression**
zero_genes_ind <- which(rowSums(logcounts(sce)==0)==num_cells)

cat('Gene ind with all zeros',zero_genes_ind)

genes_to_remove <- rownames(sce)[zero_genes_ind]

#removing those genes with all zero elements 
sce <- sce[! rownames(sce) %in% genes_to_remove,]




# Load the library for regular expressions
library(stringr)

# Define the string
string <- data_dr

# Define the regular expression pattern
pattern <- "(?<=_formatted_)(.*?)(?=\\..*?$)"

# Extract the error_str using regex
error_str <- str_match(string, pattern)[,2]


# for the case where no genes are removed 
if (length(zero_genes_ind) > 0) {
  gene_sets_array <- metadata(sce)$Pathways.with.error[-zero_genes_ind,]
} else {
  gene_sets_array <- metadata(sce)$Pathways.with.error
}

gene_sets_list <- list()

# Loop through the columns of gene_sets_array
for (i in 1:ncol(gene_sets_array)) {
    # Create GeneSet object and add it to the list
    gene_set_name <- paste("Fac_", i, sep="")
    gene_set_gene_ids <- rownames(sce)[gene_sets_array[, i] == 1]
    gene_sets_list[[gene_set_name]] <- GeneSet(setName = gene_set_name, geneIds = gene_set_gene_ids)
}
genesets<-GeneSetCollection(gene_sets_list)



sce_slalom <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = logcounts(sce))
)
num_ann <- length(genesets)
num_unann <- 1

model <- newSlalomModel(sce_slalom, genesets, n_hidden = num_unann, min_genes = 1, prune_genes=F)
model <- initSlalom(model)

start_time <- Sys.time()

model <- trainSlalom(model, nIterations = 10000)

end_time <- Sys.time()

training_time <-as.numeric(end_time - start_time)
print(paste('training_time ', training_time))
terms <-topTerms(model, unannotated_dense =T, mad_filter = 0)
terms <- terms[order(terms$term), ]
print(terms)

factor_loadings<- model$W_E1
#factor_loadings<-factor_loadings[,c(3:6,1:2)]

factor_loadings <- factor_loadings[,c((num_unann+1):(num_ann+num_unann),1:num_unann)]

latent <- model$X_E1[,c((num_unann+1):(num_ann+num_unann),1:num_unann)]
i_posterior <- model$W_gamma0[,c((num_unann+1):(num_ann+num_unann),1:num_unann)]


# re-introduce gene that were removed due ot being all zeros
factor_loadings<-insert_zeros(factor_loadings, zero_genes_ind, "row")
i_posterior<-insert_zeros(i_posterior, zero_genes_ind, "row")


save_dr<-sub('data/', "slalom_intermediate_results_R/slalom_results_", data_dr)
save_dr<-sub('.h5ad', "", save_dr)

write.csv(latent, paste0(save_dr, '_latent.csv'), row.names = FALSE)
write.csv(factor_loadings, paste0(save_dr, '_factor_loadings.csv'), row.names = FALSE)
write.csv(terms$relevance, paste0(save_dr, '_factor_importance.csv'), row.names = FALSE)
write.csv(training_time, paste0(save_dr, '_training_time.csv'), row.names = FALSE)
write.csv(i_posterior, paste0(save_dr,'_posterior.csv'), row.names = FALSE)
write.csv(model$converged, paste0(save_dr, '_converged.csv'), row.names = FALSE)



