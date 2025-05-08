#import packages
import numpy as np
import json 
import scanpy as sc
from collections import OrderedDict
import scipy 
import pandas as pd
import ast

#spectra imports 
import Spectra as spc
from Spectra import Spectra_util as spc_tl
from Spectra import K_est as kst
from Spectra import default_gene_sets
import sys
import os
import re
import timeit
import argparse








parser = argparse.ArgumentParser()
parser.add_argument('--data_directory', type=str, help='Path to data directory')

args = parser.parse_args()
data_directory = args.data_directory
print(f'Data directory: {data_directory}')



pattern = r'(\d+)'
matches = re.findall(pattern, data_directory)
num_cells, num_genes, num_facs_ann = map(int, matches[:3])
print(f'num_cells: {num_cells}, num_genes: {num_genes}, Factors: {num_facs_ann}')




adata = sc.read(data_directory)

                
adata.X  = adata.layers['logcounts']
# med_libsize = np.median(np.sum(adata.X, axis = 1))
# print("Median library size is " + str(med_libsize))
# adata_norm = sc.pp.normalize_per_cell(adata, counts_per_cell_after = med_libsize, copy = True)
# adata = sc.pp.log1p(adata_norm, copy = True)

I = adata.uns['Pathways with error']

gene_array = adata.var.index
annotations ={}
annotations['global'] ={'Gene_set_'+str(i):list(gene_array[I[:,i]==1]) for i in range(num_facs_ann)}


start = timeit.default_timer()


model_global = spc.est_spectra(adata = adata, gene_set_dictionary = annotations['global'], # because we do not use the cell types
                                                                                    # we will supply a regular dict
                                                                                    # instead of the nested dict above
                        use_highly_variable = False, cell_type_key =None, #"cell_type_annotations" ,
                        use_weights = True, lam = 0.1, 
                        delta=0.001,kappa = 0.00001, rho = 0.00001, 
                        use_cell_types = False, #set to False to not use the cell type annotations
                        n_top_vals = 25, 
                        clean_gs=True,
                        label_factors=True,
                        num_epochs=10000, #for demonstration purposes we will only run 2 epochs, we recommend 10,000 epochs
                        verbose = True
                       )

stop = timeit.default_timer()

training_time = stop - start

name_of_anndata = data_directory.split('/')[-1]


results = {}

results['factor_loadings'] = adata.uns['SPECTRA_factors']

results['latent_representation'] =  adata.obsm['SPECTRA_cell_scores']
results['factor_importance'] = None

results['training_time'] = training_time


np.save('intermediate_results/spectra_results_'+name_of_anndata.replace('.h5ad','.npy'), results, allow_pickle=True)
