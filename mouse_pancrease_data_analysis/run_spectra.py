#import packages
import numpy as np
import scipy

#spectra imports 
import Spectra as spc
from Spectra import Spectra_util as spc_tl
from Spectra import K_est as kst
import scanpy as sc

import timeit
import argparse

#pathways = [torch.tensor(i) for i in adata.uns['pathways_9000hvg']]


parser = argparse.ArgumentParser()
parser.add_argument('--data_directory', type=str, help='Path to data directory')

args = parser.parse_args()
data_directory = args.data_directory
print(f'Data directory: {data_directory}')


adata = sc.read(data_directory)


pathway_names = adata.uns['pathway_names'].tolist()
print(pathway_names)
num_facs_ann = len(pathway_names)

adata.X  = adata.X.toarray()       
adata.X  = adata.X * 1.0

med_libsize = np.median(np.sum(adata.X, axis = 1))
print("Median library size is " + str(med_libsize))
adata_norm = sc.pp.normalize_per_cell(adata, counts_per_cell_after = med_libsize, copy = True)
adata = sc.pp.log1p(adata_norm, copy = True)

I = adata.uns['pathways_5000hvg'].T #adata.uns['Pathways with error']
print(I, I.shape)
gene_array = adata.var.index
annotations ={}
annotations['global'] ={pathway_names[i]:list(gene_array[I[:,i]==1]) for i in range(num_facs_ann)}


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
                        num_epochs=10000,  #10000, #for demonstration purposes we will only run 2 epochs, we recommend 10,000 epochs
                        L = 93
                       )

stop = timeit.default_timer()

training_time = stop - start

name_of_anndata = data_directory.split('/')[-1]





results = {}

results['factor_loadings'] = adata.uns['SPECTRA_factors']

print('#' * 100)
print('results[factor_loadings]', results['factor_loadings'].shape)

results['latent_representation'] =  adata.obsm['SPECTRA_cell_scores']
results['factor_importance'] = None

results['training_time'] = training_time

adata.write('intermediate_results/spectra_adata'+name_of_anndata)

np.save('intermediate_results/spectra_results_'+name_of_anndata.replace('.h5ad','.npy'), results, allow_pickle=True)

model_global.save('intermediate_results/spectra_results_'+name_of_anndata.replace('.h5ad','') +'_model.pt')
