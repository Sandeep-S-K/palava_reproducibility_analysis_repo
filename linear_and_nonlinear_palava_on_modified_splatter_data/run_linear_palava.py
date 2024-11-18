import os

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import scvi
import subprocess
import argparse
#import seaborn as sb
from scvi import settings
import re
import json
import torch
import time
import timeit

num_facs_unann = 1

parser = argparse.ArgumentParser()
parser.add_argument('--data_directory', type=str, help='Path to data directory')


args = parser.parse_args()

# Access the data directory
data_directory = args.data_directory
print(f'Data directory: {data_directory}')



# Determine num cells, genes, factors from file name
pattern = r'(\d+)'
matches = re.findall(pattern, data_directory)
num_cells, num_genes, num_facs_ann = map(int, matches[:3])
print(f'num_cells: {num_cells}, num_genes: {num_genes}, Factors: {num_facs_ann}')



adata = sc.read(data_directory)


gene_names = adata.var



pathways = adata.uns['Pathways with error']

path_un = torch.zeros(num_genes, dtype=torch.float64)

pathways_bool = [pathways[:,i] for i in range(pathways.shape[1])] 


pathways_bool += [path_un]*num_facs_unann

mask = torch.ones(num_genes) #(ground_truth_pathways[:, 3]==0).int()
masks =[mask for _ in pathways_bool]


settings.seed = 0

plan_kwargs = {'lr' : 5e-4,'n_epochs_kl_warmup' : 400,'n_steps_kl_warmup' : None, 'max_kl_weight' : 1.0}

SCVI_palava = scvi.model.SCVI_palava

SCVI_palava.setup_anndata(adata)

lam = [0.2] * num_facs_ann + [0.2] * num_facs_unann

nonlinear  = [0 for _ in range(len(pathways_bool))]

for i in range(num_facs_ann):
    nonlinear[i] = 0



scvi_palava = SCVI_palava(adata, masks= masks, n_latent =  len(pathways_bool), nonlinear = nonlinear, other_n_latent = 20, lam = lam, pathways_bool = pathways_bool, dispersion = 'gene', palava_n_hidden = 0, palava_n_layers = 0,  n_layers = 3, use_batch_norm = 'none',use_layer_norm = 'none', momentum_train = 0.8,  weight_decay_train = 0.0001, decoder_bias= False,  non_neg_decoder =True)

start = timeit.default_timer()

scvi_palava.train(1000, plan_kwargs = plan_kwargs, batch_size = 100, gradient_clip_val=100.0)

stop = timeit.default_timer()
training_time = stop - start
print('Time: ', training_time)  

slopes_for_all_fac = scvi_palava.spline_approach()


factor_loadings = np.vstack([np.quantile(abs(slopes_for_all_fac[fac]), q = 0.95, axis = 0) for fac in range(num_facs_unann + num_facs_ann)])

name_of_anndata = data_directory.split('/')[-1]


results = {}

results['factor_loadings'] = factor_loadings

latent = scvi_palava.get_latent_representation()
results['latent_representation'] =  latent[:,0:(num_facs_ann + num_facs_unann)]
results['latent_representation_with_other'] = latent
results['factor_importance'] = scvi_palava.factor_importance(adata, num_batchs =500)

results['training_time'] = training_time
results['spline_slopes'] = np.array(slopes_for_all_fac)

#scvi_palava.save('intermediate_results/linear_palava_no_bias_results_'+name_of_anndata.replace('.h5ad','_') + 'model/scvi_model' )



np.save('intermediate_results/linear_palava_results_'+name_of_anndata.replace('.h5ad','.npy'), results, allow_pickle=True)
