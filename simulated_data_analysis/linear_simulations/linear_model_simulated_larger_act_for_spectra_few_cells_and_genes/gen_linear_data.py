import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import re
import os
import palava
import timeit
import warnings
warnings.filterwarnings('ignore')
import scanpy as sc
import torch
import pandas as pd
import argparse
import random



def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

parser = argparse.ArgumentParser()
parser.add_argument('--seed', type=str, help='random seed')


args = parser.parse_args()

# Access the data directory
seed = int(args.seed)
set_seed(seed)


# Access the data directory
set_seed(seed)


activation_matrix= torch.zeros(4, 5000) 
#q = 0.4# fraction shuffled
#frac_zero_act = 0.85
max_  = 2.5

num_cells_active = 100
#num_zero_act = int(5000 * frac_zero_act)

common_activations =  torch.zeros((5000)) 

# common activations 
activation_matrix[0, :] =  common_activations
activation_matrix[1, :] =  common_activations
activation_matrix[2, :] =  common_activations
activation_matrix[3, :] =  common_activations
#activation_matrix += torch.randn((4, 5000)) 

# random resampling for shuffled values
#rand_ind = random.sample(list(range(5000)), num_cells_active) 
activation_matrix[0, :100] =   torch.rand((num_cells_active)) * max_ 
rand_ind = random.sample(list(range(5000)), num_cells_active)
activation_matrix[1, 1000:1100] =   torch.rand((num_cells_active)) * max_ 
rand_ind = random.sample(list(range(5000)), num_cells_active)
activation_matrix[2, 2000:2100] =   torch.rand((num_cells_active)) *  max_ 
rand_ind = random.sample(list(range(5000)), num_cells_active)
activation_matrix[3, 3000:3100] =   torch.rand((num_cells_active)) * max_ 

# for i in range(4):
#     rand_ind = random.sample(list(range(5000)), num_zero_act)
#     activation_matrix[i, rand_ind] =   torch.rand((num_zero_act)) * 0.001


activation_matrix = activation_matrix.abs()


factor_loadings = torch.zeros(4, 2000)
factor_loadings[0,0:50] =  torch.rand((50)) * max_ + 0.5
factor_loadings[1,500:550] =  torch.rand((50)) * max_ + 0.5
factor_loadings[2,1000:1050] =  torch.rand((50)) * max_ + 0.5 
factor_loadings[3,1500:1550] =  torch.rand((50)) * max_ + 0.5

factor_loadings = factor_loadings.abs() 

gene_set_info = torch.zeros(4, 2000)
gene_set_info[0,0:50] =  torch.ones((50)) * 2
gene_set_info[1,500:550] =  torch.ones((50)) * 2
gene_set_info[2,1000:1050] =  torch.ones((50)) * 2
gene_set_info[3,1500:1550] =  torch.ones((50)) * 2

data = factor_loadings.T @ activation_matrix 
data += (torch.randn(2000, 5000)*(0.10)- 0.05)*max_
data = data.abs() 
data_ = data # saved for log count data
data = (data.T/data.T.sum(axis = 1, keepdim = True))  * torch.ones(5000).log_normal_(mean=11.0, std=0.2).unsqueeze(dim=1)

X_pois =  torch.poisson(data)

adata = ad.AnnData(X=X_pois.detach().numpy())


for fac in [0,1,2,3]:
    adata.obs['steps_' +str(fac+1)] = activation_matrix[fac,:].detach().numpy() 


for fac in [0,1,2,3]: # if gene_set_is_zero element then 1 or else gene_set_is_zero value
    
    adata.var['DEFac_' +str(fac+1)] = (gene_set_info.T[:, fac] + (gene_set_info.T[:, fac] == 0)*1).detach().numpy()

#X_pois = (X_pois/X_pois.sum(axis = 1, keepdim = True))  * torch.ones(5000).log_normal_(mean=11.0, std=0.2).unsqueeze(dim=1)

adata.layers['counts'] = X_pois.numpy()

adata.uns['activations'] =  activation_matrix.detach().numpy()
adata.uns['loadings'] =  factor_loadings.T.detach().numpy()
adata.obs['lb_size'] = X_pois.sum(axis =1) 
median_lb_size = torch.median(X_pois.sum(axis = 1).T)

adata.layers['logcounts'] = data_.T.detach().numpy() #torch.log(((X_pois.T / X_pois.sum(axis = 1)).T *median_lb_size  ) + 1).detach().numpy()

    

num_facs = 4
num_cells = 5000
num_genes = 2000


adata.write('data/sim_anndata_linear_model_' + str(num_cells)+'_cells_'+str(num_genes)+'_genes_'+str(num_facs)+'_factors_'+str(seed)+'_seed.h5ad')

