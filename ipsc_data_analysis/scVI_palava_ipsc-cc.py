import scanpy as sc
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import torch
import warnings
import palava
import json
import sys
import h5py

import os
import gc
from palava import settings
import matplotlib.patches as mpatches
import time 
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import argparse

def set_seed(seed):
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
set_seed(0)

num_unann = int(40)


parser = argparse.ArgumentParser()
parser.add_argument('--lambda_', type=str, help='Lambda values')

parser.add_argument('--str_ab_exp', type=str, help='String about experiment, filders will be names based on this')

parser.add_argument('--list_of_nonlin_factors', type=str, help='lst of ind of nonlinear factors')

parser.add_argument('--palava_width', type=str, help='width of palava decoder')

parser.add_argument('--lambda_marker_genes', type=str, help='width of palava decoder')


parser.add_argument('--dir', type=str, help='directory to save figs')

parser.add_argument('--seed', type=str, help='seed')

args = parser.parse_args()

seed = int(args.seed)


# Access the data directory
hyper = args.lambda_

#string about experiment
str_ab_exp = args.str_ab_exp

palava_width = int(args.palava_width)


lambda_marker_genes = float(args.lambda_marker_genes)

print('palava_width', palava_width)

list_of_nonlin_factors = args.list_of_nonlin_factors

print('nonlin_facs', list_of_nonlin_factors)


print('string about experiment: ', str_ab_exp)

fdir = 'Out_files_and_results_cc/'
fdir += str_ab_exp

fdir += '/' +args.dir
fdir += '/'

#try:
#  os.makedirs(fdir)
#except:
#  print('dir already made')

# Write all the figures to single pdf

adata = sc.read('data/iPSC_with_pathways_cell_metadata_and_stages.h5ad')
pathways = [torch.tensor(i) for i in adata.uns['pathways_9000hvg']]

gene_names = adata.var

num_genes = len(adata.var)
pathway_names = adata.uns['pathway_names']
    
settings.seed = seed


SCVI_palava = palava.model.SCVI_palava

SCVI_palava.setup_anndata(adata, layer = 'counts')

pathways_bool = pathways

num_ann = len(pathway_names)


if list_of_nonlin_factors == 'all':
    list_of_nonlin_factors = '-'.join([str(i) for i in range(num_ann+num_unann)])
    decoder_bias = True


if list_of_nonlin_factors == 'none':
    nonlin_facs = []
    decoder_bias = False

elif list_of_nonlin_factors == 'none_no_bias':
    nonlin_facs = []
    decoder_bias = False
else:
    nonlin_facs = [int(i) for i in list_of_nonlin_factors.split('-')]
    decoder_bias = True

    

nonlinear  = [0 for _ in range(num_ann + num_unann)]

print('Non linear factors')
for i in nonlin_facs:
    nonlinear[i] = 1
    
plan_kwargs = {'lr' :  1e-4,'n_epochs_kl_warmup': 400,'n_steps_kl_warmup': None, 'max_kl_weight' : 1.0}

print(plan_kwargs)


lam = float(hyper)

lam_ = []
lam_  = [lam] * 51 # 50 hallmark
lam_  = [0.15] * 3 + lam_ # marker
lam_ = lam_ + [lam/2] * num_unann # unann

lam_[16] = lambda_marker_genes
lam_[21] = lambda_marker_genes

pathway_names_plot = [pathway_names[i].replace('_', ' ' ).capitalize() + ' ['+str(i)+']'  for i in range(len(pathway_names))] + ['Unannotated factor '+ str(i + 1) + ' ['+str(i+num_ann)+']'  for i in range(num_unann)]


marker_gene_set_lst = [0 ,1, 2]

#gene-cell
scvi_palava = SCVI_palava(adata, nonlinear  = nonlinear, marker_gene_set_lst = marker_gene_set_lst, n_annotated_latent = len(pathways_bool), n_unannotated_latent = num_unann, lam_lst = lam_ , other_n_latent = 40, pathways_bool = pathways_bool, dispersion = 'gene', palava_n_hidden = palava_width, palava_n_layers = 1,  n_layers = 3, use_batch_norm = 'none',use_layer_norm = 'none',  momentum_train = 0.1,  weight_decay_train = 0.00001, decoder_bias = decoder_bias)

print(scvi_palava)


scvi_palava.train(850, gradient_clip_val = 1e4, train_size = 0.9, batch_size = 100, plan_kwargs = plan_kwargs )


scvi_palava.save(fdir + '/latent_and_slope_data/scvi_model' )


# Ws= scvi_palava.module.decoder.px_decoder.annotated_generative_model.out_Ws

# for i in range(len(pathways)):
    
#     print(torch.topk(Ws[i],10).values)

   
adata.obsm["X_scVI"] = scvi_palava.get_latent_representation()




latent = adata.obsm["X_scVI"][:,:len(pathways_bool)]

cell_stages = adata.obs['cell_differentiation'].tolist()

sorted_ind = np.argsort(adata.obsm['pseudo'])

learned_activations = np.transpose(latent[sorted_ind, :])

fig, ax = plt.subplots(figsize=(15, 7.5))

plt.imshow(learned_activations,   aspect = 'auto', cmap='bwr', interpolation = 'none', rasterized=True)
plt.title("Learned activation")
plt.xlabel('Cells')
plt.ylabel('Factors')
_ = plt.yticks([i for i in range(len(pathway_names_plot))], pathway_names_plot, fontsize = 7.5)
cbar = plt.colorbar()
plt.savefig(fdir+'image plot of latent space.png', dpi =150)
plt.close()
plt.clf()
gc.collect()
#pdf_pages.savefig(bbox_inches='tight')



pathways_names_to_plot = ['HALLMARK_G2M_CHECKPOINT', 'HALLMARK_E2F_TARGETS', 'HALLMARK_HYPOXIA'] + ['pluir_markers','mesendo_markers','defendo_markers','respiratory_gene_set']
ind_of_genes = [i for i in range(len(pathway_names)) if pathway_names[i] in  pathways_names_to_plot ]

n_lat = len(ind_of_genes) #adata.obsm["X_scVI"].shape[1]
latent = adata.obsm["X_scVI"]
# Assuming 'steps' and 'latent' are the two arrays with dimensions 5000x4

# Create subplots for each combination of columns
fig, axs = plt.subplots(n_lat, n_lat, figsize=(12, 12))

# Iterate over each combination of columns
for i in range(n_lat):
    for j in range(n_lat):
        
        i_ind = ind_of_genes[i]
        j_ind = ind_of_genes[j]   
        # Scatter plot between steps[:, i] and latent[:, j]
        x = latent[:, i_ind]
        y = latent[:, j_ind]
        category_colors = {'ips': 'red', 'mesendo': 'blue', 'defendo': 'green', 'no_stage': 'yellow'}

        # Map categories to colors
        colors = [category_colors[cell_stage] for cell_stage in cell_stages]

        cmap = plt.cm.viridis

        # Create the scatter plot
       
        axs[i, j].scatter(y, x, s = 0.5, c= colors, edgecolor = 'none')
        if i == n_lat - 1:
            axs[i, j].set_xlabel(pathway_names_plot[j_ind].replace(' ', '\n'), fontsize=10)
        if j == 0:
            axs[i, j].set_ylabel(pathway_names_plot[i_ind].replace(' ', '\n'), fontsize=10)
# Adjust the spacing between subplots
plt.tight_layout()
legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, label=cell_stages)
                   for cell_stages, color in category_colors.items()]
plt.legend(handles=legend_elements)
# Show the plot
plt.savefig(fdir+'matrix plot of some ann latent.png', dpi =150)
plt.close()
plt.clf()
gc.collect()



# Helper functions for the plots
def sort_based_on_first_row(matrix):
    sorted_indices = np.argsort(-matrix[0])
    # Rearrange the columns based on the sorted indices
    sorted_matrix = matrix[:, sorted_indices]
    return  sorted_matrix



# def pathways_with(markers):    
    
#     result = analysis.identifiers(ids=markers,page_size='-1', page='-1',p_value =1)
#     token = result['summary']['token']
#     str_ = markers + ' \n '
#     for path_n in [path['name'] for path in result['pathways']]:
#         str_ +='    ' + path_n + ' \n '
#     return str_ 



pathway_names_plot_pres = [pathway_names[i].replace('_', ' ' ).capitalize()  for i in range(len(pathway_names))] + ['Unannotated factor '+ str(i + 1)   for i in range(num_unann)]



fig, ax = plt.subplots(figsize=(15, 7.5))

plt.imshow(learned_activations,   aspect = 'auto', cmap='bwr', interpolation = 'none', rasterized=True)
plt.title("Learned activation")
plt.xlabel('Cells')
plt.ylabel('Factors')
_ = plt.yticks([i for i in range(len(pathway_names_plot_pres))], pathway_names_plot_pres, fontsize = 7.5)
cbar = plt.colorbar()
plt.savefig(fdir+'image plot of latent space pres.png', dpi =150)
plt.close()
plt.clf()
gc.collect()





pathways_names_to_plot = ['pluir_markers','mesendo_markers','defendo_markers']
ind_of_genes = [i for i in range(len(pathway_names)) if pathway_names[i] in  pathways_names_to_plot ]

n_lat = len(ind_of_genes) #adata.obsm["X_scVI"].shape[1]
latent = adata.obsm["X_scVI"]

# Assuming 'steps' and 'latent' are the two arrays with dimensions 5000x4

# Create subplots for each combination of columns
fig, axs = plt.subplots(n_lat, n_lat, figsize=(12, 12))

# Iterate over each combination of columns
for i in range(n_lat):
    for j in range(n_lat):
        
        i_ind = ind_of_genes[i]
        j_ind = ind_of_genes[j]   
        # Scatter plot between steps[:, i] and latent[:, j]
        x = latent[:, i_ind]
        y = latent[:, j_ind]
        category_colors = {'ips': 'red', 'mesendo': 'blue', 'defendo': 'green', 'no_stage': 'yellow'}

        # Map categories to colors
        colors = [category_colors[cell_stage] for cell_stage in cell_stages]

        cmap = plt.cm.viridis

        # Create the scatter plot
       
        axs[i, j].scatter(y, x, s = 1, c= colors, edgecolor = 'none')
        if i == n_lat - 1:
            axs[i, j].set_xlabel(pathway_names_plot[j_ind].replace(' ', '\n'), fontsize=10)
        if j == 0:
            axs[i, j].set_ylabel(pathway_names_plot[i_ind].replace(' ', '\n'), fontsize=10)
# Adjust the spacing between subplots
plt.tight_layout()
legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, label=cell_stages)
                   for cell_stages, color in category_colors.items()]
plt.legend(handles=legend_elements)
# Show the plot
plt.savefig(fdir+'matrix plot of marker latent.png', dpi =150)
plt.close()
plt.clf()
gc.collect()



##saving ann dec using torch.save is much easier!!


try:
  os.mkdir(fdir + '/latent_and_slope_data')
except:
  print('dir already made')

np.save(fdir +'/latent_and_slope_data/latent.npy', latent)


factor_importance = scvi_palava.factor_importance(adata, num_batchs =300)


#reverse order and top 10
index_asc_order = np.argsort(factor_importance)[::-1][:10]

plt.figure(figsize=(15, 5))
plt.bar(np.array(pathway_names_plot)[index_asc_order], factor_importance[index_asc_order] )
plt.xticks(rotation=90)
plt.savefig(fdir +'sorted factor importance num_batch 200.png', dpi = 150,bbox_inches='tight')

factor_importance_scores_dict = {'pathway_names_plot':pathway_names_plot,'factor_importance':factor_importance }

np.save(fdir+'factor_importance_scores_dict.npy', factor_importance_scores_dict, allow_pickle=True)



if True:    # Use top 20 factors
    for fac in np.argsort(factor_importance)[::-1][:20]:

        factors = [fac]

        slopes_for_all_fac_lst = scvi_palava.spline_approach(factors = factors)
        slopes_for_all_fac  = {i:j for i, j  in zip(factors,slopes_for_all_fac_lst)}
        


        data_1 = np.quantile(abs(slopes_for_all_fac[fac]), q=0.95, axis=0)
        data_1_and_de_no_error = np.concatenate((data_1[None], pathways_bool[fac].numpy()[None],gene_names.index.to_numpy()[None]), axis=0)
        data_1_and_de_no_error = sort_based_on_first_row(data_1_and_de_no_error)
    
        n_top_genes = 75
    
        y = data_1_and_de_no_error[0][:n_top_genes]
        x = np.array([i for i in range(n_top_genes)])
        
        plt.figure(figsize=(10, 5))
        plt.scatter(x, y, s = 5, edgecolor='none')
    
        # plotting pathway genes
        pathway_gene_ind = data_1_and_de_no_error[1][:n_top_genes] == 1
        #x_pathway = (data_1_and_de_no_error[1][:n_top_genes] * x)[pathway_gene_ind]
        #y_pathway = (np.ones(n_top_genes) * -0.05)[pathway_gene_ind]
        # plt.scatter(x_pathway, y_pathway, c='black', marker='|', label='In gene set',rasterized=True)
    
        plt.xlabel('Ranking')
        plt.ylabel('Summarised slope values')
        plt.title(pathway_names_plot[fac] + ' Factor ' + str(fac))
        
        #pdf_pages.savefig(bbox_inches='tight')
        
        
        top_gene_names = [element.split('_')[1] if '_' in element else '' for element in data_1_and_de_no_error[2][:n_top_genes]]
        for i, txt in enumerate(top_gene_names):
             color = ['r', 'g'][int(data_1_and_de_no_error[1][i])]      
             plt.annotate(txt, (x[i], y[i]), ha = 'right', va = 'bottom', rotation = 90, fontsize = 7.5, color = color)
        geneset = mpatches.Patch(color='g', label='Gene set genes')
        nongeneset = mpatches.Patch(color='r', label='Non gene set genes')
        
        plt.legend(handles=[geneset, nongeneset])
        #plt.show()
        plt.savefig(fdir+pathway_names_plot[fac]+ ' rank plot.png', dpi =150)
    
        plt.close()
    
        
    
        plt.clf()
        gc.collect()
         ####################################





with open(fdir+'done.txt', 'w') as file:
    # Write 'done!' to the file
    file.write('done!')
