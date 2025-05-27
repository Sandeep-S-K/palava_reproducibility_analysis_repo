import scanpy as sc

import matplotlib.pyplot as plt
import numpy as np
import torch
import palava


import gc
from palava import settings
import argparse

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

fdir = 'Out_files_and_results/'
fdir += str_ab_exp

fdir += '/' +args.dir
fdir += '/'

adata = sc.read('data/mouse_pancreas_updated_beta_with_pdx1.h5ad')
pathways = [torch.tensor(i) for i in adata.uns['pathways_5000hvg']]

gene_names = adata.var

num_genes = len(adata.var)
pathway_names = adata.uns['pathway_names']
    
settings.seed = seed


SCVI_palava = palava.model.SCVI_palava

SCVI_palava.setup_anndata(adata, layer = 'counts')

pathways_bool = pathways



num_ann = len(pathway_names)

decoder_bias = True

if list_of_nonlin_factors == 'all':
    list_of_nonlin_factors = '-'.join([str(i) for i in range(num_ann+num_unann)])


if list_of_nonlin_factors == 'none':
    decoder_bias =False

    nonlin_facs = []
elif list_of_nonlin_factors == 'none_no_bias':
    nonlin_facs = []
    decoder_bias =False
else:
    nonlin_facs = [int(i) for i in list_of_nonlin_factors.split('-')]




nonlinear  = [0 for _ in range(len(pathways_bool) + num_unann)]
print('non linear len', len(nonlinear))
print('Non linear factors')
for i in nonlin_facs:
    nonlinear[i] = 1

for i in range(num_ann):
    print(i, pathway_names[i], nonlinear[i])
    
plan_kwargs = {'lr' :  1e-4,'n_epochs_kl_warmup': 400,'n_steps_kl_warmup': None, 'max_kl_weight' : 1.0}

print(plan_kwargs)


lam = float(hyper)

#lam_ = []
lam_  = [lam] * len(pathway_names) # 50 hallmark + marker gene sets genes
lam_ = lam_ + [lam/2] * num_unann # unann

pathway_names_plot = [pathway_names[i].replace('_', ' ' ).capitalize() + ' ['+str(i)+']'  for i in range(len(pathway_names))] + ['Unannotated factor '+ str(i + 1) + ' ['+str(i+num_ann)+']'  for i in range(num_unann)]

marker_gene_set_lst = [50, 51, 52, 53]

for i in marker_gene_set_lst:
    lam_[i] = lambda_marker_genes
# 50, 51, 52 are anchor genes

scvi_palava = SCVI_palava(adata, nonlinear  = nonlinear, marker_gene_set_lst = marker_gene_set_lst, n_annotated_latent = len(pathways_bool), n_unannotated_latent = num_unann, lam_lst = lam_, other_n_latent = 40, pathways_bool = pathways_bool, dispersion = 'gene', palava_n_hidden = palava_width, palava_n_layers = 1,  n_layers = 3, use_batch_norm = 'none', use_layer_norm = 'none',  momentum_train = 0.1,  weight_decay_train = 0.0001, decoder_bias = decoder_bias, non_neg_decoder = True, just_out_weight_decay = True) #, gene_likelihood = 'nb')
print(scvi_palava)

scvi_palava.train(6000, gradient_clip_val = 1e4, train_size = 1.0, batch_size = 100, plan_kwargs = plan_kwargs )


scvi_palava.save(fdir + '/latent_and_slope_data/scvi_model' )


Ws = scvi_palava.module.decoder.px_decoder.annotated_generative_model.out_Ws


factor_importance = scvi_palava.factor_importance(adata, num_batchs =300)


#reverse order and top 10
index_asc_order = np.argsort(factor_importance)[::-1][:10]

plt.figure(figsize=(15, 5))
plt.bar(np.array(pathway_names_plot)[index_asc_order], factor_importance[index_asc_order] )
plt.xticks(rotation=90)
plt.savefig(fdir +'sorted factor importance num_batch 300.png', dpi = 150,bbox_inches='tight')

factor_importance_scores_dict = {'pathway_names_plot':pathway_names_plot,'factor_importance':factor_importance }

np.save(fdir+'factor_importance_scores_dict.npy', factor_importance_scores_dict, allow_pickle=True)






adata.obsm["X_scVI"] = scvi_palava.get_latent_representation()




latent = adata.obsm["X_scVI"][:,:len(pathways_bool)]



learned_activations = np.transpose(latent)

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




with open(fdir+'done.txt', 'w') as file:
    # Write 'done!' to the file
    file.write('done!')
