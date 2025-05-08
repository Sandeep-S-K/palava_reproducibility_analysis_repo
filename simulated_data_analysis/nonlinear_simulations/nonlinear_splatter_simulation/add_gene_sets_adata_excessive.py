
import re
import argparse

import scanpy as sc
import numpy as np



parser = argparse.ArgumentParser()
parser.add_argument('--data_directory', type=str, help='Path to data directory')
parser.add_argument('--EX', type=str, help='List of FP rate')
parser.add_argument('--MS', type=str, help='List of FN rate')

args = parser.parse_args()


MS = int(args.MS)/100
EX = int(args.EX)/100




data_directory =  args.data_directory
adata =  sc.read(data_directory)

# Determine num cells, genes, factors from file name
pattern = r'(\d+)'
matches = re.findall(pattern, data_directory)
num_cells, num_genes, num_factors = map(int, matches[:3])
print(f'num_cells: {num_cells}, num_genes: {num_genes}, Factors: {num_factors}')



gt_gene_sets = adata.var.values.copy()

# If defac is one, it is not a path gene 
gt_gene_sets[gt_gene_sets == 1] = 0  

# All other genes are path genes
gt_gene_sets[gt_gene_sets != 0] = 1
#gt_gene_sets = torch.tensor(pathways.copy())


FN_genes = np.full(gt_gene_sets.shape, False)
FP_genes = np.full(gt_gene_sets.shape, False)

gt_gene_sets_with_error = np.empty(gt_gene_sets.shape)

for fac in range(4):
    gt_gene_sets_fac = gt_gene_sets[:, fac]
    
    TP_genes_fac= np.where(gt_gene_sets_fac == 1)[0]
    
    
    FN_genes_fac = np.random.choice(TP_genes_fac, size=int(TP_genes_fac.shape[0] * MS), replace=False)

    FN_genes[FN_genes_fac, fac]  = True
    
    TN_genes_fac = np.where(gt_gene_sets_fac == 0)[0]
    
    FP_genes_fac = np.random.choice(TN_genes_fac, size=int(TP_genes_fac.shape[0] * EX), replace=False)
    
    FP_genes[FP_genes_fac, fac]  = True

    pathway_with_errors = gt_gene_sets_fac.copy()
    
    pathway_with_errors[FN_genes_fac] = 0
    
    pathway_with_errors[FP_genes_fac] = 1

    gt_gene_sets_with_error[:,fac] = pathway_with_errors
    




    

print('Total num of pathway genes (ground truth)', np.sum(gt_gene_sets, axis = 0))
print('Total num of  genes set (with erros)', np.sum(gt_gene_sets_with_error, axis = 0))


FP_genes = np.logical_and((gt_gene_sets_with_error == 1), (gt_gene_sets == 0))
FN_genes = np.logical_and((gt_gene_sets_with_error == 0), (gt_gene_sets == 1))



uns_dict = adata.uns

for i in range(gt_gene_sets_with_error.shape[1]):
    pathway = gt_gene_sets_with_error[:,i]
    print('Pathway '+ str(i))
    print('False Negatives rate: ', np.sum(gt_gene_sets_with_error[gt_gene_sets==1]==0)/np.sum(gt_gene_sets==1))
    print('False Positives rate: ', np.sum(gt_gene_sets_with_error[gt_gene_sets==0]==1)/np.sum(gt_gene_sets==0))
    uns_dict['FN rate pathway '+str(i)] = np.sum(gt_gene_sets_with_error[gt_gene_sets==1]==0)/np.sum(gt_gene_sets==1)
    uns_dict['FP rate pathway '+str(i)] = np.sum(gt_gene_sets_with_error[gt_gene_sets==0]==1)/np.sum(gt_gene_sets==0)
    print()



uns_dict.update({'FP_genes' : FP_genes, 'FN_genes' : FN_genes, 'Pathways with error' : gt_gene_sets_with_error, 'defac': adata.var.values.copy(), 'Ground truth gene sets': gt_gene_sets})




to_save_dr = data_directory.split('.')[0]+ '_gene_sets_formatted_' +'MS_'+args.MS +'_EX_'+args.EX + '.h5ad'

adata.uns = uns_dict
print(to_save_dr)
adata.write(to_save_dr)
print('End of script')
