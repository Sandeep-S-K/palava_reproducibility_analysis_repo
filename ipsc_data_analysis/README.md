Run to format data properly, it was too large for the GIT LFS, thus removed duplicate layer but  required for pipeline.

'''
adata = sc.read('data/no_counts_layer_iPSC_with_pathways_cell_metadata_and_stages.h5ad') 
adata.layers['counts'] = adata.X
adata.write('data/iPSC_with_pathways_cell_metadata_and_stages.h5ad')
'''