The original .h5ad file exceeded GitHub and Git LFS limits (>2 GB). To fix this, a duplicate layer was removed. You can recreate the required counts layer with the following code:

<pre lang="python"><code>
import scanpy as sc 
# Load reduced file (without duplicate counts layer)
adata = sc.read('data/no_counts_layer_iPSC_with_pathways_cell_metadata_and_stages.h5ad')
# Restore the 'counts' layer (required by pipeline) 
adata.layers['counts'] = adata.X 
# Save the updated file 
adata.write('data/iPSC_with_pathways_cell_metadata_and_stages.h5ad') </code></pre>