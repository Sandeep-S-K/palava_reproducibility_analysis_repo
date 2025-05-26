The original .h5ad file exceeded GitHub and Git LFS limits (>2 GB). To fix this, a duplicate layer was removed. You can recreate the required counts layer with the following code:

<pre lang="python"><code>
import scanpy as sc 
# Load reduced file (without duplicate counts layer)
adata = sc.read('data/no_counts_layer_iPSC_with_pathways_cell_metadata_and_stages.h5ad')
# Restore the 'counts' layer (required by pipeline) 
adata.layers['counts'] = adata.X 
# Save the updated file 
adata.write('data/iPSC_with_pathways_cell_metadata_and_stages.h5ad') </code></pre>

The counts layer and data of the Spectra's .h5ad result object was deleted before saving to reduce space, but the model outputs are preserved as is in `intermediate_results/no_X_no_counts_layer_spectra_adataiPSC_with_pathways_cell_metadata_and_stages.h5ad`.