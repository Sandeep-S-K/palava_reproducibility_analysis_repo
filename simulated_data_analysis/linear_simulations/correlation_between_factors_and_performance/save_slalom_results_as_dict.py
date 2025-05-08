import pandas as pd
import argparse
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('--data_directory', type=str, help='Path to data directory')


args = parser.parse_args()

# Access the data directory
data_directory = args.data_directory
print(f'Data directory: {data_directory}')



data_dr_results = data_directory.replace('data/', 'slalom_intermediate_results_R/slalom_results_')
data_dr_results = data_dr_results.replace('.h5ad', '')


latent  = pd.read_csv(data_dr_results +'_latent.csv').values
training_time  = float(pd.read_csv(data_dr_results +'_training_time.csv').values[0][0])
factor_loadings  = pd.read_csv(data_dr_results +'_factor_loadings.csv').values.T
factor_importance  = pd.read_csv(data_dr_results +'_factor_importance.csv').values.squeeze()
posterior  = pd.read_csv(data_dr_results +'_posterior.csv').values.T
conv  = pd.read_csv(data_dr_results +'_converged.csv').values
results = {}

results['factor_loadings'] = factor_loadings

results['latent_representation'] =  latent
results['factor_importance'] = factor_importance

results['training_time'] = training_time
results['posterior'] = posterior
results['conv'] = conv


name_of_anndata = data_directory.split('/')[-1]

np.save('intermediate_results/slalom_results_'+name_of_anndata.replace('.h5ad','.npy'), results, allow_pickle=True)
