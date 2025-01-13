# PALAVA Simulated Data Analysis  

## Splatter Folder  

This folder contains the Splatter package, which has been modified to produce more nonlinear trajectory simulations. The `bridge` function in the `splat-simulate.R` script was updated to enable these simulations.  

## Linear and Nonlinear PALAVA on Modified Splatter Data  

This folder contains the Snakemake pipeline for running linear and nonlinear PALAVA analyses on the modified simulated data.  The snakemake verison used is `7.22.0`.

### Setup  

1. **Create Conda Environments**  
   - Set up two Conda environments:
   -  Create a conda invironment `conda create -n r-modified-splatter-env`. After activating the environment (`conda activate r-modified-splatter-env`) install the required packages using teh command `conda install -c conda-forge -c bioconda -c defaults r-base bioconductor-scater bioconductor-zellkonverter r-devtools r-reticulate anndata`. The specific version of the packages used for our analysis can be seen in the `r_mod_splatter_env.yml` file in environments. 
     - Install the local version of the splatter package from the repository using devtools. First, ensure the r-base environment is active. Then, navigate to the repository directory in an active R session and run `devtools::install_local("splatter")`.   
     - Create another conda environment named `palava-env` and install  `PALAVA` along with Scanpy scanpy using `pip install scanpy==1.9.8`.  
   - Update the Snakemake file to reflect the names of these environments that modified splatter and PALAVA are installed.  

2. **Create a Snakemake Profile**  
   - To run the pipeline with SLURM, you need a [Snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html).  
   - Below is an example slurm profile. Creates a configuration file at `.config/snakemake/slurm/config.yaml` with:  

     ```yaml
     jobs: 50
     restart-times: 0
     rerun-incomplete: false
     max-jobs-per-second: 2
     printshellcmds: true

     cluster: "sbatch -t {resources.runtime} --account={resources.account} --partition={resources.partition} --mem={resources.mem_mb} -c {resources.cpus} -o {resources.output_dr} --job-name {resources.job_name} --gres {resources.gres}"

     default-resources: [cpus=1, mem_mb=2000, runtime=60]
     ```  

   - Refer to the Snakemake documentation for further details on setting up profiles.  
   - Update the rules and profile to reflect your system as desired (e.g partition, etc)

### Running the Pipeline  

Once the environments and profile are configured:  

- Use Snakemake to generate the data by running the pipeline. The pipeline can be run by executing  `snakemake --profile slurm --use-conda --rerun-triggers mtime`  
- The pipeline will save the results as dictionaries using NumPy.  

### Analyzing Results  

The provided notebook will load the data and create plots to visualise the results.  


### TODO 
- Remove unncessary cells from the fetal erythroposesis ananlysis notebook

