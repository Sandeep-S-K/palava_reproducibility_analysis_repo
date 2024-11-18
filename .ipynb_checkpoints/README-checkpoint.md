# PALAVA simulated data analysis

## Splatter folder

 -  This is the splatter package but it is modified to produced more nonlinear trajectory simulations, the 'bridge' function in the 'splat-simulate.R' script is modified.

## linear_and_nonlinear_palava_on_modified_splatter_data folder
 - This contains the snakemake pipeline that ruins the linear and nonlinear PALAVA analysis on the modified simulated data.
 -  You will need to set up two conda environments, one with r-base and splatter installed the other wiht PALAVA installed. The snakemake file will need to be updated to account for the the names of the two conda environments you create.
 -  You will also need to create a [snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html) to run the pipeline to run with slurm. The snakemake profile I used is below. This code is in my home directory '.config/snakemake/slurm/config.yaml'. This description is very brief.

```
jobs: 50
restart-times: 0
rerun-incomplete: false
max-jobs-per-second: 2
printshellcmds: true

cluster: "sbatch -t {resources.runtime} --account={resources.account} --partition={resources.partition} --mem={resources.mem_mb} -c {resources.cpus} -o {resources.output_dr} --job-name {resources.job_name} --gres {resources.gres}"

default-resources: [cpus=1, mem_mb=2000, runtime=60]
```
   
 -  After that you can run the data the snakemake pipeline will save the results of the analysis.
 -  The notebook in the folder will load the data and create the plots.