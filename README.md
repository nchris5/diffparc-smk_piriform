# **Snakemake workflow: diffparc-smk_piriform**

Snakemake workflow for diffusion-based parcellation of the piriform cortex on graham.computecanada

Performs: 
1. bedpostx_gpu on HCP 7T diffusion data
2. Generates target segmentation 3D nifti volumes from HCP32k surf giftis
3. probtrackx2_gpu from the piriform seed in each subject's space to targets
4. Brings connectivity data for each subject's seed voxels into template and performs spectral clustering on the concatenated feature vectors to parcellate into k regions
5. Creates node and edges tables to use as the input to the Gephi software for each of the clustering solutions

## Inputs in config/config.yml:
1. participants.tsv: Target subject IDs that have both HCP 7T diffusion and be used in creation of a subject specific template from: [ants_build_template_smk](https://github.com/akhanf/ants_build_template_smk) 
2. template_prob_seg: Probabilistic segmentation as a 3D nifti for individual piriform hemisphere to process
3. prob_seg_threshold: Value to threshold the probabilistic segmentation at (default is 0.5, e.g. 0.2 would include 80% of the probseg)
4. template: Name of template to use for filenaming within the pipeline
5. ANTS transformations from template T1w space to/from each subject's T1w, i.e. from: [ants_build_template_smk](https://github.com/akhanf/ants_build_template_smk)
 
   * ants_affine_mat
   * ants_invwarp_nii
   * ants_warp_nii
   * ants_ref_nii

6. targets_text: List of ROIs in atlas segmentations to use as targets for piriform connectivity (default hcp_mmp_sym_180)
7. max_k: Maximum number of clusters to create during spectral clustering (will create all cluster solutions <= max_k)
8. probtrack: Tractography parameters, seed_resolution defines the resolution to resample all other inputs to
9. probtrack_tractmap: Tractography parameters to generate individual tractmaps for each cluster following spectral clustering (useful for visualization)
10. in_freesurfer: Path to subject freesurfer dir from HCP_1200 release
11. dwi_preproc_dir: Path to base directory containing the HCP 7T diffusion data for all subjects

   * in_dwi_nii
   * in_dwi_bval
   * in_dwi_bvec
   * in_dwi_grad_dev
   * in_dwi_mask_nii

### Required Singularity Containers:
* singularity_neuroglia: FSL, ANTS, NiftyReg, gnu-parallel
* Freesurfer
* Connectome workbench
* FSL6 with CUDA container (for running gpu versions of bedpostx and probtrackx2)

### Authors:
* Ali Khan @ akhanf
* Nick Christidis @ nchris5

## Usage:
If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

#### Step 1: Obtain a copy of this workflow
1. Create a new github repository using this workflow as a template.
2. Clone the newly created repository to your local system, into the place where you want to perform the data analysis.

#### Step 2: Configure workflow
Configure the workflow according to your needs via editing the files in the config/ folder. Adjust config.yml to configure the workflow execution, and participants.tsv to specify your subjects.

#### Step 3: Install and Activate Snakemake
Install Snakemake using conda: ```conda create -c bioconda -c conda-forge -n snakemake snakemake```
For installation details, see the instructions in the Snakemake documentation.
Activate the conda environment: ```conda activate snakemake``` or add this line into your bash startup

#### Step 4: Execute workflow
* OptionA)
Execute the workflow locally via: ```snakemake --use-conda --use-singularity --cores $N```
using $N cores or run it in a cluster environment via
```
snakemake --use-conda --use-singularity --cluster qsub --jobs 100
```
or
```
snakemake --use-conda --use-singularity --drmaa --jobs 100
```
If you are using Compute Canada, you can use the cc-slurm profile, which submits jobs and takes care of requesting the correct resources per job (including GPUs). Once it is set-up with cookiecutter, run:
```
snakemake --profile cc-slurm
```
Or, with neuroglia-helpers can get a 1-GPU, 8-core, 32gb node and run locally there. First, get a node with a GPU (default 8-core, 32gb, 3 hour limit):
```
regularInteractive -g
```
Then, run:
```
snakemake --use-conda --use-singularity --cores 8 --resources gpu=1 mem=32000
```
See the Snakemake documentation for further details.

#### Step 5: Investigate results
After successful execution, you can create a self-contained interactive HTML report with all results via:

```
snakemake --report report.html
```
This report can, e.g., be forwarded to your collaborators. An example (using some trivial test data) can be seen here.

#### Step 6: Commit changes
Whenever you change something, don't forget to commit the changes back to your github copy of the repository:
```
git commit -a
git push
```

#### Step 7: Obtain updates from upstream 
Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: ```git remote add -f upstream git@github.com:snakemake-workflows/ants_build_template_smk.git or git remote add -f upstream https://github.com/snakemake-workflows/ants_build_template_smk.git``` if you do not have setup ssh keys.
2. Update the upstream version: ```git fetch upstream```.
3. Create a diff with the current version: ```git diff HEAD upstream/master workflow > upstream-changes.diff```.
4. Investigate the changes: ```vim upstream-changes.diff```.
5. Apply the modified diff via: ```git apply upstream-changes.diff```.
6. Carefully check whether you need to update the config files: ```git diff HEAD upstream/master config```. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.

#### Step 8: Contribute back
In case you have also changed or added steps, please consider contributing them back to the original repository:

Fork the original repo to a personal or lab account.
Clone the fork to your local system, to a different place than where you ran your analysis.
Copy the modified files from your analysis to the clone of your fork, e.g., ```cp -r workflow path/to/fork```. Make sure to not accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
Commit and push your changes to your fork.
Create a pull request against the original repository.
