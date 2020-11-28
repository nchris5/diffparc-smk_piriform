# **Snakemake workflow: diffparc-smk_piriform**

Snakemake workflow overview for diffusion-based parcellation of the piriform cortex:
1. ```rule hcp7T_Diffusion_bedpostx_gpu```: Generates white matter fibre orientations within each brain voxel for each subject's HCP 7T diffusion data via FSL's ```bedpostx_gpu```. 
2. ```rule hcp_mmk```: Generates target segmentation 3D nifti volumes (bilateral) from HCP32k surf gifti files (180 regions)
3. ```rule diffparc_smk```: Resamples the piriform seed/hcp-mmk180_targets-->HCP7TDiffusionResolution and subsequently performs probabilistic tractography from the piriform seed in each subject's space to hcp-mmk180_targets via ```probtrackx2_gpu```

    3b. Brings connectivity data for each subject's seed voxels into template and performs spectral clustering on the concatenated feature vectors to parcellate into k regions
    
    3c. Creates node and edges tables to use as the input to the Gephi software for each of the clustering solutions
    
4. ```rule tractmap``` Performs probabilistic tractography on each individual cluster to generate tractography maps for each cluster following spectral clustering (useful for visualization of each cluster's individual connectivity)

## Inputs in config/config.yml:
1. ```participants.tsv```: Target subject IDs that have both HCP 7T diffusion and be used in creation of a subject specific template from: [ants_build_template_smk](https://github.com/akhanf/ants_build_template_smk) , (see input 5)
2. ```template_prob_seg```: Probabilistic segmentation as a 3D nifti for individual piriform hemisphere to process
3. ```prob_seg_threshold```: Value to threshold the probabilistic segmentation at (default is 0.5, e.g. 0.2 would include 80% of the probseg)
4. ```template```: Name of template to use for filenaming within the pipeline
5. ANTS transformations from template T1w space to/from each subject's T1w, i.e. from: [ants_build_template_smk](https://github.com/akhanf/ants_build_template_smk)
 
   * ```ants_affine_mat```: Forward affine from subject-->template
   * ```ants_invwarp_nii```: Warp from template-->subject
   * ```ants_warp_nii```: Warp from subject-->template
   * ```ants_ref_nii```: Study-specific T1w created from: [ants_build_template_smk]

6. ```targets_txt```: List of ROIs in atlas segmentations to use as targets for piriform connectivity (default hcp_mmp_sym_180)
7. ```max_k```: Maximum number of clusters to create during spectral clustering (will create all cluster solutions <= max_k)
8. ```probtrack```: Tractography parameters, seed_resolution defines the resolution to resample all other inputs to
9. ```probtrack_tractmap```: Tractography parameters to generate individual tractmaps for each cluster following spectral clustering (useful for visualization)
10. ```in_freesurfer```: Path to subject freesurfer dir from HCP_1200 release
11. ```dwi_preproc_dir```: Path to base directory containing the HCP 7T diffusion data for all subjects

   * ```in_dwi_nii```
   * ```in_dwi_bval```
   * ```in_dwi_bvec```
   * ```in_dwi_grad_dev```
   * ```in_dwi_mask_nii```

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


Activate the conda environment: ```conda activate snakemake``` or add this line into your bash startup

For installation details, see the instructions in the Snakemake documentation.

#### Step 4: Execute workflow
* Option A)
Execute the workflow locally via: ```snakemake --use-conda --use-singularity --cores $N```
using $N cores, or, run it in a cluster environment via ```snakemake --use-conda --use-singularity --cluster qsub --jobs 100``` or ```snakemake --use-conda --use-singularity --drmaa --jobs 100```

* Option B)
If you are using Compute Canada, you can use the cc-slurm profile, which submits jobs and takes care of requesting the correct resources per job (including GPUs). Once it is set-up with cookiecutter, run: ```snakemake --profile cc-slurm```

* Option C)
Or, with neuroglia-helpers can get a 1-GPU, 8-core, 32gb node and run locally there. First, get a node with a GPU (default 8-core, 32gb, 3 hour limit): ```regularInteractive -g```. Then, run: ```snakemake --use-conda --use-singularity --cores 8 --resources gpu=1 mem=32000```

##### Step 4b: Executing workflow on Compute Canada via the cc-slurm profile in iterations
If you have a large number of subjects, we recommend you run the workflow in the following stages:
1. bedpostx_gpu: via ```snakemake --profile cc-slurm --until perform_bedpostx_gpu```
2. generate target segs and probtrackx2_gpu: via ```snakemake --profile cc-slurm --until run_probtrack```
3. piriform diffusion parcellation with spectral_clustering in template space on concatenated subjects and creation of Gephi input nodes/edges tables: via ```snakemake --profile cc-slurm --until create_gephi_input_nodes_and_edge_tables```
4. probtrackx2_gpu run again on each cluster to generate visualizations of each cluster's connectivity: via ```snakemake --profile cc-slurm --until vote_tractmap_template```

Be sure to run ```snakemake --profile cc-slurm --until <iteration name previously run> -np``` after running each iteration to ensure that all outputs were sucessfully created, before moving on to the next iteration.
