
rule all_bedpostx_gpu:
        input:
                mean_s0samples = expand('results/Diffusion_7T/{subject}.bedpostX/mean_S0samples.nii.gz',subject=subjects),
                merged_f1samples = expand('results/Diffusion_7T/{subject}.bedpostX/merged_f1samples.nii.gz',subject=subjects),
                merged_f2samples = expand('results/Diffusion_7T/{subject}.bedpostX/merged_f2samples.nii.gz',subject=subjects),
                merged_f3samples = expand('results/Diffusion_7T/{subject}.bedpostX/merged_f3samples.nii.gz',subject=subjects),
                merged_ph1samples = expand('results/Diffusion_7T/{subject}.bedpostX/merged_ph1samples.nii.gz',subject=subjects),
                merged_ph2samples = expand('results/Diffusion_7T/{subject}.bedpostX/merged_ph2samples.nii.gz',subject=subjects),
                merged_ph3samples = expand('results/Diffusion_7T/{subject}.bedpostX/merged_ph3samples.nii.gz',subject=subjects),
                merged_th1samples = expand('results/Diffusion_7T/{subject}.bedpostX/merged_th1samples.nii.gz',subject=subjects),
                merged_th2samples = expand('results/Diffusion_7T/{subject}.bedpostX/merged_th2samples.nii.gz',subject=subjects),
                merged_th3samples = expand('results/Diffusion_7T/{subject}.bedpostX/merged_th3samples.nii.gz',subject=subjects)

#Copy 7T Diffusion files to the snakemake results directory so bedpost will subsequently save it there
rule copy_hcp7T_Diffusion:
	input:
		dwi7T_nii = join(config['dwi_preproc_dir'],config['in_dwi_nii']),
                dwi7T_bval = join(config['dwi_preproc_dir'],config['in_dwi_bval']),
                dwi7T_bvec = join(config['dwi_preproc_dir'],config['in_dwi_bvec']),
                dwi7T_grad_dev = join(config['dwi_preproc_dir'],config['in_dwi_grad_dev']),
                dwi7T_mask_nii = join(config['dwi_preproc_dir'],config['in_dwi_mask_nii'])
	output:
		outdir_copy = directory('results/Diffusion_7T/{subject}'),
		dwi7T_nii_copy = 'results/Diffusion_7T/{subject}/data.nii.gz',
		dwi7T_bval_copy = 'results/Diffusion_7T/{subject}/bvals',
		dwi7T_bvec_copy = 'results/Diffusion_7T/{subject}/bvecs',
		dwi7T_grad_dev_copy = 'results/Diffusion_7T/{subject}/grad_dev.nii.gz',
		dwi7T_mask_nii_copy = 'results/Diffusion_7T/{subject}/nodif_brain_mask.nii.gz'
	group: 'bedpost'
	shell:
		'cp {input.dwi7T_nii} {output.dwi7T_nii_copy} && '
		'cp {input.dwi7T_bval} {output.dwi7T_bval_copy} && '
		'cp {input.dwi7T_bvec} {output.dwi7T_bvec_copy} && '
		'cp {input.dwi7T_grad_dev} {output.dwi7T_grad_dev_copy} && '
		'cp {input.dwi7T_mask_nii} {output.dwi7T_mask_nii_copy}'

rule perform_bedpostx_gpu:
	input:	
		indir = 'results/Diffusion_7T/{subject}'
	params:
		container = config['singularity_cuda']
	output:
		outdir = directory('results/Diffusion_7T/{subject}.bedpostX'),
		mean_s0samples = 'results/Diffusion_7T/{subject}.bedpostX/mean_S0samples.nii.gz'
	threads: 8
	resources:
		mem_mb = 32000,
		time = 360,    #6 hours
		gpus = 1       #1 gpu
	log: 'logs/perform_bedpostx_gpu/{subject}.log'
	group: 'bedpost'
	shell:
		'mkdir -p {output.outdir} && singularity exec -e --nv {params.container} bedpostx_gpu {input.indir} -g &> {log}'
