#Warp group-based clusters back to each subject's T1wSpace
rule transform_clus_to_subj:
    input:
        cluster_k = expand('results/clustering/group_space-{template}_seed-{seed}_method-spectralcosine_k-{k}_cluslabels.nii.gz',k=range(2,config['max_k']+1),allow_missing=True),
        affine = config['ants_affine_mat'],
        invwarp = config['ants_invwarp_nii'],
        ref = rules.combine_lr_hcp.output.lh_rh
    output:
        cluster_k = expand('results/tractmap/transform_clus_to_subject/sub-{subject}/Native_space_from-{template}_seed-{seed}_method-spectralcosine_k-{k}_cluslabels.nii.gz',k=range(2,config['max_k']+1),allow_missing=True)
    envmodules: 'ants'
    container: config['singularity_neuroglia']
    log: 'logs/transform_clus_to_subject/sub-{subject}_template-{template}_{seed}.log'
    group: 'tractmap'
    threads: 8
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 parallel  --jobs {threads} antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {{1}} -o {{2}}  -r {input.ref} -t [{input.affine},1] -t {input.invwarp} &> {log} :::  {input.cluster_k} :::+ {output.cluster_k}'

rule resample_brainmask_tractmaps:
    input:
        dwi = 'results/Diffusion_7T/{subject}.bedpostX/mean_S0samples.nii.gz'
    params:
        seed_resolution = config['probtrack_tractmap']['seed_resolution']
    output:
        mask = 'results/tractmap/masks/sub-{subject}/brain_mask_dwi.nii.gz',
        mask_res = 'results/tractmap/masks/sub-{subject}/brain_mask_dwi_resampled_tractmap.nii.gz'
    container: config['singularity_neuroglia']
    log: 'logs/resample_brainmask_tractmap/sub-{subject}.log'
    group: 'tractmap'
    shell:
        'fslmaths {input.dwi} -bin {output.mask} &&'
        'mri_convert {output.mask} -vs {params.seed_resolution} {params.seed_resolution} {params.seed_resolution} {output.mask_res} -rt nearest &> {log}'

rule resample_clus_seed:
    input:
        cluster_k = 'results/tractmap/transform_clus_to_subject/sub-{subject}/Native_space_from-{template}_seed-{seed}_method-spectralcosine_k-{k}_cluslabels.nii.gz',
        mask_res = 'results/tractmap/masks/sub-{subject}/brain_mask_dwi_resampled_tractmap.nii.gz'
    output:
        cluster_k_res = 'results/tractmap/resample_clus_seed/sub-{subject}/Native_space_from-{template}_seed-{seed}_method-spectralcosine_k-{k}_cluslabels_resampled.nii.gz'
    container: config['singularity_neuroglia']
    log: 'logs/resample_clus_seed/sub-{subject}_template-{template}_{seed}_k-{k}.log'
    group: 'tractmap'
    shell:
        'reg_resample -flo {input.cluster_k} -res {output.cluster_k_res} -ref {input.mask_res} -NN 0 &> {log}'

#Split cluster segmentations into binary masks
rule subj_split_clus_to_binary_masks:
    input:
        cluster_k_res = rules.resample_clus_seed.output.cluster_k_res
    params:
        mask_file = lambda wildcards,output: expand('{}/cluster{{kiter:02d}}_binaryMask.nii.gz'.format(output.cluster_k_splitdir),subject=wildcards.subject,kiter=range(0,int(wildcards.k)+1)),
        mask_bg = lambda wildcards,output: expand('{}/cluster00_binaryMask.nii.gz'.format(output.cluster_k_splitdir),subject=wildcards.subject) #we remove this file
    output:
        cluster_k_splitdir = directory('results/tractmap/clus_subj_seeds_binarized/sub-{subject}/sub-{subject}_{seed}_{template}_k-{k}_individual_seeds')
    container: config['singularity_neuroglia']
    log: 'logs/subj_split_clus_to_binary_masks/sub-{subject}_template-{template}_{seed}_k-{k}.log'
    group: 'tractmap'
    threads: 8
    shell:
        #use c3d's split command to go from discrete seg image to multiple binary images.. we remove image 00 since that is the background label
        'mkdir {output.cluster_k_splitdir} && c3d {input.cluster_k_res} -split -oo {params.mask_file} &>> {log}  && rm -f {params.mask_bg}'

#perform tracking from each cluster in subj space to get tract maps
# check bids-deriv dwi - tractography type?
# space-T1w, res-?
rule run_probtrack_from_clusters:
    input:
        cluster_k_splitdir = 'results/tractmap/clus_subj_seeds_binarized/sub-{subject}/sub-{subject}_{seed}_{template}_k-{k}_individual_seeds',
        mask = 'results/tractmap/masks/sub-{subject}/brain_mask_dwi.nii.gz',
        target_txt = rules.gen_targets_txt.output.target_txt,
        target_seg_dir = rules.split_targets.output.target_seg_dir
    params:
        seeds = lambda wildcards,input: expand('{}/cluster{{k_index:02d}}_binaryMask.nii.gz'.format(input.cluster_k_splitdir),subject=wildcards.subject,k_index=range(1,int(wildcards.k)+1)),
        bedpost_merged = 'results/Diffusion_7T/{subject}.bedpostX/merged',
        probtrack_opts = config['probtrack_tractmap']['opts'],
        nsamples = config['probtrack_tractmap']['nsamples'],
        container = config['singularity_cuda'],
        out_track_dirs = lambda wildcards,output: expand('{}/cluster{{k_index:02d}}'.format(output.probtrack_dir),k_index=range(1,int(wildcards.k)+1)),
        out_target_seg = lambda wildcards,output: expand('{}/cluster{{k_index:02d}}/seeds_to_{{target}}.nii.gz'.format(output.probtrack_dir),target=targets,k_index=range(1,int(wildcards.k)+1),allow_missing=True)
    output:
        probtrack_dir = directory('results/tractmap/probtrack_from_clusters/sub-{subject}/sub-{subject}_{seed}_{template}_k-{k}_probtrack')
    threads: 8
    resources:
        mem_mb = 32000,
        time = 30, #30 mins
        gpus = 1 #1 gpu
    log: 'logs/probtrack_from_clusters/sub-{subject}_template-{template}_{seed}_k-{k}.log'
    group: 'tractmap'
    shell:
        #Runs probtrackx2_gpu for each seed created from split cluster segmentations into binary masks
        'mkdir -p {params.out_track_dirs} && parallel --jobs 1 singularity exec -e --nv {params.container} '
        'probtrackx2_gpu --samples={params.bedpost_merged} --mask={input.mask} --seed={{1}} --seedref={{1}} --nsamples={params.nsamples} --dir={{2}} {params.probtrack_opts} -V 2 &> {log} '
        ' ::: {params.seeds} :::+ {params.out_track_dirs}'

# check bids-deriv dwi - tractography ?
# space-T1w, res-?
rule combine_tractmaps:
    input:
        probtrack_dir = rules.run_probtrack_from_clusters.output.probtrack_dir
    params:
        tractmaps = lambda wildcards,input: expand('{}/cluster{{k_index:02d}/fdt_paths.nii.gz'.format(input.probtrack_dir),k_index=range(1,int(wildcards.k)+1)),
    output:
        tractmaps_4d = 'results/tractmap/probtrack_from_clusters/sub-{subject}/sub-{subject}_NativeSpace_{seed}_{template}_k-{k}_tractmaps_4d.nii.gz'
    log: 'logs/combine_tractmaps/sub-{subject}_NativeSpace_template-{template}_{seed}_k-{k}.log'
    container: config['singularity_neuroglia']
    group: 'tractmap'
    resources:
        mem_mb = 32000
    shell:
        'fslmerge -t {output.tractmaps_4d} {params.tractmaps} &> {log}'

#transform tract maps back to template space
# space-{template}, tractogrpahy?
rule transform_tractmaps_to_template:
    input:
        probtrack_dir = rules.run_probtrack_from_clusters.output.probtrack_dir,
        affine =  config['ants_affine_mat'],
        warp =  config['ants_warp_nii'],
        ref = config['ants_ref_nii']
    params:
        in_tractmaps = lambda wildcards,input: expand('{}/cluster{{k_index:02d}}/fdt_paths.nii.gz'.format(input.probtrack_dir),k_index=range(1,int(wildcards.k)+1)),
        out_tractmaps = lambda wildcards,output: expand('{}/cluster{{k_index:02d}}_tractmap.nii.gz'.format(output.probtrack_dir),k_index=range(1,int(wildcards.k)+1))
    output:
        probtrack_dir = directory('results/tractmap/transform_tractmaps_to_template/sub-{subject}/sub-{subject}_{seed}_Space-{template}_k-{k}_probtrack')
    container: config['singularity_neuroglia']
    log: 'logs/transform_tractmaps_to_template/sub-{subject}_{seed}_Space-{template}_k-{k}.log'
    group: 'tractmap'
    threads: 8
    resources:
        mem_mb = 32000,
        time = 120
    shell:
        'mkdir -p {output} && ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 parallel  --jobs {threads} antsApplyTransforms -d 3 --interpolation Linear -i {{1}} -o {{2}}  -r {input.ref} -t {input.warp} -t {input.affine} &> {log} :::  {params.in_tractmaps} :::+ {params.out_tractmaps}' 

###If you want to have the individual clusters --> targets probtrack output, then need to develop this code, it is currently incomplete
#
#rule transform_conn_clusters_to_template:
#    input:
#        probtrack_dir = rules.run_probtrack_from_clusters.output.probtrack_dir,
#        affine =  config['ants_affine_mat'],
#        warp =  config['ants_warp_nii'],
#        ref = config['ants_ref_nii']
#    params:
#        in_connmap_clusters_3d = lambda wildcards,input: expand('{}/cluster{{k_index:02d}}/seeds_to_{{target}}.nii.gz'.format(input.probtrack_dir),subject=wildcards.subject,target=targets,k_index=range(1,int(wildcards.k)+1),allow_missing=True),
#        out_connmap_clusters_3d = lambda wildcards,output: expand('{}/cluster{{k_index:02d}}/seeds_to_{{target}}_space-{{template}}.nii.gz'.format(output.connmap_clusters_dir),subject=wildcards.subject,target=targets,template=config['template'],k_index=range(1,int(wildcards.k)+1),allow_missing=True),
#        clusterdirs = lambda wildcards: expand('cluster0{kran}',kran=range(1,int(wildcards.k)+1))
#    output:
#        connmap_clusters_dir = directory('results/tractmap/transform_connmap_clusters_to_template/sub-{subject}/sub-{subject}_{seed}_Space-{template}_k-{k}_probtrack')
#    container: config['singularity_neuroglia']
#    threads: 32
#    resources:
#        mem_mb = 128000
#    log: 'logs/transform_connmap_cluster_to_template/sub-{subject}_{seed}_{template}_k-{k}.log'
#    group: 'tractmap'
#    shell:
#        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 parallel  --jobs {threads} antsApplyTransforms -d 3 --interpolation Linear -i {{1}} -o {{2}}  -r {input.ref} -t {input.warp} -t {input.affine} &> {log} :::  {params.in_connmap_clusters_3d} :::+ {params.out_connmap_clusters_3d}'


# space-{template}, tractography 4d?
rule combine_tractmaps_warped:
    input:
        probtrack_dir = rules.transform_tractmaps_to_template.output.probtrack_dir
    params:
        tractmaps = lambda wildcards,input: expand('{}/cluster{{k_index:02d}}_tractmap.nii.gz'.format(input.probtrack_dir),k_index=range(1,int(wildcards.k)+1))
    output:
        tractmaps_4d = 'results/tractmap/transform_tractmaps_to_template/sub-{subject}/sub-{subject}_templateSpace-{template}_{seed}_k-{k}_tractmaps_4d.nii.gz'
    log: 'logs/combine_tractmaps_warped/sub-{subject}_templateSpace-{template}_{seed}_k-{k}.log'
    container: config['singularity_neuroglia']
    resources:
        mem_mb = 32000 
    group: 'tractmap'
    shell:
        'fslmerge -t {output.tractmaps_4d} {params.tractmaps} &> {log}'



#now, average the tract maps (over subjects) in template space
# space-template, desc-avg , tractography 4d?
rule avg_tractmaps_template:
    input: 
        tractmaps_4d = expand('results/tractmap/transform_tractmaps_to_template/sub-{subject}/sub-{subject}_templateSpace-{template}_{seed}_k-{k}_tractmaps_4d.nii.gz',subject=subjects,allow_missing=True)
    output: 
        average = 'results/tractmap/transform_tractmaps_to_template/templateSpace-{template}_{seed}_k-{k}_tractmaps_Averaged_4d.nii.gz'
    container: config['singularity_neuroglia']
    group: 'group_tractmap'
    threads: 8
    resources:
        mem_mb = 32000 
    shell:
        'AverageImages 4 {output} 0 {input}'
             

#use voting to get a discrete segmentation of tract map 
# space-template, desc-avgtractmap, dseg
rule vote_tractmap_template:
    input: 
        tractmaps = rules.avg_tractmaps_template.output.average
    params:
        bg_th = 100   # set only if avg streamline count > bg_th 
    output:
        vote_seg = 'results/tractmap/transform_tractmaps_to_template/templateSpace-{template}_{seed}_k-{k}_tractmaps_Averaged_dseg_FINAL.nii.gz'
    container: config['singularity_neuroglia']
    log: 'logs/vote_tractmap_template/{template}_{seed}_{k}.log'
    group: 'group_tractmap'
    threads: 8
    resources:
        mem_mb = 32000
    shell: 
        #get first vol, set it to bg_th value, this becomes image 0
        # then load all the tractmaps as images 1-k
        # then vote amongst the 0-k images, those where bg (bg_th) is highest get set to 0, otherwise set to the index (from 1-k)
        'c4d -verbose {input.tractmaps} -slice w 0:0 -threshold -inf +inf {params.bg_th} 0 {input.tractmaps} -slice w 0:-1  -vote -type uchar {output.vote_seg} &> {log}' 

