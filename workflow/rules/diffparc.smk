
rule combine_lr_hcp:
    input:
        lh = 'results/hcp_mmp/sub-{subject}/lh.native.hcp-mmp.nii.gz',
        rh = 'results/hcp_mmp/sub-{subject}/rh.native.hcp-mmp.nii.gz'
    output:
        lh_rh = 'results/diffparc/sub-{subject}/masks/lh_rh.native.hcp-mmp.nii.gz'
    container: config['singularity_neuroglia']
    log: 'logs/combine_lr_hcp/{subject}.log'
    group: 'hcp_mmp_subj'
    shell: 'fslmaths {input.lh} -max {input.rh} {output.lh_rh} &> {log}'


rule get_template_seed:
    input: 
        seed = lambda wildcards: config['template_prob_seg'][wildcards.seed],
    output: 'results/template_masks/sub-{template}_desc-{seed}_probseg.nii.gz'
    container: config['singularity_neuroglia']
    log: 'logs/get_template_seed/{template}_{seed}.log'
    group: 'group0'
    shell:
        'cp -v {input} {output} &> {log}'
 

rule binarize_template_seed:
    input: rules.get_template_seed.output
    params:
        threshold = config['prob_seg_threshold']
    output: 'results/template_masks/sub-{template}_desc-{seed}_thr02_binarized.nii.gz'
    container: config['singularity_neuroglia']
    log: 'logs/binarize_template_seed/{template}_{seed}.log'
    group: 'group0'
    shell:
        'fslmaths {input} -thr {params.threshold} -bin {output} &> {log}'


#transform probabilistic seed to subject
rule transform_to_subject:
    input: 
        seed = rules.get_template_seed.output,
        affine =  config['ants_affine_mat'],
        invwarp =  config['ants_invwarp_nii'],
        ref = rules.combine_lr_hcp.output.lh_rh
    output: 'results/diffparc/sub-{subject}/masks/seed_from-{template}_{seed}.nii.gz'
    envmodules: 'ants'
    container: config['singularity_neuroglia']
    log: 'logs/transform_to_subject/{template}_sub-{subject}_{seed}.log'
    group: 'hcp_mmp_subj'
    threads: 8
    shell:
        #Spline interp here now, since probabilistic seg
        'antsApplyTransforms -d 3 --interpolation BSpline -i {input.seed} -o {output} -r {input.ref} -t [{input.affine},1] -t {input.invwarp} &> {log}'


###Use this version if you would like to perform probtrack at a resolution different than the HCP_7T_Diffusion data as provided by probtrack seed_resolution in the config file
#rule binarize_dwi_brainmask:
#    input:
#        dwi = 'results/Diffusion_7T/{subject}.bedpostX/mean_S0samples.nii.gz'
#    params:
#        seed_resolution = config['probtrack']['seed_resolution']
#    output:
#        mask = 'results/diffparc/sub-{subject}/masks/brain_mask_dwi.nii.gz',
#        mask_res = 'results/diffparc/sub-{subject}/masks/brain_mask_resampled_dwi_resolution.nii.gz'
#    container: config['singularity_neuroglia']
#    log: 'logs/resample_brainmask/sub-{subject}.log'
#    group: 'hcp_mmp_subj'
#    shell:
#        'fslmaths {input.dwi} -bin {output.mask} && '
#        'mri_convert {output.mask} -vs {params.seed_resolution} {params.seed_resolution} {params.seed_resolution} {output.mask_res} -rt nearest &> {log}'

#rule resample_targets:
#    input:
#        mask = rules.binarize_dwi_brainmask.output.mask_res,
#        targets = rules.combine_lr_hcp.output.lh_rh
#    output:
#        targets_res = 'results/diffparc/sub-{subject}/masks/lh_rh_targets_resampled_dwi_resolution.nii.gz'
#    envmodules: 'ants'
#    container: config['singularity_neuroglia']
#    log: 'logs/resample_targets/sub-{subject}.log'
#    group: 'hcp_mmp_subj'
#    shell:
#        #Nearest interp since targets are labels
#        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.targets} -o {output.targets_res} -r {input.mask} &> {log}'

#rule resample_seed:
#    input:
#        seed = rules.transform_to_subject.output,
#        mask = rules.binarize_dwi_brainmask.output.mask_res
#    output:
#        seed_res = 'results/diffparc/sub-{subject}/masks/seed_from-{template}_{seed}_resampled_dwi_resolution.nii.gz',
#    envmodules: 'ants'
#    container: config['singularity_neuroglia']
#    log: 'logs/resample_seed/{template}_sub-{subject}_{seed}.log'
#    group: 'hcp_mmp_subj'
#    shell:
#        #Spline interp here now, since probabilistic seg
#        'antsApplyTransforms -d 3 --interpolation BSpline -i {input.seed} -o {output.seed_res} -r {input.mask} &> {log}'
###End of commented version if wanting probtrack at resolution different than HCP_7T_Diffusion data


###DEFAULT Use this version if want to perform probtrack at the HCP_7T_Diffusion resolution
rule binarize_dwi_brainmask:
    input:
        dwi = 'results/Diffusion_7T/{subject}.bedpostX/mean_S0samples.nii.gz'
    output:
        mask = 'results/diffparc/sub-{subject}/masks/brain_mask_dwi.nii.gz',
    container: config['singularity_neuroglia']
    log: 'logs/binarize_dwi_brainmask/sub-{subject}.log'
    group: 'hcp_mmp_subj'
    shell:
        'fslmaths {input.dwi} -bin {output.mask} &> {log}'

rule resample_targets:
    input:
        mask = rules.binarize_dwi_brainmask.output.mask,
        targets = rules.combine_lr_hcp.output.lh_rh
    output:
        targets_res = 'results/diffparc/sub-{subject}/masks/lh_rh_targets_resampled_dwi_resolution.nii.gz'
    envmodules: 'ants'
    container: config['singularity_neuroglia']
    log: 'logs/resample_targets/sub-{subject}.log'
    group: 'hcp_mmp_subj'
    shell:
        #Nearest interp since targets are labels
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.targets} -o {output.targets_res} -r {input.mask} &> {log}'


rule resample_seed:
    input: 
        seed = rules.transform_to_subject.output,
        mask = rules.binarize_dwi_brainmask.output.mask
    output:
        seed_res = 'results/diffparc/sub-{subject}/masks/seed_from-{template}_{seed}_resampled_dwi_resolution.nii.gz',
    envmodules: 'ants'
    container: config['singularity_neuroglia']
    log: 'logs/resample_seed/{template}_sub-{subject}_{seed}.log'
    group: 'hcp_mmp_subj'
    shell:
        #Spline interp here now, since probabilistic seg
        'antsApplyTransforms -d 3 --interpolation BSpline -i {input.seed} -o {output.seed_res} -r {input.mask} &> {log}'
###End of DEFAULT version if want to perofrm probtrack at the HCP_7T_Diffusion resolution


rule binarize_subject_seed:
    input: 
        seed_res = rules.resample_seed.output.seed_res
    params:
        threshold = config['prob_seg_threshold']
    output: 'results/diffparc/sub-{subject}/masks/seed_from-{template}_{seed}_resampled_dwi_resolution_thr02_binarized.nii.gz'
    container: config['singularity_neuroglia']
    log: 'logs/binarize_subject_seed/{template}_sub-{subject}_{seed}.log'
    group: 'hcp_mmp_subj'
    shell:
        'fslmaths {input} -thr {params.threshold} -bin {output} &> {log}'
         


rule split_targets:
    input: 
        targets_res = 'results/diffparc/sub-{subject}/masks/lh_rh_targets_resampled_dwi_resolution.nii.gz',
    params:
        target_nums = lambda wildcards: [str(i) for i in range(len(targets))],
        target_seg = expand('results/diffparc/sub-{subject}/targets/{target}.nii.gz',target=targets,allow_missing=True)
    output:
        target_seg_dir = directory('results/diffparc/sub-{subject}/targets')
    container: config['singularity_neuroglia']
    log: 'logs/split_targets/sub-{subject}.log'
    threads: 32 
    group: 'hcp_mmp_subj'
    shell:
        'mkdir -p {output} && parallel  --jobs {threads} fslmaths {input.targets_res} -thr {{1}} -uthr {{1}} -bin {{2}} &> {log} ::: {params.target_nums} :::+ {params.target_seg}'

rule gen_targets_txt:
    input:
        target_seg_dir = rules.split_targets.output.target_seg_dir
    params:
        target_seg = expand('results/diffparc/sub-{subject}/targets/{target}.nii.gz',target=targets,allow_missing=True)
    output:
        target_txt = 'results/diffparc/sub-{subject}/target_images.txt'
    log: 'logs/get_targets_txt/sub-{subject}.log'
    group: 'hcp_mmp_subj'
    run:
        f = open(output.target_txt,'w')
        for s in params.target_seg:
            f.write(f'{s}\n')
        f.close()



rule run_probtrack:
    input:
        seed_res = rules.binarize_subject_seed.output,
        target_txt = rules.gen_targets_txt.output.target_txt,
        mask = rules.binarize_dwi_brainmask.output.mask,
        target_seg_dir = rules.split_targets.output.target_seg_dir
    params:
        bedpost_merged = 'results/Diffusion_7T/{subject}.bedpostX/merged',
        probtrack_opts = config['probtrack']['opts'],
        out_target_seg = expand('results/diffparc/sub-{subject}/probtrack_{template}_{seed}/seeds_to_{target}.nii.gz',target=targets,allow_missing=True),
        nsamples = config['probtrack']['nsamples'],
        container = config['singularity_cuda']
    output:
        probtrack_dir = directory('results/diffparc/sub-{subject}/probtrack_{template}_{seed}')
    threads: 8
    resources: 
        mem_mb = 32000, 
        time = 30, #30 mins
        gpus = 1 #1 gpu
    log: 'logs/run_probtrack/{template}_sub-{subject}_{seed}.log'
    shell:
        'mkdir -p {output.probtrack_dir} && singularity exec -e --nv {params.container} probtrackx2_gpu --samples={params.bedpost_merged} --mask={input.mask} --seed={input.seed_res} --targetmasks={input.target_txt} --seedref={input.seed_res} --nsamples={params.nsamples} --dir={output.probtrack_dir} {params.probtrack_opts} -V 2  &> {log}' 



rule transform_conn_to_template:
    input:
        probtrack_dir = rules.run_probtrack.output.probtrack_dir,
        affine =  config['ants_affine_mat'],
        warp =  config['ants_warp_nii'],
        ref = config['ants_ref_nii']
    params:
        in_connmap_3d = expand('results/diffparc/sub-{subject}/probtrack_{template}_{seed}/seeds_to_{target}.nii.gz',target=targets,allow_missing=True),
        out_connmap_3d = expand('results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped/seeds_to_{target}_space-{template}.nii.gz',target=targets,allow_missing=True)
    output:
        connmap_dir = directory('results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped')
    envmodules: 'ants'
    container: config['singularity_neuroglia']
    threads: 32
    resources:
        mem_mb = 128000
    log: 'logs/transform_conn_to_template/sub-{subject}_{seed}_{template}.log'
    group: 'post_track'
    shell:
        'mkdir -p {output} && ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1 parallel  --jobs {threads} antsApplyTransforms -d 3 --interpolation Linear -i {{1}} -o {{2}}  -r {input.ref} -t {input.warp} -t {input.affine} &> {log} :::  {params.in_connmap_3d} :::+ {params.out_connmap_3d}' 



rule save_connmap_template_npz:
    input:
        mask = 'results/template_masks/sub-{template}_desc-{seed}_thr02_binarized.nii.gz', #HCP_Template_Piriform_Mask
        probtrack_dir = 'results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped'
    params:
        connmap_3d = expand('results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped/seeds_to_{target}_space-{template}.nii.gz',target=targets,allow_missing=True),
    output:
        connmap_npz = 'results/diffparc/sub-{subject}/connmap/sub-{subject}_space-{template}_seed-{seed}_connMap.npz'
    log: 'logs/save_connmap_to_template_npz/sub-{subject}_{seed}_{template}.log'
    group: 'post_track'
    conda: '../envs/sklearn.yml'
    script: '../scripts/save_connmap_template_npz.py'



rule gather_connmap_group:
    input:
        connmap_npz = expand('results/diffparc/sub-{subject}/connmap/sub-{subject}_space-{template}_seed-{seed}_connMap.npz',subject=subjects,allow_missing=True)
    output:
        connmap_group_npz = 'results/connmap/group_space-{template}_seed-{seed}_connMap.npz'
    log: 'logs/gather_connmap_group/{seed}_{template}.log'
    conda: '../envs/sklearn.yml'
    group: 'group_post_track'
    script: '../scripts/gather_connmap_group.py'
     
rule spectral_clustering:
    input:
        connmap_group_npz = 'results/connmap/group_space-{template}_seed-{seed}_connMap.npz'
    params:
        max_k = config['max_k']
    output:
        cluster_k = expand('results/clustering/group_space-{template}_seed-{seed}_method-spectralcosine_k-{k}_cluslabels.nii.gz',k=range(2,config['max_k']+1),allow_missing=True)
    log: 'logs/spectral_clustering/{seed}_{template}.log'
    conda: '../envs/sklearn.yml'
    group: 'group_post_track'
    script: '../scripts/spectral_clustering.py'



rule save_connmapClusters_template_npz:
    input:
        probtrack_dir = 'results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped',
        clusters_template_space = expand('results/clustering/group_space-{template}_seed-{seed}_method-spectralcosine_k-{k}_cluslabels.nii.gz',k=range(2,config['max_k']+1),allow_missing=True)
    params:
        connmap_3d = expand('results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped/seeds_to_{target}_space-{template}.nii.gz',target=targets,allow_missing=True),
        max_k = config['max_k'],
        connmap_Clusters_npz = directory('results/connmapClusters/sub-{subject}/connmap_{seed}_{template}/sub-{subject}_space-{template}_seed-{seed}')
    output:
        connmap_Clusters_npz_depend = directory(expand('results/connmapClusters/sub-{subject}/connmap_{seed}_{template}/sub-{subject}_space-{template}_seed-{seed}_k-{k}',k=range(2,config['max_k']+1),allow_missing=True))
    log: 'logs/save_connmapClusters_template_npz/sub-{subject}_{seed}_{template}.log'
    conda: '../envs/sklearn.yml'
    group: 'post_track'
    script: '../scripts/save_connmapClusters_template_npz.py'

rule gather_connmapClusters_group:
    input:
        connmap_Clusters_npz_dir = expand('results/connmapClusters/sub-{subject}/connmap_{seed}_{template}/sub-{subject}_space-{template}_seed-{seed}_k-{k}',subject=subjects,k=range(2,config['max_k']+1),allow_missing=True)
    params:
        max_k = config['max_k'],
        connmap_Clusters_group_npz_dir_base = directory('results/connmapClusters/connmapGroup_{seed}_{template}')
    output:
        connmap_Clusters_group_npz_dir = directory(expand('results/connmapClusters/connmapGroup_{seed}_{template}_k-{k}',k=range(2,config['max_k']+1),allow_missing=True))
    log: 'logs/gather_connmapClusters_group/{seed}_{template}.log'
    conda: '../envs/sklearn.yml'
    group: 'groupC_post_track'
    script: '../scripts/gather_connmapClusters_group.py'



rule create_gephi_input_nodes_and_edge_tables:
    input:
        connmap_Clusters_group_dir = expand('results/connmapClusters/connmapGroup_{seed}_{template}_k-{k}',seed=seeds,template=config['template'],k=range(2,config['max_k']+1))
    params:
        seeds_py = list(config['template_prob_seg'].keys()),
        max_k = config['max_k'],
        connmapClusters_dir = 'results/connmapClusters',
        GlasserTable_22regions_colours = 'resources/Glasser_2016_Table_BoldSections_Manually_Identified_convcsv.csv'
    log: 'logs/create_gephi_input_nodes_and_edge_tables/Gephi_Input_NodesandEdges.log'
    conda: '../envs/sklearn.yml'
    group: 'groupGephipost_track'
    output:
        Gephi_Nodes_and_Edges_Tables_dir = directory('results/gephi_input_NodesandEdges_tables')
    script: '../scripts/create_connmapClusters_group_GephiNodes_and_Edges_Tables.py'


