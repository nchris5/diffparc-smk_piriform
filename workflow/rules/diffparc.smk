#
#When lhrh_targets_independent_condition set to TRUE: Analyze connectivity of piriform segmentation within each hemisphere independently (considers each lh and rh region to be separate/different between hemispheres)
###Note that lh = 1-179 and rh = 180-358 and the values 110-179 and 290-358 were shifted 1 from initial positions for continuity of label values after the piriform 110 was removed

#When lhrh_targets_independent_condition set to FALSE: Analyze connectivity of piriform segmentation irrespective of hemispheres (consdiers each lh and rh region to be a combined single region across hemispheres)
###Note that lhANDrh both = 1-179 and the values 110-179 were shifted 1 from initial positions for continuity of label values after the piriform 110 was removed
if config['lhrh_targets_independent_condition']:
    rule combine_lr_hcp:
        input:
            lh = 'results/hcp_mmp/sub-{subject}/lh_removepir110_relabelled_1-179.native.hcp-mmp.nii.gz',
            rh = 'results/hcp_mmp/sub-{subject}/rh_removepir290_relabelled_180-358.native.hcp-mmp.nii.gz',
        output:
            lh_rh = 'results/diffparc/sub-{subject}/masks/lh_rh_removepir110_290_relabelled_1-358.native.hcp-mmp.nii.gz'
        threads: 8
	container: config['singularity_neuroglia']
        log: 'logs/combine_lr_hcp/{subject}.log'
        group: 'hcp_mmp_subj'
        shell: 'fslmaths {input.lh} -max {input.rh} {output.lh_rh} &> {log}'
else:
    rule combine_lr_hcp:
        input:
            lh = 'results/hcp_mmp/sub-{subject}/lh_removepir110_relabelled_1-179.native.hcp-mmp.nii.gz',
            rh = 'results/hcp_mmp/sub-{subject}/rh_removepir110_relabelled_1-179.native.hcp-mmp.nii.gz',
        output:
            lh_rh = 'results/diffparc/sub-{subject}/masks/lh_rh_removepir110_relabelled_1-179.native.hcp-mmp.nii.gz'
        threads: 8
	container: config['singularity_neuroglia']
        log: 'logs/combine_lr_hcp/{subject}.log'
        group: 'hcp_mmp_subj'
        shell: 'fslmaths {input.lh} -max {input.rh} {output.lh_rh} &> {log}'


###Get and binarize probabilistic HCP7TTemplate Piriform seed (275x327x275, 0.7mm iso). Note, the binarized seed is not used until much later.
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
    output: 'results/template_masks/sub-{template}_desc-{seed}_thr_binarized.nii.gz'
    container: config['singularity_neuroglia']
    log: 'logs/binarize_template_seed/{template}_{seed}.log'
    group: 'group1'
    shell:
        'fslmaths {input} -thr {params.threshold} -bin {output} &> {log}'



###This HCP7T_DWI_BRAINMASK_SPACE will (173x207x173, 1.05mm iso) determine the resolution to perform probtrack
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


###TRANSFORM AND RESAMPLE PIRIFORM SEED AND HCP-MMP TARGETS TO HCP7T_DWI_SPACE_AND_RESOLUTION
#Transform HCP7TTemplate probabilistic seed (275x327x275, 0.7mm iso) to HCP7TSubjectSpace (256x256x256, 1.0mm iso) then resample to HCP7TSubjectSpace_dwi_brainmask resolution (173x207x173, 1.05mm iso)
rule transform_seed_to_subject:
    input:
        seed = rules.get_template_seed.output,
        affine =  config['ants_affine_mat'],
        invwarp =  config['ants_invwarp_nii'],
        ref = rules.combine_lr_hcp.output.lh_rh
    output: 'results/diffparc/sub-{subject}/masks/seed_from-{template}_{seed}_space-HCPnative.nii.gz'
    envmodules: 'ants'
    container: config['singularity_neuroglia']
    log: 'logs/transform_seed_to_subject/{template}_sub-{subject}_{seed}.log'
    group: 'hcp_mmp_subj'
    threads: 8
    shell:
        #Linear interp here now, since probabilistic seg
        'antsApplyTransforms -d 3 --interpolation Linear -i {input.seed} -o {output} -r {input.ref} -t [{input.affine},1] -t {input.invwarp} &> {log}'

rule resample_subject_seed:
    input:    
        seed = rules.transform_seed_to_subject.output,
        mask = rules.binarize_dwi_brainmask.output.mask
    output:
        seed_res = 'results/diffparc/sub-{subject}/masks/seed_from-{template}_{seed}_space-HCPnative_resampled_dwi_resolution.nii.gz',
    envmodules: 'ants'
    container: config['singularity_neuroglia']
    log: 'logs/resample_seed/{template}_sub-{subject}_{seed}.log'
    group: 'hcp_mmp_subj'
    shell:
        #Linear interp here now, since probabilistic seg
        'antsApplyTransforms -d 3 --interpolation Linear -i {input.seed} -o {output.seed_res} -r {input.mask} &> {log}' 

#Binarize seed in HCP7TSubjectSpace_dwi_brainmask resolution (173x207x173, 1.05mm iso)
rule binarize_subject_seed:
    input:
        seed_res = rules.resample_subject_seed.output
    params:
        threshold = config['prob_seg_threshold']
    output: 'results/diffparc/sub-{subject}/masks/seed_from-{template}_{seed}_space-HCPnative_resampled_dwi_resolution_thr_binarized.nii.gz'
    container: config['singularity_neuroglia']
    log: 'logs/binarize_subject_seed/{template}_sub-{subject}_{seed}.log'
    group: 'hcp_mmp_subj'
    shell:
        'fslmaths {input} -thr {params.threshold} -bin {output} &> {log}'

#Transform T1w space targets created in rules/hcp_mmk.smk (256x256x256, 1.0mm iso) to HCP7TSubjectSpace_dwi_brainmask resolution (173x207x173, 1.05mm iso)
if config['lhrh_targets_independent_condition']:
    rule resample_targets:
        input:
            targets = 'results/diffparc/sub-{subject}/masks/lh_rh_removepir110_290_relabelled_1-358.native.hcp-mmp.nii.gz',
            mask = rules.binarize_dwi_brainmask.output.mask,
        output:
            targets_res = 'results/diffparc/sub-{subject}/masks/lh_rh_targets_from-native_label-hcpmmp_space-native_resampled_dwi_resolution.nii.gz'
        envmodules: 'ants'
        container: config['singularity_neuroglia']
        log: 'logs/resample_targets/sub-{subject}.log'
        group: 'hcp_mmp_subj'
        shell:
            #Nearest interp since targets are labels
            'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.targets} -o {output.targets_res} -r {input.mask} &> {log}'
else:
    rule resample_targets:
        input:
            targets = 'results/diffparc/sub-{subject}/masks/lh_rh_removepir110_relabelled_1-179.native.hcp-mmp.nii.gz',
            mask = rules.binarize_dwi_brainmask.output.mask,
        output:
            targets_res = 'results/diffparc/sub-{subject}/masks/lh_rh_targets_from-native_label-hcpmmp_space-native_resampled_dwi_resolution.nii.gz'
        envmodules: 'ants'
        container: config['singularity_neuroglia']
        log: 'logs/resample_targets/sub-{subject}.log'
        group: 'hcp_mmp_subj'
        shell:
            #Nearest interp since targets are labels
            'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.targets} -o {output.targets_res} -r {input.mask} &> {log}'


#The piriform seed overlaps adjacent atlas targets, so in subject space the target voxels that overlapped seed voxels are removed. First the binarized HCP7TSubjectSpace_dwi_resolution piriform seed is given a large value (1000) and then the maximum between the seed and HCP7TSubjectSpace_dwi_resolution hcp-mmp targets is taken and the piriform is then removed from the targets.
rule remove_piriform_target_overlap:
    input:
        seed_res_bin = rules.binarize_subject_seed.output,
        targets_res = rules.resample_targets.output.targets_res
    output:
        targets_pir_overlap_removed = 'results/diffparc/sub-{subject}/masks/lh_rh_targets_from-{template}_{seed}_space-native_resampled_dwi_resolution_pir-overlap-removed.nii.gz'
    container: config['singularity_neuroglia']
    log: 'logs/remove_piriform_target_overlap/sub-{subject}_{template}_{seed}.log'
    group: 'hcp_mmp_subj'
    shell:
        'fslmaths {input.seed_res_bin} -mul 1000 -max {input.targets_res} -uthr 999 {output.targets_pir_overlap_removed} &> {log}'
        

#Split HCP7TSubjectSpace_dwi_resolution hcp-mmp targets (173x207x173, 1.05mm iso) into individual files for each target
#Note that in Ali's diffparc the target nums is [str(i) for i in range(len(targets))], which skips V1 since it has a value of 1 and this numbering starts at 0 in this form
rule split_targets:
    input:
        targets_res = rules.remove_piriform_target_overlap.output.targets_pir_overlap_removed,
    params:
        target_nums = lambda wildcards: [str(i) for i in range(1,len(targets)+1)],
        target_seg = expand('results/diffparc/sub-{subject}/lh_rh_targets_from-{template}_{seed}_space-native_resampled_dwi_resolution_pir-overlap-removed/{target}.nii.gz',target=targets,allow_missing=True)
    output:
        target_seg_dir = directory('results/diffparc/sub-{subject}/lh_rh_targets_from-{template}_{seed}_space-native_resampled_dwi_resolution_pir-overlap-removed')
    container: config['singularity_neuroglia']
    log: 'logs/split_targets/sub-{subject}_{template}_{seed}.log'
    threads: 32
    group: 'hcp_mmp_subj'
    shell:
        'mkdir -p {output} && parallel  --jobs {threads} fslmaths {input.targets_res} -thr {{1}} -uthr {{1}} -bin {{2}} &> {log} ::: {params.target_nums} :::+ {params.target_seg}'


#Generate a text file listing the hcp-mmp target names as a column vector
rule gen_targets_txt:
    input:
        target_seg_dir = rules.split_targets.output.target_seg_dir
    params:
        target_seg = expand('results/diffparc/sub-{subject}/lh_rh_targets_from-{template}_{seed}_space-native_resampled_dwi_resolution_pir-overlap-removed/{target}.nii.gz',target=targets,allow_missing=True)
    output:
        target_txt = 'results/diffparc/sub-{subject}/target_images_from-{template}_{seed}_overlap_removed.txt'
    log: 'logs/get_targets_txt/sub-{subject}_{template}_{seed}.log'
    group: 'hcp_mmp_subj'
    run:
        f = open(output.target_txt,'w')
        for s in params.target_seg:
            f.write(f'{s}\n')
        f.close()



###RUN PROBABILISTIC TRACTOGRAPHY
#If using 378 regions, use probtrackx2 due to memory issues encountered with gpu at this number of targets (submitted as per subject group jobs with everything previous); If using 179 regions, use probtrackx2_gpu (submitted as per subject jobs independent of everything previous)
#Piriform seed (seed_res), hcp-mmp targets (target_seg_dir), and dwi brainmask (mask) are at the dwi size and resolution (i.e. 173x207x173, 1.05mm iso)
if config['lhrh_targets_independent_condition']:
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
            container = config['singularity_neuroglia']
        output:
            probtrack_dir = directory('results/diffparc/sub-{subject}/probtrack_{template}_{seed}')
        threads: 16
        resources:
            mem_mb = 64000,
            time = 180 #3hrs if running probtrackx2 without gpu (test ~30 mins so 3hr for safety)
        log: 'logs/run_probtrack/{template}_sub-{subject}_{seed}.log'
        group: 'run_probtrack'
        shell:
            'mkdir -p {output.probtrack_dir} && singularity exec -e --nv {params.container} probtrackx2 --samples={params.bedpost_merged} --mask={input.mask} --seed={input.seed_res} --targetmasks={input.target_txt} --seedref={input.seed_res} --nsamples={params.nsamples} --dir={output.probtrack_dir} {params.probtrack_opts} -V 1  &> {log}'
else:
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
            time = 30, #30 mins if running probtrackx2_gpu
            gpus = 1 #1 gpu if running probtrackx2_gpu
        log: 'logs/run_probtrack/{template}_sub-{subject}_{seed}.log'
        group: 'run_probtrack'
        shell:
            'mkdir -p {output.probtrack_dir} && singularity exec -e --nv {params.container} probtrackx2_gpu --samples={params.bedpost_merged} --mask={input.mask} --seed={input.seed_res} --targetmasks={input.target_txt} --seedref={input.seed_res} --nsamples={params.nsamples} --dir={output.probtrack_dir} {params.probtrack_opts} -V 2  &> {log}' 



###TRANSFORM PIRIFORM CONNECTIVITY MAPS (173x207x173, 1.05mm iso) BACK TO HIGH RES TEMPLATE SPACE (275x327x275; 0.7mm iso) AND THEN PERFORM SPECTRAL CLUSTERING
rule transform_conn_to_template:
    input:
        probtrack_dir = rules.run_probtrack.output.probtrack_dir,
        affine =  config['ants_affine_mat'],
        warp =  config['ants_warp_nii'],
        ref = config['ants_ref_nii']
    params:
        in_connmap_3d = expand('results/diffparc/sub-{subject}/probtrack_{template}_{seed}/seeds_to_{target}.nii.gz',target=targets,allow_missing=True),
        out_connmap_3d = expand('results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped/seeds_to_{target}_space-{template}.nii.gz',target=targets,allow_missing=True),
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


#Save connmaps in template space for each subject as npz files
rule save_connmap_template_npz:
    input:
        mask = 'results/template_masks/sub-{template}_desc-{seed}_thr_binarized.nii.gz', #HCP_Template_Piriform_Mask
        probtrack_dir = 'results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped',
        dependrun = rules.transform_conn_to_template.output.connmap_dir
    params:
        connmap_3d = expand('results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped/seeds_to_{target}_space-{template}.nii.gz',target=targets,allow_missing=True),
    output:
        connmap_npz = 'results/diffparc/sub-{subject}/connmap/sub-{subject}_space-{template}_seed-{seed}_connMap.npz'
    log: 'logs/save_connmap_to_template_npz/sub-{subject}_{seed}_{template}.log'
    threads: 8
    resources:
        mem_mb = 32000
    group: 'post_track'
    conda: '../envs/sklearn.yml'
    script: '../scripts/save_connmap_template_npz.py'


#Group all subjects connmaps into a single npz file
rule gather_connmap_group:
    input:
        connmap_npz = expand('results/diffparc/sub-{subject}/connmap/sub-{subject}_space-{template}_seed-{seed}_connMap.npz',subject=subjects,allow_missing=True),
        dependrun = expand(rules.save_connmap_template_npz.output.connmap_npz,subject=subjects,allow_missing=True)
    output:
        connmap_group_npz = 'results/connmap/group_space-{template}_seed-{seed}_connMap.npz'
    log: 'logs/gather_connmap_group/{seed}_{template}.log'
    threads: 8
    resources:
        mem_mb = 32000
    conda: '../envs/sklearn.yml'
    group: 'group_post_track'
    script: '../scripts/gather_connmap_group.py'


     
rule spectral_clustering:
    input:
        connmap_group_npz = 'results/connmap/group_space-{template}_seed-{seed}_connMap.npz',
        dependrun = rules.gather_connmap_group.output.connmap_group_npz
    params:
        max_k = config['max_k']
    output:
        cluster_k = expand('results/clustering/group_space-{template}_seed-{seed}_method-spectralcosine_k-{k}_cluslabels.nii.gz',k=range(2,config['max_k']+1),allow_missing=True)
    log: 'logs/spectral_clustering/{seed}_{template}.log'
    threads: 8
    resources:
        mem_mb = 32000
    conda: '../envs/sklearn.yml'
    group: 'group_post_track'
    script: '../scripts/spectral_clustering.py'



###SAVE PIRIFORM CONNECTIVITY FOR EACH SPECIFIC CLUSTER DETERMINED AND PLOT CONNECTIVITY OF EACH CLUSTER
#After switching to 358 regions, this is very computationally expensive, so greatly increased the resources for this
rule save_connmapClusters_template_npz:
    input:
        probtrack_dir = 'results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped',
        clusters_template_space = expand('results/clustering/group_space-{template}_seed-{seed}_method-spectralcosine_k-{k}_cluslabels.nii.gz',k=range(2,config['max_k']+1),allow_missing=True),
    params:
        connmap_3d = expand('results/diffparc/sub-{subject}/probtrack_{template}_{seed}_warped/seeds_to_{target}_space-{template}.nii.gz',target=targets,allow_missing=True),
        max_k = config['max_k'],
        connmap_Clusters_npz = directory('results/connmapClusters/sub-{subject}/connmap_{seed}_{template}/sub-{subject}_space-{template}_seed-{seed}')
    output:
        connmap_Clusters_npz_depend = directory(expand('results/connmapClusters/sub-{subject}/connmap_{seed}_{template}/sub-{subject}_space-{template}_seed-{seed}_k-{k}',k=range(2,config['max_k']+1),allow_missing=True))
    log: 'logs/save_connmapClusters_template_npz/sub-{subject}_{seed}_{template}.log'
    conda: '../envs/sklearn.yml'
    threads: 16
    resources:
        mem_mb = 64000,
	time = 180
    group: 'post_track_clust'
    script: '../scripts/save_connmapClusters_template_npz.py'


if config['lhrh_targets_independent_condition']:
    rule gather_connmapClusters_group:
        input:
            connmap_Clusters_npz_dir = expand('results/connmapClusters/sub-{subject}/connmap_{seed}_{template}/sub-{subject}_space-{template}_seed-{seed}_k-{k}',subject=subjects,k=range(2,config['max_k']+1),allow_missing=True),
        params:
            max_k = config['max_k'],
            connmap_Clusters_group_npz_dir_base = directory('results/connmapClusters/connmapGroup_{seed}_{template}'),
            num_targets = 358
        output:
            connmap_Clusters_group_npz_dir = directory(expand('results/connmapClusters/connmapGroup_{seed}_{template}_k-{k}',k=range(2,config['max_k']+1),allow_missing=True))
        log: 'logs/gather_connmapClusters_group/{seed}_{template}.log'
        conda: '../envs/sklearn.yml'
        threads: 16
        resources:
            mem_mb = 64000,
            time = 180
        group: 'group_post_track_clust'
        script: '../scripts/gather_connmapClusters_group.py'
else:
    rule gather_connmapClusters_group:
        input:
            connmap_Clusters_npz_dir = expand('results/connmapClusters/sub-{subject}/connmap_{seed}_{template}/sub-{subject}_space-{template}_seed-{seed}_k-{k}',subject=subjects,k=range(2,config['max_k']+1),allow_missing=True),
        params:
            max_k = config['max_k'],
            connmap_Clusters_group_npz_dir_base = directory('results/connmapClusters/connmapGroup_{seed}_{template}'),
            num_targets = 179
        output:
            connmap_Clusters_group_npz_dir = directory(expand('results/connmapClusters/connmapGroup_{seed}_{template}_k-{k}',k=range(2,config['max_k']+1),allow_missing=True))
        log: 'logs/gather_connmapClusters_group/{seed}_{template}.log'
        conda: '../envs/sklearn.yml'
        threads: 16
        resources:
            mem_mb = 64000,
            time = 180
        group: 'group_post_track_clust'
        script: '../scripts/gather_connmapClusters_group.py' 



###CREATE INPUT NODES AND EDGES TABLES TO BE USED IN GEPHI FOR CONNECTIVITY NETWORK VISUALIZATION/ANALYSIS IN GEPHI
if config['lhrh_targets_independent_condition']:
    rule create_gephi_input_nodes_and_edge_tables:
        input:
            connmap_Clusters_group_dir = expand('results/connmapClusters/connmapGroup_{seed}_{template}_k-{k}',template=config['template'],k=range(2,config['max_k']+1),allow_missing=True),
        params:
            seeds_py = list(config['template_prob_seg'].keys()),
            max_k = config['max_k'],
            connmapClusters_dir = 'results/connmapClusters',
            num_targets = 358,
            GlasserTable_22regions_colours = 'resources/Glasser_2016_Table_BoldSections_Manually_Identified_removePiriform_LRseparate_380_Regions_relabelled_L1-179_R180-378_convcsv.csv',
        log: 'logs/create_gephi_input_nodes_and_edge_tables/Gephi_Input_NodesandEdges_{seed}.log'
        conda: '../envs/sklearn.yml'
        threads: 16
        resources:
            mem_mb = 64000,
            time = 180
        group: 'groupGephipost_track'
        output:
            Gephi_Nodes_and_Edges_Tables_dir = directory('results/gephi_input_NodesandEdges_tables_{seed}')
        script: '../scripts/create_connmapClusters_group_GephiNodes_and_Edges_Tables.py'
else:
    rule create_gephi_input_nodes_and_edge_tables:
        input:
            connmap_Clusters_group_dir = expand('results/connmapClusters/connmapGroup_{seed}_{template}_k-{k}',template=config['template'],k=range(2,config['max_k']+1),allow_missing=True),
        params:
            seeds_py = list(config['template_prob_seg'].keys()),
            max_k = config['max_k'],
            connmapClusters_dir = 'results/connmapClusters',
            num_targets = 179,
            GlasserTable_22regions_colours = 'resources/Glasser_2016_Table_BoldSections_Manually_Identified_removePiriform110_convcsv.csv'
        log: 'logs/create_gephi_input_nodes_and_edge_tables/Gephi_Input_NodesandEdges_{seed}.log'
        conda: '../envs/sklearn.yml'
        threads: 16
        resources:
            mem_mb = 64000,
            time = 180
        group: 'groupGephipost_track_clust'
        output:
            Gephi_Nodes_and_Edges_Tables_dir = directory('results/gephi_input_NodesandEdges_tables_{seed}')
        script: '../scripts/create_connmapClusters_group_GephiNodes_and_Edges_Tables.py'

