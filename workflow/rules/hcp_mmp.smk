#
#rule all_hcpmmk:
#    input: expand('results/hcp_mmp/sub-{subject}/{hemi}.native.hcp-mmp.nii.gz',subject=subjects,hemi=hemis)

hemi = hemis

rule convert_to_gifti:
    input: join(config['in_freesurfer'],'surf','{hemi}.{surfname}')
    output: 'results/hcp_mmp/sub-{subject}/{hemi}_space-fsaverage.{surfname}.surf.gii'
    params: 
        license = config['fs_license']
    container: config['singularity_freesurfer']
    log: 'logs/convert_to_gifti/sub-{subject}_{hemi}_{surfname}.log'
    group: 'hcp_mmp_subj'
    shell: 'FS_LICENSE={params.license} mris_convert {input} {output} &> {log}'

rule convert_to_nifti:
    input: join(config['in_freesurfer'],'mri','{volname}.mgz')
    output: 'results/hcp_mmp/sub-{subject}/{volname}.nii.gz'
    params: 
        license = config['fs_license']
    container: config['singularity_freesurfer']
    log: 'logs/convert_to_nifti/sub-{subject}_{volname}.log'
    group: 'hcp_mmp_subj'
    shell: 'FS_LICENSE={params.license} mri_convert {input} {output} &> {log}'

rule get_tkr2scanner:
    input: 
        t1 = 'results/hcp_mmp/sub-{subject}/T1.nii.gz'
    output:
        tkr2scanner = 'results/hcp_mmp/sub-{subject}/tkr2scanner.xfm'
    params: 
        license = config['fs_license']
    container: config['singularity_freesurfer']
    log: 'logs/get_tkr2scanner/sub-{subject}.log'
    group: 'hcp_mmp_subj'
    shell: 'FS_LICENSE={params.license} mri_info {input.t1} --tkr2scanner > {output.tkr2scanner} 2> {log}'
     
rule apply_surf_tkr2scanner:
    input: 
        surf = 'results/hcp_mmp/sub-{subject}/{hemi}_space-fsaverage.{surfname}.surf.gii',
        tkr2scanner = 'results/hcp_mmp/sub-{subject}/tkr2scanner.xfm'
    output: 
        surf = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native.{surfname}.surf.gii'
    threads: 8
    container: config['singularity_connectome_workbench']
    log: 'logs/apply_surf_tkr2scanner/sub-{subject}_{hemi}_{surfname}.log'
    group: 'hcp_mmp_subj'
    shell: 'wb_command -surface-apply-affine {input.surf} {input.tkr2scanner} {output.surf} &> {log}'


rule gen_midthickness:
    input:
        white_fsaverage = 'results/hcp_mmp/sub-{subject}/{hemi}_space-fsaverage.white.surf.gii',
        pial_fsaverage = 'results/hcp_mmp/sub-{subject}/{hemi}_space-fsaverage.pial.surf.gii',
	white_native = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native.white.surf.gii',
	pial_native = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native.pial.surf.gii'
    output: 
        midthickness_fsaverage = 'results/hcp_mmp/sub-{subject}/{hemi}_space-fsaverage.midthickness.surf.gii',
	midthickness_native = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native.midthickness.surf.gii'
    container: config['singularity_connectome_workbench']
    threads: 8
    log: 'logs/gen_midthickness/sub-{subject}_{hemi}.log'
    group: 'hcp_mmp_subj'
    shell: 'wb_command -surface-average {output.midthickness_fsaverage} -surf {input.white_fsaverage} -surf {input.pial_fsaverage} && '
           'wb_command -surface-average {output.midthickness_native} -surf {input.white_native} -surf {input.pial_native} &> {log}'
   

rule resample_subj_to_fsaverage_sphere:
    input: 
        surf = 'results/hcp_mmp/sub-{subject}/{hemi}_space-fsaverage.midthickness.surf.gii',
        current_sphere = 'results/hcp_mmp/sub-{subject}/{hemi}_space-fsaverage.sphere.reg.surf.gii',
        new_sphere = lambda wildcards: 'resources/standard_mesh_atlases/resample_fsaverage/'
                                        'fs_LR-deformed_to-fsaverage.{H}.sphere.32k_fs_LR.surf.gii'.format(
                                                H = hemi_to_H[wildcards.hemi])
    params:
        method = 'BARYCENTRIC'
    output:
        surf = 'results/hcp_mmp/sub-{subject}/{hemi}_space-32k_fs_LR.midthickness.surf.gii'
    container: config['singularity_connectome_workbench']
    threads: 8
    log: 'logs/resample_subj_to_fsaverage_sphere/sub-{subject}_{hemi}.log'
    group: 'hcp_mmp_subj'
    shell: 'wb_command -surface-resample {input.surf} {input.current_sphere} {input.new_sphere} {params.method} {output.surf} &> {log}'


rule resample_labels_to_subj_sphere:
    input:
        label = 'resources/standard_mesh_atlases/{hemi}.hcp-mmp.32k_fs_LR.label.gii',
        current_sphere = lambda wildcards: 'resources/standard_mesh_atlases/resample_fsaverage/'
                                            'fs_LR-deformed_to-fsaverage.{H}.sphere.32k_fs_LR.surf.gii'.format(
                                                H = hemi_to_H[wildcards.hemi]),
        new_sphere = 'results/hcp_mmp/sub-{subject}/{hemi}_space-fsaverage.sphere.reg.surf.gii',
        current_surf = 'results/hcp_mmp/sub-{subject}/{hemi}_space-32k_fs_LR.midthickness.surf.gii',
        new_surf = 'results/hcp_mmp/sub-{subject}/{hemi}_space-fsaverage.midthickness.surf.gii'
    params:
        method = 'ADAP_BARY_AREA'
    output:
        label = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native_label-hcpmmp.label.gii'
    container: config['singularity_connectome_workbench']
    threads: 8
    log: 'logs/resample_labels_to_subj_sphere/sub-{subject}_{hemi}.log'
    group: 'hcp_mmp_subj'
    shell: 
        'wb_command -label-resample {input.label} {input.current_sphere} {input.new_sphere}'
        ' {params.method} {output.label}'
        ' -area-surfs {input.current_surf} {input.new_surf} &> {log}'



rule map_labels_to_volume_ribbon:
    input: 
        label = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native_label-hcpmmp.label.gii',
        surf = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native.midthickness.surf.gii',
        vol_ref = 'results/hcp_mmp/sub-{subject}/T1.nii.gz',
        white_surf = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native.white.surf.gii',
        pial_surf = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native.pial.surf.gii',
    output:
        label_vol = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native_label-hcpmmp.nii.gz',
    container: config['singularity_connectome_workbench']
    threads: 8
    log: 'logs/map_labels_to_volume_ribbon/sub-{subject}_{hemi}.log'
    group: 'hcp_mmp_subj'
    shell:
        'wb_command -label-to-volume-mapping {input.label} {input.surf} {input.vol_ref} {output.label_vol}'
        ' -ribbon-constrained {input.white_surf} {input.pial_surf}'
        ' -greedy &> {log}'



#We now have label maps for each hemisphere that range from values 1=V1 --> 180=p24, however label 110 pertains to the piriform cortex poorly segmented on the Glasser hcp_mmp atlas, and is thus removed from each hemisphere to yield 179 regions per hemisphere
rule remove_piriform_label_110:
    input:
        label_vol = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native_label-hcpmmp.nii.gz'
    output:
        label_vol_removepir110 = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native_label-hcpmmp_removepir110.nii.gz',
    container: config['singularity_neuroglia']
    threads: 8
    log: 'logs/remove_piriform_label_110/sub-{subject}_{hemi}.log'
    group: 'hcp_mmp_subj'
    shell:
        'c3d {input.label_vol} -replace 110 0 -o {output.label_vol_removepir110} &> {log}'


#Since the piriform label 110 was removed, a value of 1 was subtracted from labels 111-->180 for continuity in each hemisphere. This yields labels in each hemisphere from 1-->179, with what used to be labels 111-->180 becoming relabelled as 110-->179.
### A) Done by first thresholding the piriform removed volume below 111 such that 1-110=0 and 111-180 retain original values, then 1 is subtracted from this volume such that values 111(AVI)-180(p24) become 110(AVI)-179(p24), note that the -1 values will be removed later when taking maximum of this image with following image. 
### B) Additionally, the original piriform removed volume was also thresholded above 109 such that 110-180=0 and 1(V1)-109(MI) retain original values. 
#These maximum of these two images A) 1-109(MI)=-1 and 110(AVI)-179=110=179 B) 1-109=1-109 and 110-180=0, then yields the relabelled piriform vol for each hemisphere from 1-179 without skipping any values due to removal of the piriform
rule relabel_atlas_without_piriform:
    input:
        label_vol_removepir110 = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native_label-hcpmmp_removepir110.nii.gz'
    output:
        tempA_thrBelow111_sub1 = 'results/hcp_mmp/sub-{subject}/tempA_{hemi}_space-native_label-hcpmmp_removepir110_thrBelow111_subj1_vals110-179remain.nii.gz',
        tempB_thrAbove109 = 'results/hcp_mmp/sub-{subject}/tempB_{hemi}_space-native_label-hcpmmp_removepir110_thrAbove109_vals1-109remain.nii.gz',
        relabelled_vol_removepir110 = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native_label-hcpmmp_removepir110_relabelled_1-179.nii.gz'
    container: config['singularity_neuroglia']
    threads: 8
    log: 'logs/relabel_atlas_without_piriform/sub-{subject}_{hemi}.log'
    group: 'hcp_mmp_subj'
    shell:
        'fslmaths {input.label_vol_removepir110} -thr 111 -sub 1 {output.tempA_thrBelow111_sub1} && '
        'fslmaths {input.label_vol_removepir110} -uthr 109 {output.tempB_thrAbove109} && '
        'fslmaths {output.tempA_thrBelow111_sub1} -max {output.tempB_thrAbove109} {output.relabelled_vol_removepir110} &> {log}'

#Since we want to analyze the connectivity of each hemisphere's piriform to each hemisphere independently, we add 179 to the rh_relabelledmap to yield 378 final regions, with regions lh = 1-->179 and rh = 180-->358. Note that since the background of rh has a value of zero, its value becomes 179 after addition and thus the image is thresholded below 180 to revert to the original background label and not overlap with the lh values
rule relabel_rh_labels_to_separate_hemis_without_piriform:
    input:
        relabelled_rh_vol_removepir110 = 'results/hcp_mmp/sub-{subject}/rh_space-native_label-hcpmmp_removepir110_relabelled_1-179.nii.gz' 
    output:
        relabelled_rh_vol_removepir110_add179 = 'results/hcp_mmp/sub-{subject}/rh_space-native_label-hcpmmp_removepir290_relabelled_180-358.nii.gz'
    container: config['singularity_neuroglia']
    threads: 8
    log: 'logs/relabel_rh_labels_to_separate_hemis_without_piriform/sub-{subject}_rh.log'
    group: 'hcp_mmp_subj'
    shell:
        'fslmaths {input.relabelled_rh_vol_removepir110} -add 179 -thr 180 {output.relabelled_rh_vol_removepir110_add179} &> {log}'



#currently optional
rule map_labels_to_volume_wmboundary:
    input: 
        label = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native_label-hcpmmp.label.gii',
        surf = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native.white.surf.gii',
        vol_ref = 'results/hcp_mmp/sub-{subject}/T1.nii.gz',
        white_surf = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native.white.surf.gii',
        pial_surf = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native.pial.surf.gii',
    params:
        nearest_vertex = '{wmbdy}' 
    output:
        label_vol = 'results/hcp_mmp/sub-{subject}/{hemi}_space-native_label-hcpmmp_wmbound{wmbdy}.nii.gz',
    container: config['singularity_connectome_workbench']
    threads: 8
    log: 'logs/map_labels_to_volume_wmboundary/sub-{subject}_{hemi}_wmbound-{wmbdy}.log'
    group: 'hcp_mmp_subj'
    shell:
        'wb_command -label-to-volume-mapping {input.label} {input.surf} {input.vol_ref} {output.label_vol}'
        ' -nearest-vertex {params.nearest_vertex} &> {log}'
 


