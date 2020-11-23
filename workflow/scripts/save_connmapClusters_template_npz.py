import nibabel as nib
import numpy as np
import os

#kclusters_nii_list_names = os.listdir(path=snakemake.input.clusters_template_space_dir)
#kclusters_nii_list = list()
#for filename in kclusters_nii_list_names:
#    fullpath = os.path.join(snakemake.input.clusters_template_space_dir, filename)
#    kclusters_nii_list.append(fullpath)

kclusters_nii_list = list()
for a,thiskclustervol in enumerate(snakemake.input.clusters_template_space):
    kclusters_nii_list.append(thiskclustervol)
print('kclusters_nii_list = ',kclusters_nii_list)

cluster_range = range(2,snakemake.params.max_k+1) #Range of the clusters [2, 3, 4]
ntargets = len(snakemake.params.connmap_3d) #Number of targets = 180

for i,this_k in enumerate(cluster_range): #For 0,2 1,3 2,4
    label = 1 #Lowest label value is 1
    newClusterDir = '{snakeout}_k-{thisk}'.format(snakeout=snakemake.params.connmap_Clusters_npz,thisk=this_k)
    if not os.path.exists(newClusterDir):
        os.mkdir(newClusterDir)
    while label <= this_k: #For each individual label in that cluster file, stop after label hits that current k value
        print("Label val for k-%s prior to creating indices = %s" % (this_k, label))
        clustmask_nib = nib.load(kclusters_nii_list[i]) #Load in cluster 2, 3, 4
        print("This k-%scluster = %s" % (this_k, kclusters_nii_list[i]))
        clustmask_vol = clustmask_nib.get_fdata()
        clustmask_indices = clustmask_vol==label #Indices within the clustmask_vol for label=label (e.g. indices in k-2 colume where label = 1)
        maskedCount = clustmask_vol[clustmask_indices] #Number of voxels in a list that contain this label
        nvoxels = len(maskedCount) #Count of the number of voxels
        print("nvoxels for k-%s and label %s = %s" % (this_k, label, nvoxels))
        maskedVoxelsIndividualCluster = clustmask_vol
        maskedVoxelsIndividualCluster[maskedVoxelsIndividualCluster!=label] = 0 #Any of the voxels that do not contain the label value become 0 (e.g. if this_k = 2, and label = 1, then label = 2 becomes 0)
        conn = np.zeros((nvoxels,ntargets)) #Create conn matrix of size nvoxels in this cluster,ntargets
        savestring = '{newClustDir}/cluster{k_i:02d}_connmap.npz'.format(newClustDir=newClusterDir,k_i=label)
        print("savestring for k-%s and label %s = %s" % (this_k, label, savestring))
        for j,conn_file in enumerate(snakemake.params.connmap_3d):
            targvol = nib.load(conn_file).get_fdata()
            maskedtargvol = targvol[clustmask_indices].T
            conn[:,j] = maskedtargvol
        print("Conn shape for k-%s and label %s =  %s \n" % (this_k, label, conn.shape))
        np.savez(savestring, conn=conn,mask=maskedVoxelsIndividualCluster,affine=clustmask_nib.affine)
        label +=1 #Label becomes 2 and loop while loop continues, then label becomes 3 and the loop exits
