import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy import stats

targetnames = pd.read_table('results/diffparc/sub-100610/target_images.txt', names=['TargetNames'], header=None)
targetnamesList = targetnames.TargetNames.to_list()
targetnamesList = [ line.replace('results/diffparc/sub-100610/targets/','') for line in targetnamesList ] #Remove that string from each line
targetnamesList = [ line.replace('_ROI.nii.gz','') for line in targetnamesList ] #Remove this string from each
cluster_range = range(2,snakemake.params.max_k+1)
for i,this_k in enumerate(cluster_range): #For 0,2 1,3 2,4
    label = 1
    thiskdirslist = list()
    thiskString = str(this_k)
    for j,thiskdir in enumerate(snakemake.input.connmap_Clusters_npz_dir):
        if thiskdir.endswith(thiskString):
            thiskdirslist.append(thiskdir)
    newClusterGroupDir = '{snakeout}_k-{thisk}'.format(snakeout=snakemake.params.connmap_Clusters_group_npz_dir_base,thisk=this_k)
    if not os.path.exists(newClusterGroupDir):
        os.mkdir(newClusterGroupDir)
    fig = plt.figure(figsize=[6.4*this_k,4.8])
    while label <= this_k:
        thislabel_namesList = list()
        labelString = str(label) + '_connmap.npz'
        print('labelString = ',labelString)
        for l,eachthiskdir in enumerate(thiskdirslist):
            thislabel_names = [f for f in os.listdir(path=eachthiskdir) if f.endswith(labelString)]
            thislabel_namesString = ''.join(thislabel_names)
            thislabel_fullpath = eachthiskdir + '/' + thislabel_namesString
            thislabel_namesList.append(thislabel_fullpath)
        thisk_thislabelonesubj = thislabel_namesList[0]
        print('thisk_thislabelonesubj = ',thisk_thislabelonesubj)
        data_thisk_thislabelonesubj = np.load(thisk_thislabelonesubj)
        affine = data_thisk_thislabelonesubj['affine']
        mask = data_thisk_thislabelonesubj['mask']
        conn_shape = data_thisk_thislabelonesubj['conn'].shape 
        nsubjects = len(thislabel_namesList)
        print('nsubjects = ',nsubjects)
        conn_group = np.zeros([nsubjects,conn_shape[0],conn_shape[1]])
        savestring = '{newClusterGroupD}/cluster{k_i:02d}_connmap_group.npz'.format(newClusterGroupD=newClusterGroupDir,k_i=label)
        for sub,npz in enumerate(thislabel_namesList):
            data_thisk_thislabel_eachsub = np.load(npz)
            conn_group[sub,:,:] = data_thisk_thislabel_eachsub['conn']
        np.savez(savestring, conn_group=conn_group,mask=mask,affine=affine)
        connmatrix_for_each_target_avggroup = conn_group.mean(axis=0) #Average across rows (subjects) to get connmap matrix average across all subjects
        print('connmatrix_for_each_target_avggroup shape = ',connmatrix_for_each_target_avggroup.shape)
        avg_conn_for_each_target_avggroup = connmatrix_for_each_target_avggroup.mean(axis=0) #Average across rows (voxels in piriform, indexed) to get connmap matrix average across all subjects
        sum_avg_conn_for_each_target_avggroup = np.sum(avg_conn_for_each_target_avggroup)
        quotients = [number / sum_avg_conn_for_each_target_avggroup for number in avg_conn_for_each_target_avggroup]
        #Create dataframe that includes columns target names and the values of the targets
        df = pd.DataFrame({'targetnames': targetnamesList, 'ConnectivityScore': avg_conn_for_each_target_avggroup, 'ConnectivityScoreNormalized': quotients}, index=range(0,179))
        dfsorted = df.sort_values(by='ConnectivityScoreNormalized', ascending=False) #Sort the target regions descending
        print("Targets sorted based on avg conn for k-%s and label %s = %s" % (this_k, label, dfsorted))
        #Plot the 10 targets with the highest connectivity to the piriform
        x = np.linspace(1,10,10) #start,stop,number of ticks
        y = np.linspace(0,1,11) #Makes it incerement by 0.1
        ax = fig.add_subplot(1,this_k,label)
        ax.set_title("k-%s Cluster-0%s" % (this_k, label))
        ax.set_xticks(x)
        ax.set_yticks(y)
        ax.set_ylim([0,0.5])
        ax.set_xticklabels(dfsorted.targetnames[0:10], rotation='vertical')
        ax.set_xlabel('Target Name')
        ax.set_ylabel('Average Connectivity to Target')
        ax.bar(x,dfsorted.ConnectivityScoreNormalized[0:10]) #Gives the target regions with the top 10 highest average connectivites from the group averaged piriform seed
        label +=1
    figuresavestring = '{newClusterGroupDir}/k-{thisk}_Individual_Clusters_Top10connectivity_to_targets.png'.format(newClusterGroupDir=newClusterGroupDir,thisk=this_k)
    print('figuresavestring is : ',figuresavestring)
    plt.savefig(figuresavestring)
