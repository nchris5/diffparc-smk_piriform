import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import nibabel as nib
from functools import reduce
import os


#Set the name of the seed you want to create the Gephi input for
seeds = snakemake.params.seeds_py

#Set the k parcellation you want to create the Gephi input for
#this_k = 3
cluster_range = range(2,snakemake.params.max_k+1)

#Set the path to the location of the directory containing the connmapGroup directories for each seed and each k 
#(i.e. connmapClusters dir from diffparc located in /scratch/nchris5/HCP_template_piriform_segmentations/diffparc/results)
path = snakemake.params.connmapClusters_dir


#Read in the Glasser_2016_Table with the 'SectionsBOLD' column manually included, which indicates which of the 22 overarching regions each of the 180 regions belongs to; 
#and the 'ColourRGB' column manually included, which indicates the RGB for that specific SectionsBOLD
GlasserTable = pd.read_csv(snakemake.params.GlasserTable_22regions_colours, header=0, usecols=range(1,8))
GlasserTable.index = range(1,180)
print('Original GlasserTable : \n', GlasserTable)
#Sort by the section values from 1-->22
GlasserTable_sorted_PrimarySection = GlasserTable.sort_values(by='SectionsBOLD', ascending=True)
#Change the actual index to 1-->179
GlasserTable_sorted_PrimarySection.index = range(1,180)
#Insert a column that describes the sorted dataframe new index from 1-->180
GlasserTable_sorted_PrimarySection.insert(loc=1, column='Sorted_PrimarySection_Index', value=range(1,180))
#Create the list of target ID's with leading zeros so all ID's have 3 digits (i.e. 001, 010, 100), as for gephi to read in proper order
GlasserTable_sorted_PrimarySection_Index_ThreeDigits = [str(indID).zfill(3) for indID in GlasserTable_sorted_PrimarySection.Sorted_PrimarySection_Index]
print('GlasserTable_sorted_PrimarySection : \n', GlasserTable_sorted_PrimarySection)

print('\n ------ Creating Nodes df ------ \n ')
#Includes the Id, Label, color, and SectionsBOLD for each of 180 Glasser parcellations (1-180)
Targets = pd.DataFrame({'Id': GlasserTable_sorted_PrimarySection_Index_ThreeDigits, 'Label': GlasserTable_sorted_PrimarySection.AreaName, 'color': GlasserTable_sorted_PrimarySection.ColourRGB, 'Modularity Class': GlasserTable_sorted_PrimarySection.SectionsBOLD}, index=range(1,180))
print('Targetsdf : \n', Targets)

for seed in seeds:
    for i,this_k in enumerate(cluster_range): #For 0,2 1,3 2,4
        #Set the RBG colours for each of the clusters within this_k to the same colour that is present in the clustering output in itksnap (i.e. diffparc spectral_clustering output)
        if this_k==2:
            this_k_ColourRGB = ['255,0,0','0,255,0'] #Cluster01=Red, #Cluster02=Green
        elif this_k==3:
            this_k_ColourRGB = ['255,0,0','0,255,0','0,0,255'] #Cluster01=Red, #Cluster02=Green, #Cluster03=Blue
        elif this_k==4:
            this_k_ColourRGB = ['255,0,0','0,255,0','0,0,255', '255,255,0'] #Cluster01=Red, #Cluster02=Green, #Cluster03=Blue #Yellow
        #Create the Id's and labels for the individual clusters  for and get the basename for the individual clusters (e.g. Cluster01, Cluster02, Cluster03), assigning Id values above 180 (since 180 targets) e.g. 181=Cluster01, 182=Cluster02, 183=Cluster03
        IdlistClusters = list()
        labellist = list()
        for i,this_label in enumerate(range(1,this_k+1)):
            IdlistClusters.append(180+i)
            labellist.append('Cluster{this_label:02d}'.format(this_label=this_label))
        #Set the Modularity Class for each individual cluster to be distinct (i.e. Cluster01=2, Cluster02=3, Cluster03=4
        ModularityClassClusters = list(range(23,this_k+23))
        #Create the dataframe for the clusters to add to the list of targets
        individualclusterstoadd = pd.DataFrame({'Id': IdlistClusters, 'Label': labellist, 'color': this_k_ColourRGB, 'Modularity Class': ModularityClassClusters}, index=range(179,179+this_k))
        print('\n individualclusterstoadd : \n', individualclusterstoadd)

        #Append the cluster labels dataframe to the bottom of the targets list dataframe to yield the final Nodesdf
        Nodesdf = Targets.append(individualclusterstoadd)
        print('\n Nodesdf : \n ', Nodesdf)




        print('\n ------ Creating Edges df ------ \n ')
        print(print('------ Creating Edges df SourceId and SourceName columns ------ \n '))

        #Multiply the targets df by number of clusters (i.e. first 180 will represent cluster 1, next 180 will represent cluster 2, etc.)
        EdgesSourceId = np.zeros((179,this_k), dtype=int)
        EdgesSourceNames = list()
        for i,this_label in enumerate(range(1,this_k+1)):
            EdgesSourceId[:,i] = this_label+179
            tempSourceNames = list()
            for j in range(179):
                tempSourceNames.append('Cluster{this_label:02d}'.format(this_label=this_label))
            EdgesSourceNames.append(tempSourceNames)
        #Use functools reduce to remove iterative lists into a single list i.e. [179xCluster01], [179xCluster02], [179xCluster03] to [179xCluster01, 179xCluster02, 179xCluster03]
        EdgesSourceNames = reduce(lambda x,y: x+y, EdgesSourceNames)
        #Reshape the edges source col current with shape (179,numClusters) to 179*numClusters with order 'F' (i.e. 179x'181' --> 179x'182' --> 179x'183')
        EdgesSourceId = np.reshape(EdgesSourceId, (179*this_label), order='F')

        TargetsMultiplied = pd.concat([Targets]*this_k, ignore_index=True)
        print('TargetsMultiplied Shape = ', TargetsMultiplied.shape)
        print('TargetsMultiplied = \n', TargetsMultiplied)

        print('\n ------ Creating Edges df Columns:   Source, SLabel, Target, TLabel, ConnectivityScore, ConnectivityScoreNormalized, Rank_in_Cluster, Weight ------ \n ')

        #Create list of paths to individual clusternpz files for this_k
        connmapClustersDirList = list()
        thiskString = str(this_k)
        #Only select the directories that have the specified seed and also end with k=this_k
        for connmapClustersDir in os.listdir(path):
            if ((connmapClustersDir.find(seed) != -1) and (connmapClustersDir.endswith(thiskString))):
                connmapClustersDirList.append(connmapClustersDir)
                newpath = path + '/' + connmapClustersDirList[0]
                individualClusters_thisk_npzList = list()
                for individualClusters_thisk_npz in os.listdir(newpath):
                    if individualClusters_thisk_npz.endswith('.npz'):
                        fullpath = newpath + '/' + individualClusters_thisk_npz
                        individualClusters_thisk_npzList.append(fullpath)
                print(individualClusters_thisk_npzList)

        #Initiate an empty dataframe that will be appended to include the values from each cluster
        df_ALLCLUSTERS_InsertRank = pd.DataFrame(columns = ['targetnames', 'ConnectivityScore', 'ConnectivityScoreNormalized', 'SortedRank'])

        #Create a list of ranks from 1(lowest connectivity) to 180(highest connectivity)
        for i,thisk_thisCluster_npz in enumerate(individualClusters_thisk_npzList):
            data_thisCluster = np.load(thisk_thisCluster_npz)
            conngroup_thisCluster = data_thisCluster['conn_group']

            connAvg_across_subs_thisCluster = conngroup_thisCluster.mean(axis=0) #Shape of input is (#sub,#voxels,#targets)
            targetAvg_across_voxels_thisCluster = connAvg_across_subs_thisCluster.mean(axis=0) #Shape of input is (#voxels,#targets)
            sum_targetAvg_across_voxels_thisCluster = np.sum(targetAvg_across_voxels_thisCluster) #Sum all the targets to get sum of conn scores
            quotients = [number / sum_targetAvg_across_voxels_thisCluster for number in targetAvg_across_voxels_thisCluster]

            df_thisCluster = pd.DataFrame({'targetnames': GlasserTable.AreaName, 'ConnectivityScore': targetAvg_across_voxels_thisCluster, 'ConnectivityScoreNormalized': quotients}, index=range(1,180))
            #Insert a column the SectionsBold column from the GlasserTable
            df_thisCluster.insert(loc=1, column='SectionsBOLD', value=GlasserTable.SectionsBOLD)
            print('\n df_thisCluster_InsertedSectionsBold k-%s and cluster-%s : \n %s' % (this_k, i+1, df_thisCluster))
            df_thisCluster_sorted_PrimarySection =  df_thisCluster.sort_values(by='SectionsBOLD', ascending=True)
            df_thisCluster_sorted_PrimarySection.index = range(1,180)
            #Insert a column that describes the sorted dataframe new index from 1-->180
            df_thisCluster_sorted_PrimarySection.insert(loc=1, column='Sorted_PrimarySection_Index', value=range(1,180))
            print('\n df_thisCluster_sorted_PrimarySection k-%s and cluster-%s : \n %s' % (this_k, i+1, df_thisCluster_sorted_PrimarySection))

            #Assign weight based on ConnectivityScore
            df_sorted_thisCluster = df_thisCluster_sorted_PrimarySection.sort_values(by='ConnectivityScoreNormalized', ascending=False)
            print('\n df_sorted_thisCluster k-%s and cluster-%s : \n %s' % (this_k, i+1, df_sorted_thisCluster))

            df_sorted_thisCluster_InsertRank = df_sorted_thisCluster
            df_sorted_thisCluster_InsertRank['SortedRank'] = pd.Series(reversed(range(1,180)),index=(df_sorted_thisCluster.index))

            df_thisCluster_sorted_PrimarySection_InsertRank = df_thisCluster_sorted_PrimarySection
            df_thisCluster_sorted_PrimarySection_InsertRank['SortedRank'] = pd.Series(reversed(range(1,180)),index=(df_sorted_thisCluster.index))
            print('\n df_sorted_thisCluster k-%s and cluster-%s with Rank : \n %s' % (this_k, i+1, df_sorted_thisCluster_InsertRank))
            print('\n df_thisCluster k-%s and cluster -%s with Rank : \n %s' % (this_k, i+1, df_thisCluster_sorted_PrimarySection_InsertRank))

            df_ALLCLUSTERS_InsertRank = df_ALLCLUSTERS_InsertRank.append(df_thisCluster_sorted_PrimarySection_InsertRank, ignore_index=True)
            print('CURRENT df_ALLCLUSTERS_InsertRank Shape = ', df_ALLCLUSTERS_InsertRank.shape)
        print('FINAL df_ALLCLUSTERS_InsertRank Shape = ', df_ALLCLUSTERS_InsertRank.shape)
        print('FINAL df_ALLCLUSTERS_InsertRank : \n', df_ALLCLUSTERS_InsertRank)

        #Add in the edge weight column
        #Multiply the SortedRank*ConnectivityScoreNormalized to get the Weight of each edge, round to integer, then add 1 to ensure that lowest edge weight value is > 0
        df_ALLCLUSTERS_InsertRank['Weight'] = df_ALLCLUSTERS_InsertRank['SortedRank'] * df_ALLCLUSTERS_InsertRank['ConnectivityScoreNormalized']
        df_ALLCLUSTERS_InsertRank['Weight'] = df_ALLCLUSTERS_InsertRank['Weight'].apply(lambda x: round(x, 0) + 1)

        Edgesdf = pd.DataFrame({'Source': EdgesSourceId, 'SLabel': EdgesSourceNames, 'Target': TargetsMultiplied.Id, 'TLabel': TargetsMultiplied.Label, 'ConnectivityScore': df_ALLCLUSTERS_InsertRank.ConnectivityScore, 'ConnectivityScoreNormalized': df_ALLCLUSTERS_InsertRank.ConnectivityScoreNormalized, 'Rank_in_Cluster': df_ALLCLUSTERS_InsertRank.SortedRank, 'Weight': df_ALLCLUSTERS_InsertRank.Weight}, index=range(0,this_k*179))
        print('\n Shape of this Edges is %s and Edgesdf = \n %s  ' % (Edgesdf.shape, Edgesdf))


        print('\n ------ Saving Nodesdf and Edgees df to csv ------ \n ')

        if not os.path.exists('results/gephi_input_NodesandEdges_tables_' + seed):
            os.mkdir('results/gephi_input_NodesandEdges_tables_' + seed)

        #Create savepaths
        savepathNodesdf = 'results/gephi_input_NodesandEdges_tables_' + seed + '/GephiNodes_' + seed + '_k-' + str(this_k) + '_Sorted_Overarching22Categories_WithRGBColours.csv'
        savepathEdgesdf = 'results/gephi_input_NodesandEdges_tables_' + seed + '/GephiEdges_' + seed + '_k-' + str(this_k) + '_Sorted_Overarching22Categories_NormalizedWeight.csv'
        
        print('savepathNodesdf = ', savepathNodesdf)
        print('savepathEdgesdf = ', savepathEdgesdf)

        #Write to Nodesdf and Edgesdf to CSV files
        Nodesdf.to_csv(savepathNodesdf, index = False, header=True)
        Edgesdf.to_csv(savepathEdgesdf, index = False, header=True)
    

