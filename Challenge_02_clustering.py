#!/usr/bin/env python
import argparse
import logging as lg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
from sklearn import decomposition
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
import sys



if __name__ == '__main__':
    lg.basicConfig(filename='Challenge_02.log', filemode='w', format='%(name)s - %(levelname)s - %(message)s')
    logger = lg.getLogger('Challenge_02')
    logger.setLevel(lg.INFO)
    
    stdout_handler = lg.StreamHandler(sys.stdout)
    stdout_handler.setLevel(lg.INFO)
    logger.addHandler(stdout_handler)
    
    parser = argparse.ArgumentParser(prog = 'Challenge_02', 
                                     description = 'Program to perform the second task in the Amphora Code Challenge')
    parser.add_argument('-g', '--geno_file')
    parser.add_argument('-c', '--coord_file')
    parser.add_argument('-o', '--out_file')
    args = parser.parse_args()
    
    if args.geno_file == None:
        logger.error(f'ERROR: genotype table is required!!!')
    elif args.coord_file == None:
        logger.error(f'ERROR: coordination file is required!!!')
    else:
        logger.info(f'Reading genotype table from file: {args.geno_file}...')
        data = pd.read_csv(args.geno_file, low_memory=False)
        data.index = data.CHROM.astype(str) + ':' + data.POS.astype(str)
        data = data.dropna(axis=0).drop(['CHROM', 'POS', 'REF', 'ALT'], axis=1).astype('Int64').transpose()
        logger.info(f'Processing {data.shape[0]} individuals and {data.shape[1]} variants')

        logger.info(f'Reading coordination table from file {args.coord_file}...')
        data = data.join(pd.read_csv(args.coord_file, sep='\t', lineterminator='\n').set_index('UUID'))
        data.loc[data['Superpopulation code'].isnull(), 'Superpopulation code'] = 'Unknown'
        coord_df = pd.DataFrame(data={'Superpopulation code':data['Superpopulation code']})
        data = data.drop('Superpopulation code', axis=1)
        coord_classes, numeric_classes = np.unique(coord_df['Superpopulation code'], return_inverse=True)
        coord_df['numeric_classes'] = numeric_classes

        ##################################################
        ### compute PCA
        logger.info(f'computing PCA...')
        pca = decomposition.PCA(n_components=10)
        pca.fit(data)
        pca_val = pca.transform(data)
        
        ##################################################
        ### compute hierarchical clustering
        logger.info(f'computing hierarchical clustering...')
        hierarchical_cluster = AgglomerativeClustering(n_clusters=5, metric='euclidean', linkage='ward')
        coord_df['hierarchical'] = hierarchical_cluster.fit_predict(data)
        
        ##################################################
        ### compute k means
        logger.info(f'computing k means clustering...')
        km = KMeans(n_clusters=5, init='random', n_init=10, max_iter=300, tol=1e-04, random_state=0)
        coord_df['kmeans'] = km.fit_predict(data)
        
        #############################################
        ### saving
        if args.out_file == None:
            logger.warning(f'WARNING: Output file not specified. Results were not saved')
        else:
            coord_df.to_csv(args.out_file)
            logger.info(f'clustering results saved in file {args.out_file}')
            
        
        #############################################
        ### ploting
        logger.info(f'Ploting...')
        cols = ['gray',  '#440154FF', '#3B528BFF', '#21908CFF', '#5DC863FF', '#FDE725FF']

        fig = plt.figure(1, figsize=(12, 9))
        ### pca plot
        sub_pca = fig.add_subplot(1, 3, 1)
        sub_pca.title.set_text('PCA vals')
        sub_pca.scatter(pca_val[coord_df['Superpopulation code'] == 'Unknown', 0], 
                        pca_val[coord_df['Superpopulation code'] == 'Unknown', 1], color=cols[0], edgecolor=None, alpha=0.5, s=15)
        sub_pca.scatter(pca_val[coord_df['Superpopulation code'] == 'AFR', 0], 
                        pca_val[coord_df['Superpopulation code'] == 'AFR', 1], color=cols[1], edgecolor=None, alpha=0.5, s=15)
        sub_pca.scatter(pca_val[coord_df['Superpopulation code'] == 'AMR', 0], 
                        pca_val[coord_df['Superpopulation code'] == 'AMR', 1], color=cols[2], edgecolor=None, alpha=0.5, s=15)
        sub_pca.scatter(pca_val[coord_df['Superpopulation code'] == 'EAS', 0], 
                        pca_val[coord_df['Superpopulation code'] == 'EAS', 1], color=cols[3], edgecolor=None, alpha=0.5, s=15)
        sub_pca.scatter(pca_val[coord_df['Superpopulation code'] == 'EUR', 0], 
                        pca_val[coord_df['Superpopulation code'] == 'EUR', 1], color=cols[4], edgecolor=None, alpha=0.5, s=15)
        sub_pca.scatter(pca_val[coord_df['Superpopulation code'] == 'SAS', 0], 
                        pca_val[coord_df['Superpopulation code'] == 'SAS', 1], color=cols[5], edgecolor=None, alpha=0.5, s=15)
        sub_pca.legend(coord_df['Superpopulation code'])

        ### hierarchical plot
        sub_hier = fig.add_subplot(1, 3, 2)
        sub_hier.title.set_text('Hierarchical clustering')
        sub_hier.scatter(pca_val[coord_df['hierarchical'] == 0, 0],
                        pca_val[coord_df['hierarchical'] == 0, 1], color=cols[2], edgecolor=None, alpha=0.5, s=15)
        sub_hier.scatter(pca_val[coord_df['hierarchical'] == 1, 0],
                        pca_val[coord_df['hierarchical'] == 1, 1], color=cols[1], edgecolor=None, alpha=0.5, s=15)
        sub_hier.scatter(pca_val[coord_df['hierarchical'] == 2, 0],
                        pca_val[coord_df['hierarchical'] == 2, 1], color=cols[3], edgecolor=None, alpha=0.5, s=15)
        sub_hier.scatter(pca_val[coord_df['hierarchical'] == 3, 0],
                        pca_val[coord_df['hierarchical'] == 3, 1], color=cols[5], edgecolor=None, alpha=0.5, s=15)
        sub_hier.scatter(pca_val[coord_df['hierarchical'] == 4, 0],
                        pca_val[coord_df['hierarchical'] == 4, 1], color=cols[4], edgecolor=None, alpha=0.5, s=15)
        # sub_hier.legend(coord_df['Superpopulation code'])

        ### kmeans plot
        sub_kms = fig.add_subplot(1, 3, 3)
        sub_kms.title.set_text('KMeans clustering')
        sub_kms.scatter(pca_val[coord_df['kmeans'] == 0, 0],
                        pca_val[coord_df['kmeans'] == 0, 1], color=cols[3], edgecolor=None, alpha=0.5, s=15)
        sub_kms.scatter(pca_val[coord_df['kmeans'] == 1, 0],
                        pca_val[coord_df['kmeans'] == 1, 1], color=cols[4], edgecolor=None, alpha=0.5, s=15)
        sub_kms.scatter(pca_val[coord_df['kmeans'] == 2, 0],
                        pca_val[coord_df['kmeans'] == 2, 1], color=cols[2], edgecolor=None, alpha=0.5, s=15)
        sub_kms.scatter(pca_val[coord_df['kmeans'] == 3, 0],
                        pca_val[coord_df['kmeans'] == 3, 1], color=cols[1], edgecolor=None, alpha=0.5, s=15)
        sub_kms.scatter(pca_val[coord_df['kmeans'] == 4, 0],
                        pca_val[coord_df['kmeans'] == 4, 1], color=cols[5], edgecolor=None, alpha=0.5, s=15)
        # sub_kms.legend(coord_df['Superpopulation code'])
                        
        plt.show()
