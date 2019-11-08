#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 11:24:49 2018

@author: rhou
"""

import warnings
warnings.filterwarnings("ignore")
import argparse
import pandas as pd
import os, sys, glob
from scipy import io

def CalCords(savefolder, em, visMethod):
    cordFile = os.path.join(savefolder, 'cords_%s.csv' % visMethod)
    if os.path.exists(cordFile):
        cords = pd.read_csv(cordFile, index_col=0, header=0)
        return cords
        
    x = em.T.values
    from sklearn.preprocessing import StandardScaler
    x = StandardScaler().fit_transform(x)
    from sklearn.decomposition import PCA
    pca = PCA(n_components=10)
    principalComponents = pca.fit_transform(x)
    
    if visMethod == 'PCA':
        cords = pd.DataFrame(data = principalComponents[:,:2], columns = ['x', 'y'], index = em.columns)
    elif visMethod == 'tSNE':
        from sklearn.manifold import TSNE
        tsneComponents = TSNE(n_components=2).fit_transform(principalComponents)
        cords = pd.DataFrame(data = tsneComponents, columns = ['x', 'y'], index = em.columns)
    elif visMethod == 'UMAP':
        import umap
        umapComponents = umap.UMAP(n_components=2).fit_transform(principalComponents)
        cords = pd.DataFrame(data = umapComponents, columns = ['x', 'y'], index = em.columns)
    cords.index.name = 'barcode'    
    cords.to_csv(cordFile, index=True, header=True)
    return cords

def DrawScatters(savefolder, annoFile, visMethod, cords, annos):
    import plotly
    import plotly.graph_objs as go
    annText = os.path.basename(annoFile).split('.')[0]
        
    for kind in ['cell type', 'top sample']:
        if kind not in annos.columns:
            continue
        annotationList = sorted(list(set(annos.ix[:,kind])))
        
        import seaborn as sns 
        colorList = sns.hls_palette(n_colors=len(annotationList))
        
        data = []
        annoLen = 0
        for annoIdx in range(len(annotationList)):
            annoNames = annotationList[annoIdx]
            if len(annoNames) > annoLen:
                annoLen = len(annoNames)
            indicesOfAnno = annos[kind]==annoNames
            text = []
            for idx in annos.index[indicesOfAnno]:
                show_text = '%s: %s, barcode: %s' % (kind, annoNames, idx)
                text.append(show_text)
            trace = go.Scatter(
                x = cords.ix[annos.index[indicesOfAnno],'x'],
                y = cords.ix[annos.index[indicesOfAnno],'y'],
                name = annoNames,
                mode = 'markers',
                marker=dict(
                    color='rgb(%s, %s, %s)' % colorList[annoIdx],
                    size=5,
                    symbol='circle',
                    line=dict(
                        color='rgb(204, 204, 204)',
                        width=1
                    ),
                    opacity=0.9
                ),
                text = text,
             )
            data.append(trace)
        if annoLen < 35:
            layout = go.Layout(legend=dict(orientation="v"),autosize=True,showlegend=True)
        else:
            layout = go.Layout(legend=dict(orientation="v"),autosize=True,showlegend=False)
        fig = go.Figure(data=data, layout=layout)
        fn = os.path.join(savefolder, '%s_%s_%s.html' % (annText, kind.replace(' ', '_'), visMethod))
        print('##########saving plot: %s' % fn)
        plotly.offline.plot(fig, filename=fn)

#start to visualise test dataset
def main(testFormat, testDS, annoFile, visMethod):
    #load test data
    print('##########loading test data')
    if testFormat == '10x':
        fileItem = glob.glob(os.path.join(testDS, "matrix.mtx"))[0]
        em = io.mmread(fileItem)
        em = em.tocsr().toarray()
        row = pd.read_table(fileItem[:-10]+"features.tsv", header=None, index_col=None)
        col = pd.read_table(fileItem[:-10]+"barcodes.tsv", header=None, index_col=None)
        em = pd.DataFrame(em, index=row.T.values[1], columns=col.T.values[0])
        savefolder = testDS
    else:
        em = pd.read_csv(testDS, index_col=0, header=0)
        savefolder = testDS[:-4]
        
    print('##########reducing dimensions')
    cords = CalCords(savefolder, em, visMethod)
    annos = pd.read_csv(annoFile, index_col=0, header=0)
    commonIdx = set(cords.index).intersection(set(annos.index))
    cords = cords.ix[commonIdx,]
    annos = annos.ix[commonIdx,]
    
    print('##########darwing the scatter plots in the folder: %s' % savefolder)
    DrawScatters(savefolder, annoFile, visMethod, cords, annos)
    print('##########DONE!')

if __name__ == "__main__":
    
    #process arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--testDS', required=True, help='path to the folder of test dataset if dtype is 10x, otherwise, the path to the file')
    parser.add_argument('--dFormat', default='10x', help='10x (default) | csv')
    parser.add_argument('--annoFile', required=True, help='path to the annotation file generated by scMatch')
    parser.add_argument('--visMethod', default='p', help='p[ca] (default) | t[sne] | u[map]')
    
    opt = parser.parse_args()
    
    #check visMethod
    if opt.visMethod.lower() == "p":
        visMethod = 'PCA'
    elif opt.visMethod.lower() == "t":
        visMethod = 'tSNE'
    elif opt.visMethod.lower() == "u":
        visMethod = 'UMAP'
    else:
        sys.exit("The visMethod can only be 'p', 't' or 'u'.")
    
    #check dFormat and testDS
    if opt.dFormat.lower() == 'csv':
        if not os.path.exists(opt.testDS):
            sys.exit("The test dataset file does not exist.")
        if opt.testDS[-4:].lower() != '.csv':
            sys.exit("The test dataset file is not a CSV file.")
    elif opt.dFormat.lower() == '10x':
        if not os.path.exists(opt.testDS):
            sys.exit("The folder of test dataset does not exist.")
        if not os.path.exists(os.path.join(opt.testDS, 'matrix.mtx')):
            sys.exit("Cannot find 'matrix.mtx' file in the folder of test dataset.")
        if not os.path.exists(os.path.join(opt.testDS, 'features.tsv')):
            sys.exit("Cannot find 'features.tsv' file in the folder of test dataset.")
        if not os.path.exists(os.path.join(opt.testDS, 'barcodes.tsv')):
            sys.exit("Cannot find 'barcodes.tsv' file in the folder of test dataset.")
            
    #check annotation file 
    if not os.path.exists(opt.annoFile):
        sys.exit("Cannot find annotation file for single cells.")
        
    #pass argument check, show input data
    print('===================================================')
    print('Input data:')
    print('The format of the test dataset: %s' % opt.dFormat)
    if opt.dFormat == '10x':
        print('Test data are in the folder: %s' % opt.testDS)
    else:
        print('Test data are in the file: %s' % opt.testDS)
    print('Annotation file: %s' % opt.annoFile)
    print('The visualisation method: %s' % visMethod)
    print('===================================================')
    
    #start to visualise test dataset
    main(opt.dFormat, opt.testDS, opt.annoFile, visMethod)
    #python visAnno.py --dFormat csv --testDS GSE81861_Cell_Line_COUNT.csv --annoFile GSE81861_Cell_Line_COUNT/annotation_result_keep_all_genes/human_Spearman_top_ann.csv
