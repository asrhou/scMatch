#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 19:59:33 2018

@author: rhou
"""
import warnings
warnings.filterwarnings("ignore")
import argparse
import pandas as pd
import numpy as np
import os, sys
from scipy.stats import spearmanr, pearsonr
import glob
from scipy import io
import multiprocessing
from functools import partial

def SPMAnno(refDB, keepZeros, testMethod, testCol):
    firstLayerHeader = [testCol.columns[0], testCol.columns[0]]
    
    #find common genes between test data and ref data
    testRows = set(testCol.index)
    refRows = set(refDB.index)
    if keepZeros:
        commonRows = list(refRows)
    else:
        commonRows = list(refRows.intersection(testRows))
    commonRowsLen = len(commonRows)
    thirdLayerHeader = [commonRowsLen, commonRowsLen]
    
    #only keep non-zero genes
    testCol = testCol.loc[commonRows, ].fillna(0.0)
    testrefDB = refDB.loc[commonRows, ].fillna(0.0)
    
    if testMethod == 'Spearman':
        spr_correlation = testrefDB.apply(lambda col: spearmanr(col, testCol)[0], axis=0)
        spr_correlation = spr_correlation.to_frame().fillna(0).round(10).reset_index()
        spr_correlation.columns = ['sample name', 'Spearman correlation coefficient']
        secondLayerHeader = ['sample name', 'Spearman correlation coefficient']
        spr_correlation = spr_correlation.sort_values(by=['Spearman correlation coefficient'], ascending=False)
        spr_correlation = spr_correlation.reset_index(drop=True)
        return (firstLayerHeader, thirdLayerHeader, secondLayerHeader, spr_correlation, testCol.columns[0], spr_correlation.iloc[0,0], spr_correlation.iloc[0,1])
    elif testMethod == 'Pearson':
        pes_correlation = testrefDB.apply(lambda col: pearsonr(col.to_frame(), testCol)[0], axis=0)
        pes_correlation = pes_correlation.fillna(0).round(10).T.ix[:,0].reset_index()
        pes_correlation.columns = ['sample name', 'Pearson correlation coefficient']
        secondLayerHeader = ['sample name', 'Pearson correlation coefficient']
        pes_correlation = pes_correlation.sort_values(by=['Pearson correlation coefficient'], ascending=False)
        pes_correlation = pes_correlation.reset_index(drop=True)
        return (firstLayerHeader, thirdLayerHeader, secondLayerHeader, pes_correlation, testCol.columns[0], pes_correlation.iloc[0,0], pes_correlation.iloc[0,1])

#transfer given species gene symbols to hids
def TransferToHids(refDS, species, geneList):
    hidCol = 0  # universal id for homology genes
    taxidCol = 1  # species id
    geneSymbolCol = 3  # for jordan's data, it only has gene symbol

    homoDF = pd.read_csv(os.path.join(refDS, 'homologene.data'), sep='\t', index_col=None, header=None)
    # reduce the list to genes of the given species
    speciesDF = homoDF.ix[homoDF[taxidCol] == int(species),].set_index(geneSymbolCol)

    geneDF = speciesDF.loc[speciesDF.index.isin(geneList)]
    notnaGenes = geneDF[geneDF[hidCol].notnull()]
    notnaGenes = notnaGenes.ix[~notnaGenes[hidCol].duplicated(keep='first'),]
    notnaGenesSymbols = list(notnaGenes.index)
    notnaGenesHIDs = list(notnaGenes.ix[:, hidCol])
    notnaGenesHIDs = [int(i) for i in notnaGenesHIDs]
    return notnaGenesSymbols, notnaGenesHIDs

def AnnSCData(testType, em, refDS, refType, refTypeName, keepZeros, testMethod, coreNum, savefolder):
    #if across species need to replace symbols to HIDs
    if testType == refType:
        #load reference database
        refDB = pd.read_csv(os.path.join(refDS, '%s_symbol.csv' % refType), index_col=0, header=0)
    else:
        #load reference database
        refDB = pd.read_csv(os.path.join(refDS, '%s_HID.csv' % refType), index_col=0, header=0)
            
        #for different species, transfer symbols to hids
        geneList = list(em.index)
        geneList, hidList = TransferToHids(refDS, testType, geneList)
        em = em.ix[~em.index.duplicated(keep='first'),]
        em = em.ix[geneList, ]
        hidList = [str(i) for i in hidList]
        em.index = hidList
        hidList = [str(i) for i in refDB.index]
        refDB.index = hidList
        
    #remove duplicate indices
    refDB = refDB.ix[~refDB.index.duplicated(keep='first'),]
    print('reference dataset shape: %s genes, %s samples' % refDB.shape)
    em = em.ix[~em.index.duplicated(keep='first'),]
    
    #split expression matrix to single-cell expression profiles
    eps = np.split(em, len(em.columns), axis=1)
    
    #annotate single-cell expression profiles in parallel
    p = multiprocessing.Pool(coreNum)
    func = partial(SPMAnno, refDB, keepZeros, testMethod)
    resultList = p.map(func, eps)
    p.close()
    p.join()
    
    #generate final sample annotation expression matrix
    firstLayerHeader = [item for i in resultList for item in i[0]]
    secondLayerHeader = [item for i in resultList for item in i[1]]
    thirdLayerHeader = [item for i in resultList for item in i[2]]
    merged = pd.concat([i[3] for i in resultList], axis=1)
    arrays = [firstLayerHeader, secondLayerHeader, thirdLayerHeader]
    tuples = list(zip(*arrays))
    headers = pd.MultiIndex.from_tuples(tuples, names=['identifier', 'tested genes', 'annotation'])
    merged.columns = headers
    
    #prepare folder to save annotation result
    if not os.path.exists(savefolder):
        os.mkdir(savefolder)
        
    rstFolder = 'annotation_result'
    if keepZeros:
        rstFolder = rstFolder + '_keep_all_genes'
    else:
        rstFolder = rstFolder + '_keep_expressed_genes'
    
    savefolder = os.path.join(savefolder, rstFolder)
    if not os.path.exists(savefolder):
        os.mkdir(savefolder)
        
    #save file
    print('##########saving annotation results in the folder: %s' % savefolder)
    if len(merged.columns) > 8190:
        colNum = len(merged.columns)
        parts = int(colNum / 8000)
        for partIdx in range(parts):
            subMerged = merged.iloc[:,8000*partIdx:8000*(partIdx+1)]
            subMerged.to_excel(os.path.join(savefolder, refTypeName+"_%s_Part%04d.xlsx" % (testMethod, partIdx+1)))
        subMerged = merged.iloc[:,8000*parts:]
        subMerged.to_excel(os.path.join(savefolder, refTypeName+"_%s_Part%04d.xlsx" % (testMethod, parts+1)))
    else:
        merged.to_excel(os.path.join(savefolder, refTypeName+"_%s.xlsx" % testMethod))
        
    #save top ann result
    topAnn = pd.DataFrame({'cell':[i[4] for i in resultList], 'cell type':[''] * len(resultList), 'top sample':[i[5] for i in resultList], 'top correlation score':[i[6] for i in resultList]})
    mapData = pd.read_csv(os.path.join(refDS, '%s_map.csv' % refType), index_col=0, header=0)
    for idx in topAnn.index:
        topAnn.ix[idx, 'cell type'] = mapData.ix[topAnn.ix[idx, 'top sample'], 'cell type']
    saveNameP = os.path.join(savefolder, refTypeName+"_%s_top_ann.csv" % (testMethod))
    topAnn.to_csv(saveNameP, index=False, columns = ['cell', 'cell type', 'top sample', 'top correlation score'])
    print('##########DONE!')
    
def SortAnno(testItem):
    oldcols = testItem.columns
    cell = oldcols.get_level_values(0)[0]
    testItem.columns = ["name", "coefficient"]
    testItem = testItem.sort_values(by=['coefficient'], ascending=False)
    topAnn = testItem.iloc[0,0]
    topCoff = testItem.iloc[0,1]
    testItem.columns = oldcols
    testItem = testItem.reset_index(drop=True)
    return (testItem, cell, topAnn, topCoff)

#start to annotate test dataset
def main(testType, testFormat, testDS, testGenes, refDS, refTypeList, keepZeros, testMethodList, coreNum):
    
    #load test data
    print('##########loading test data')
    if testFormat == '10x':
        fileItem = glob.glob(os.path.join(testDS, "matrix.mtx"))[0]
        em = io.mmread(fileItem)
        em = em.tocsr().toarray()
        row = pd.read_table(fileItem[:-10]+"genes.tsv", header=None, index_col=None)
        col = pd.read_table(fileItem[:-10]+"barcodes.tsv", header=None, index_col=None)
        em = pd.DataFrame(em, index=row.T.values[1], columns=col.T.values[0])
        savefolder = testDS
    else:
        em = pd.read_csv(testDS, index_col=0, header=0)
        savefolder = testDS[:-4]
        
    # remove unrelated genes
    geneList = set(em.index)
    if testGenes != 'none':
        testGeneList = set(pd.read_csv(testGenes,index_col=None, header=0).columns)
        geneList = testGeneList.intersection(testGeneList)
        em = em.ix[~em.index.duplicated(keep='first'),]
        em = em.ix[geneList,]
    print('test dataset shape: %s genes, %s samples' % em.shape)
    
    hidSpecDict = {'9606':'human', '10090':'mouse'}
    
    for refType in refTypeList:
        print('##########annotating test data with %s data' % hidSpecDict[refType])
        for testMethod in testMethodList:
            AnnSCData(testType, em, refDS, refType, hidSpecDict[refType], keepZeros, testMethod, coreNum, savefolder)
            
    if len(refTypeList) == 2:
        print('##########merging annotation data')
        rstFolder = 'annotation_result'
        if keepZeros:
            rstFolder = rstFolder + '_keep_all_genes'
        else:
            rstFolder = rstFolder + '_keep_expressed_genes'
        savefolder = os.path.join(savefolder, rstFolder)
        for testMethod in testMethodList:
            tomergeList = glob.glob(os.path.join(savefolder, "human*%s*.xlsx" % testMethod))
            tomergeList = [i for i in tomergeList if "_Avg" not in i]
            topAnnList = []
            for tomerge in tomergeList:
                #merge results
                humanAnn = pd.read_excel(tomerge, header=[0,1,2])
                humanAnn.columns = pd.MultiIndex.from_arrays([list(humanAnn.columns.get_level_values(0)),list(humanAnn.columns.get_level_values(2))], names=('identifier', 'annotation'))
                mouseAnn = pd.read_excel(tomerge.replace('human_','mouse_'), header=[0,1,2])
                mouseAnn.columns = pd.MultiIndex.from_arrays([list(mouseAnn.columns.get_level_values(0)),list(mouseAnn.columns.get_level_values(2))], names=('identifier', 'annotation'))
                mergedAnn = pd.concat([humanAnn, mouseAnn])
                
                #sort merged results
                #split merged results to single-cell merged results
                mergedAnnList = np.split(mergedAnn, len(mergedAnn.columns)/2, axis=1)
                
                #annotate single-cell expression profiles in parallel
                p = multiprocessing.Pool(coreNum)
                resultList = p.map(SortAnno, mergedAnnList)
                p.close()
                p.join()
                
                mergedAnn = pd.concat([i[0] for i in resultList], axis=1)
                mergedAnn.to_excel(tomerge.replace('human_','combined_'))
                
                #save top ann result
                topAnn = pd.DataFrame({'cell':[i[1] for i in resultList], 'cell type':[''] * len(resultList), 'top sample':[i[2] for i in resultList], 'top correlation score':[i[3] for i in resultList]})
                topAnnList.append(topAnn)
            mergedAnn = pd.concat(topAnnList)
            mapData1 = pd.read_csv(os.path.join(refDS, '9606_map.csv'), index_col=0, header=0)
            mapData2 = pd.read_csv(os.path.join(refDS, '10090_map.csv'), index_col=0, header=0)
            mapData = pd.concat([mapData1, mapData2])
            for idx in mergedAnn.index:
                mergedAnn.ix[idx, 'cell type'] = mapData.ix[mergedAnn.ix[idx, 'top sample'], 'cell type']
            saveNameP = os.path.join(savefolder, "combined_%s_top_ann.csv" % (testMethod))
            mergedAnn.to_csv(saveNameP, index=False, columns = ['cell', 'cell type', 'top sample', 'top correlation score'])
            
        print('##########DONE!')

if __name__ == "__main__":
    
    #process arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--refType', default='human', help='human (default) | mouse | both ') 
    parser.add_argument('--testType', default='human', help='human (default) | mouse | rat | chimpanzee |zebrafish') 
    parser.add_argument('--refDS', required=True, help='path to the folder of reference dataset(s)')
    parser.add_argument('--dFormat', default='10x', help='10x (default) | csv')
    parser.add_argument('--testDS', required=True, help='path to the folder of test dataset if dtype is 10x, otherwise, the path to the file')
    parser.add_argument('--testMethod', default='s', help='s[pearman] (default) | p[earson] | both')
    parser.add_argument('--testGenes', default='none', help='optional, path to the csv file whose first row is the genes used to calculate the correlation coefficient')
    parser.add_argument('--keepZeros', default='yes', help='y[es] (default) | n[o]')
    parser.add_argument('--coreNum', type=int, default=1, help='number of the cores to use, default=1')
    
    opt = parser.parse_args()
    
    #check testType
    avaSpecDict = {'human':'9606', 'mouse':'10090', 'rat':'10116', 'chimpanzee':'9598', 'zebrafish':'7955'}
    if opt.testType.lower() not in avaSpecDict.keys():
        sys.exit("The testType can only be 'human', 'mouse', 'rat', 'chimpanzee' or 'zebrafish'.")
    else:
        testType = avaSpecDict[opt.testType]
        
    #check testMethod
    if opt.testMethod.lower() == "p":
        testMethodList = ['Pearson']
    elif opt.testMethod.lower() == "s":
        testMethodList = ['Spearman']
    elif opt.testMethod.lower() == "both":
        testMethodList = ['Pearson', 'Spearman']
    else:
        sys.exit("The testMethod can only be 's', 'p' or 'both'.")
    
    #check keepZeros
    if opt.keepZeros.lower()[0] == "n":
        keepZeros = False
    elif opt.keepZeros.lower()[0] == "y":
        keepZeros = True
    else:
        sys.exit("The keepZeros can only be 'y[es]' or 'n[o]'.")
    
    #check refType
    if opt.refType.lower() == "both":
        refTypeList = ['9606', '10090']
    elif opt.refType.lower() == "human":
        refTypeList = ['9606']
    elif opt.refType.lower() == "mouse":
        refTypeList = ['10090']
    else:
        sys.exit("The refType can only be 'human', 'mouse' or 'both'.")
        
    #check refDS
    if not os.path.exists(opt.refDS):
        sys.exit("The folder of reference dataset does not exist.")
    for refType in refTypeList:
        if refType == testType:
            if not os.path.exists(os.path.join(opt.refDS, '%s_symbol.csv' % refType)):
                sys.exit("The reference dataset '%s_symbol.csv' dose not exist." % refType)
        else:
            if not os.path.exists(os.path.join(opt.refDS, '%s_HID.csv' % refType)):
                sys.exit("The reference dataset '%s_HID.csv' dose not exist." % refType)
        
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
        if not os.path.exists(os.path.join(opt.testDS, 'genes.tsv')):
            sys.exit("Cannot find 'genes.tsv' file in the folder of test dataset.")
        if not os.path.exists(os.path.join(opt.testDS, 'barcodes.tsv')):
            sys.exit("Cannot find 'barcodes.tsv' file in the folder of test dataset.")
            
    #check gene list 
    if opt.testGenes != 'none' and opt.testGenes.endswith('.csv'):
        if not os.path.exists(opt.testGenes):
            sys.exit("Cannot find CSV file to obtain the test gene list.")
    
    #check coreNum
    maxCoreNum = multiprocessing.cpu_count()
    if opt.coreNum > maxCoreNum:
        sys.exit("There are only %s cores availble, less than %s cores." % (maxCoreNum, opt.coreNum))
        
    #pass argument check, show input data
    print('===================================================')
    print('Input data:')
    print('Reference species: %s' % opt.refType)
    print('The folder of reference dataset(s): %s' % opt.refDS)
    print('Test species: %s' % opt.testType)
    print('The format of the test dataset: %s' % opt.dFormat)
    if opt.dFormat == '10x':
        print('Test data are in the folder: %s' % opt.testDS)
    else:
        print('Test data are in the file: %s' % opt.testDS)
    print('The correlation method: %s' % opt.testMethod)
    if keepZeros:
        print('Keep genes with zero expreesion level: yes')
    else:
        print('Keep genes with zero expreesion level: no')
    if opt.testGenes != 'none':
        print('Test gene list: %s' % opt.testGenes)
    print('The number of cores to use: %s' % opt.coreNum)
    print('===================================================')
    
    #start to annotate test dataset
    main(testType, opt.dFormat, opt.testDS, opt.testGenes, opt.refDS, refTypeList, keepZeros, testMethodList, opt.coreNum)
    #python scMatch.py --refDS FANTOM5 --dFormat csv --testDS GSE81861_Cell_Line_COUNT.csv --coreNum 4
