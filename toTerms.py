#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 10:37:14 2018

@author: rhou
"""
import warnings
warnings.filterwarnings("ignore")
import argparse
import pandas as pd
import numpy as np
import os, sys
import glob
import json
import multiprocessing
from functools import partial

def AvgTest(mapDict, barcodeDF):
    testID = barcodeDF.columns.get_level_values(0).values[0]
    barcodeDF.columns = ["spl", "score"]
    barcodeDF = barcodeDF.set_index('spl')
    termDict = {'Avg Score':[], 'Ont Term':[]}
    
    #get all sample names
    allSplSet = set(barcodeDF.index.values)
    
    for termItem in mapDict.keys():
        ffList = mapDict[termItem]
        if len(ffList) > 1:
            #number of top sample names in the list
            commonTerms = list(set(ffList).intersection(allSplSet))
            if len(commonTerms) > 0:
                totalCorr = 0.0
                for commonTerm in commonTerms:
                    totalCorr += barcodeDF.ix[commonTerm, 'score']
                avgScore = totalCorr/len(commonTerms)
                
                termDict['Ont Term'].append(termItem)
                termDict['Avg Score'].append(avgScore)
    
    termDF = pd.DataFrame.from_dict(termDict).ix[:, ['Ont Term', 'Avg Score']]
    termDF = termDF.sort_values(by=['Avg Score'], ascending=False)
    termDF = termDF.reset_index(drop=True)
    
    topScore = termDF.ix[:, 'Avg Score'][0]
    topAnns = list(termDF.ix[termDF['Avg Score']==topScore, 'Ont Term'].values)
    topAnns = sorted([i.split('.CNhs')[0] for i in topAnns])
    
    arrays = [[testID, testID], ['Ont Term', 'Avg Score']]
    tuples = list(zip(*arrays))
    headers = pd.MultiIndex.from_tuples(tuples, names=['identifier', 'annotation'])
    termDF.columns = headers
    
    return (termDF, (testID, ', '.join(topAnns), topScore))
    
#####average analysis
def MapAvg(currAnn, mapDict, splFile, refDS, spsType, coreNum):
    #split annotation file to single-cell vectors
    idIndex = range(0, len(currAnn.columns), 2)
    idIndex = currAnn.columns.get_level_values(0).values[idIndex]
    totalCounter = len(idIndex)
    annfiles = np.split(currAnn, totalCounter, axis=1)
    
    p = multiprocessing.Pool(coreNum)
    func = partial(AvgTest, mapDict)
    termAnnList = p.map(func, annfiles)
    p.close()
    p.join()
    
    merged = pd.concat([i[0] for i in termAnnList], axis=1)
    saveNameP = splFile[:-5]+"_Avg.xlsx"
    merged.to_excel(saveNameP)
    
    #save ann result
    merged = pd.DataFrame({'cell':[i[1][0] for i in termAnnList], 'avg annotation':[i[1][1] for i in termAnnList], 'top average correlation score':[i[1][2] for i in termAnnList]})
    saveNameP = splFile[:-5]+"_Avg_top_ann.csv"
    merged.to_csv(saveNameP, index=False, columns = ['cell', 'avg annotation', 'top average correlation score'])

#start to transfer original sample names
def main(splFileList, refDS, coreNum):
    for splFile in sorted(splFileList):
        print('#####processing %s' % splFile)
            
        if 'combined' in splFile:
            currAnn = pd.read_excel(splFile, header=[0,1])
        else:
            currAnn = pd.read_excel(splFile, header=[0,1,2])
            currAnn.columns = pd.MultiIndex.from_arrays([list(currAnn.columns.get_level_values(0)),list(currAnn.columns.get_level_values(2))], names=('identifier', 'annotation'))
        
        if 'human' in splFile:
            spsType = 'human'
            with open(os.path.join(refDS, 'human_samples_oto.txt')) as json_file:  
                mapDict = json.load(json_file)
        if 'mouse' in splFile:
            spsType = 'mouse'
            with open(os.path.join(refDS, 'mouse_samples_oto.txt')) as json_file:  
                mapDict = json.load(json_file)
        if 'combi' in splFile:
            spsType = 'hgmm'
            with open(os.path.join(refDS, 'hgmm_samples_oto.txt')) as json_file:  
                mapDict = json.load(json_file)
            
        #start test
        MapAvg(currAnn, mapDict, splFile, refDS, spsType, coreNum)
        print('#####DONE!')

if __name__ == "__main__":
    
    #process arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--splF', required=True, help='path to the origianl sample annotation folder which contains original sample annotation data')
    parser.add_argument('--refDS', required=True, help='path to the folder of reference dataset(s)')
    parser.add_argument('--coreNum', type=int, default=1, help='number of the cores to use, default is 1')
    
    opt = parser.parse_args()
    
    #check splFolder
    if os.path.isdir(opt.splF):
        splFileList = sorted(glob.glob(os.path.join(opt.splF, '*.xlsx')))
        splFileList = [i for i in splFileList if "Avg" not in i]
        if len(splFileList) == 0:
            sys.exit("Cannot find the original sample annotation folder.")
    else:
        sys.exit("Cannot find the original sample annotation folder.")
        
    #check refDS
    if not os.path.exists(opt.refDS):
        sys.exit("The folder of reference dataset does not exist.")    
                
    #check if the oto file exists
    speciesSet = set([os.path.basename(filename)[:5] for filename in splFileList])
    if 'human' in speciesSet:
        if len(glob.glob(os.path.join(opt.refDS, 'human_samples_oto.txt'))) == 0:
            sys.exit("The reference dataset folder's 'human_samples_oto.txt' dose not exist.")
    if 'mouse' in speciesSet:
        if len(glob.glob(os.path.join(opt.refDS, 'mouse_samples_oto.txt'))) == 0:
            sys.exit("The reference dataset folder's 'mouse_samples_oto.txt' dose not exist.")
    if 'combi' in speciesSet:
        if len(glob.glob(os.path.join(opt.refDS, 'hgmm_samples_oto.txt'))) == 0:
            sys.exit("The reference dataset folder's 'hgmm_samples_oto.txt' dose not exist.")
        
    #check coreNum
    maxCoreNum = multiprocessing.cpu_count()
    if opt.coreNum > maxCoreNum:
        sys.exit("There are only %s cores availble, less than %s cores." % (maxCoreNum, opt.coreNum))
    
    #pass argument check, show input data
    print('===================================================')
    print('Input data:')
    if os.path.isdir(opt.splF):
        print('The folder of original sample annotation data: %s' % opt.splF)
    else:
        print('The FANTOM5 sample annotation file: %s' % opt.splF)
    print('The folder of reference dataset(s): %s' % opt.refDS)
    print('The number of cores to use: %s' % opt.coreNum)
    print('===================================================')
    
    #start to transfer original sample names
    main(splFileList, opt.refDS, opt.coreNum)
    #python toTerms.py --splF GSE81861_Cell_Line_COUNT/annotation_result_keep_all_genes --refDS FANTOM5 --coreNum 4