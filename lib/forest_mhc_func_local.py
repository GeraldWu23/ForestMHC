# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 13:46:31 2020

@author: wukak
"""

import sys
import os
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier
import numpy as np
try:
    from lib.extractor_local import Extractor
    from lib.filter_local import Filter
    from lib.auxiliary import condition_allele, get_in_list
except:
    from extractor_local import Extractor
    from filter_local import Filter
    from auxiliary import condition_allele, get_in_list
import pandas


#######
##
## user must insert path for ForestMHC main directory below
PATH = '/ForestMHC_origin/ForestMHC_master'
##
##
#######

def forest_mhc(alleles, in_list):
    '''
    alleles = a list of alleles
    in_list = input peptides
    
    '''
    
    orig_peptides = in_list
    
    # filter other illegal peptide
    with open(os.path.join(PATH, 'lib', 'aa.txt')) as AAf:
        AA = AAf.read().split(',')  # legal amino acid list
        AA = [aa[:1] for aa in AA]  # get rid of '\n' and other characters

    peptides = []  # legal peptides list
    fail_peptides = []  # illegal peptides list
    filt = Filter(removeC = False, aa = AA, len_set = [8,9,10,11])
    print(AA)
    for peptide in orig_peptides:
        if filt.may_pass(peptide):
            peptides.append(peptide)
            if peptide == 'LVDGYLNTY':
                print(filt.may_pass('LVDGYLNTY'))
        else:
            if peptide == 'LVDGYLNTY':
                print('fail peptides')
            fail_peptides.append(peptide)
    del filt, orig_peptides
    print(len(peptides), len(fail_peptides))
    
    # encode the peptides to vectors
    feature_list = ['hydropathy','is_aromatic','sparse','mass']  # order of the features(the format to encode a peptide) must be the same as that in the training data of classifier
    extractor = Extractor(feature_list, AA)
    vectors = extractor.extract_from_list(peptides)  
    del extractor
    
    
    # predict scores for peptides allele by allele
    allele_list = []  # processed alleles
    fail_allele = []  # failly processed alleles
    allele_scores = []
    for raw_allele in alleles:
        allele = condition_allele(raw_allele)
        try:
            clf = joblib.load(os.path.join(PATH, 'forests', 'mono_HASM', (str(allele)+'.pkl')))  # load classifier from local file
            allele_list.append(allele)
            score = np.array(clf.predict_proba(vectors)).T[1]  # predict the score for each vector accroding to allele
            allele_scores.append(score)  # the scores of peptides for an allele are appended as a list
            print('processed ' + str(allele))  # print allele was being processing
        except:
            print('no classifier for {}.'.format(raw_allele))
            fail_allele.append(raw_allele)
    
    allele_scores = np.array(allele_scores).T  # a peptide per line       
    
        
    # end script if no allele actually has been processed
    if fail_allele == alleles:
        print('no allele has been processed')
        return None, None, None

    
    # predict which allele this a peptide binds to(according to score)
    allele_predicted = [allele_list[np.argmax(s)] for s in allele_scores]
        
    # form all info
    all_data = zip(peptides, allele_predicted, *allele_scores.T)
    sorted_data = sorted(all_data, key = lambda x:max(x[2:]), reverse = True)    
    
        
    return allele_list, sorted_data, fail_peptides






if __name__ == '__main__':
#    import time
#    start = time.time()
#    in_list = get_in_list('d:/ForestMHC_origin/ForestMHC_master/forests/mono_HASM/input')
#    allele, data, fail_peptides = forest_mhc(['A0101', 'C0501'], get_in_list('d:/ForestMHC_origin/ForestMHC_master/forests/mono_HASM/input'))
#    print(time.time() - start)
    
    # get allele list
    with open('D:/ForestMHC_origin/ForestMHC_master/outterdata/zeynap_allele.txt', 'r') as f_allele:
        alleles_raw = f_allele.readlines()
        alleles = list(set([condition_allele(rawallele) for rawallele in alleles_raw]))
    
    # get peptides and labels
    peptides = []
    true_labels = []
    with open('D:/ForestMHC_origin/ForestMHC_master/outterdata/zeynap_input.txt', 'r') as f_pep:
        lines = f_pep.readlines()
        for line in lines:
            try:
                peptide, lab = line.split()
            except:
                print(line)
            peptides.append(peptide)
            true_labels.append(lab)
    
    # classification
    allele_list, sorted_data, fail_peptides = forest_mhc(alleles, peptides)
    pep_dict = {}
    for line in sorted_data:
        pep = line[0]
        pep_dict[pep] = max(line[2:])
    
    

    
    
    
    
    