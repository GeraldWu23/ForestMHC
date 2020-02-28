# ForestMHC 
# Copyright 2018, Kevin Michael Boehm
# Use is subject to license as agreed to upon download
# Englander Institute for Precision Medicine, Weill Cornell Medical College

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
    
    with open(os.path.join(PATH, 'lib', 'aa.txt')) as aa_f:
        AA = aa_f.read().split(',')  # AA = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'] standard amino acids
    orig_peptides = in_list  # input raw peptides from input files(a line as a peptide)

    # remove peptides of len != 9 or with non-standard amino acids
    fil = Filter(removeC=False, aa=AA, len_set=[9])
    print(AA)
    peptides = []
    failed_peptides = []
    for peptide in orig_peptides:
        if fil.may_pass(peptide):
            peptides.append(peptide)
        else:
            failed_peptides.append(peptide)
    del orig_peptides, fil
    print(len(peptides), len(failed_peptides))
    
    # extract features and score each peptide
    features_to_use = ['hydropathy','is_aromatic','sparse','mass']
    extractor = Extractor(features_to_use, AA)
    features = extractor.extract_from_list(peptides)  # a list of featrows with each featrows of a vector for a peptide
                                                      # it represents every peptides    
    del extractor


    # prediction
    scores = []
    removal_list = []
    for raw_allele in alleles:
        allele = condition_allele(raw_allele)  # remove formatting in HLA allele(such as '-', '*', 'HLA', ':')
        clf_fn = os.path.join(PATH, 'forests', 'mono_HASM', '.'.join([allele, 'pkl']))
        try:
            clf = joblib.load(clf_fn)  # type(class sklearn.ensemble.forest.RandomForestClassifier)
            print(allele)
        except IOError:
            print('No predictions available for {}'.format(raw_allele))
            removal_list.append(raw_allele)  # no clf for this allele
            continue
        score = clf.predict_proba(features)  # npeptides * (not-binding, binding)
        scores.append(zip(*score)[1])  # nalleles * npeptides
        
        
        
    # end script if no predictions available for any alleles
    for allele in removal_list:
        alleles.remove(allele)
    if alleles == []:
        return None, None, None
        
    # predict corresponding allele according to score
    predicted = []
    for scores_for_peptide in zip(*scores):  # scores for peptide contains list of peptides with n prediction score for each allele
        ind = np.argmax(scores_for_peptide)
        predicted.append(alleles[ind])  # allele corresponding to the best score
    # sort by maximal probability
    all_data = zip(peptides, predicted, *scores)
    maxes = [max(el[2:]) for el in all_data]  # highest probability score of alleles
    temp = zip(maxes, all_data)
    temp.sort(reverse=True)
    sorted_data = [el[1] for el in temp]

    return alleles, sorted_data, failed_peptides




if __name__ == '__main__':
    in_list = get_in_list('d:/ForestMHC_origin/ForestMHC_master/forests/mono_HASM/input')
    alleles, sorted_data, failed_peptides_ori = forest_mhc(['A0101', 'C0501'], get_in_list('d:/ForestMHC_origin/ForestMHC_master/forests/mono_HASM/input'))