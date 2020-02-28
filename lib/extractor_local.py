# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:29:00 2020

@author: wukak
"""

'''
this file is for extracting values(vector) in classification from attributes

'''

with open('D://ForestMHC_origin/ForestMHC_master/lib/aa.txt') as aafile:
    aali = aafile.read()
    aali = aali.split(',')

hydropathy_score = {
'G':-0.4, 'A':1.8, 'P':-1.6, 'V':4.2, 'L':3.8, 'I':4.5, 'M':1.9,
'F':2.8, 'Y':-1.3, 'W':-0.9, 'S':-0.8, 'T':-0.7, 'C':2.5, 'N':-3.5,
'Q':-3.5, 'K':-3.9, 'H':-3.2, 'R':-4.5, 'D':-3.5, 'E':-3.5
                    }

molar_mass = {
'G':75, 'A':89, 'P':115, 'V':117, 'L':131, 'I':131, 'M':149,
'F':165, 'Y':181, 'W':204, 'S':105, 'T':119, 'C':121, 'N':132,
'Q':146, 'K':146, 'H':155, 'R':174, 'D':133, 'E':147
              }

aromatics = ['F', 'Y', 'W']


def sparse(pep, aa):
    '''
    pep = a peptide
    '''
    onehot = []
    for a in pep:
        for r in aa:
            onehot.append(1 if a==r else 0)
    
    return onehot # length = (number of amino-acid) * 20


def hydropathy(pep, aa):
    '''
    pep = a peptide
    '''
    hydro_out = []
    for a in pep:
        hydro_out.append(hydropathy_score[a])
    
    return hydro_out

def mass(pep, aa):
    '''
    pep = a peptide
    '''
    mass_out = []
    for a in pep:
        mass_out.append(float(molar_mass[a])/(max(molar_mass.values())))
    
    return mass_out


def is_aromatics(pep, aa):
    '''
    pep = a peptide
    '''
    onehot = []
    for a in pep:
        onehot.append(1 if a in aromatics else 0)  # only answer whether there is aromatics with this amino-acid
    
    return onehot

    
class Extractor():
    def __init__(self, featlist, aa):
        self.featlist = featlist  # feature list that used to clf
        self.aa = aa  # amino-acid dict
        self.menu = {
                'hydropathy':hydropathy,
                'mass':mass,
                'is_aromatic':is_aromatics,
                'sparse':sparse
                }
        
        self.funcli = [self.menu.get(feat) if feat in self.menu.keys() else lambda a,b :"no such feature" for feat in featlist]
        
    def extract_from_list(self, pepls):
        '''
        pepls = peptide list
        
        '''
        output = []
        for pep in pepls:
            pepvector = []
            for func in self.funcli:
                pepvector.extend(func(pep, self.aa))
            output.append(pepvector)
        return output
    
    

if __name__ == '__main__':
    
    testExtractor = Extractor(['hydropathy','mass','sparse','is_aromatic'],aali)  # aali from extractor_local.py
    result_local = testExtractor.extract_from_list(['SADEPMTTF','DSETRRLYY','SSDRKGGSY'])  # 3*(9*(1+1+20+1))   


