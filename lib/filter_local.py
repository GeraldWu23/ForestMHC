# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:40:11 2020

@author: wukak
"""

class Filter:
    def __init__(self, removeC, aa, threshold = None, len_set = [9]):
        self.removeC = removeC
        self.aa = aa
        self.threshold = threshold
        self.len_set = len_set
        
        self.callcount = 0
        self.aacount = 0
        self.lencount = 0
        self.removeCcount = 0
        
    
    def may_pass(self, pep):
        self.callcount += 1
    
        for a in str(pep):
            #print(type(pep))
            if a not in self.aa:
                self.aacount += 1
                return False # forbidden amino acid
        
        if len(pep) not in self.len_set:
            self.lencount += 1
            return False # forbidden peptide length
        
        if self.removeC and 'C' in pep:
            self.removeCcount += 1
            return False # cystein exists
        else:
            return True
        
    
    
    def report(self):
        try:
            print('times of callings are ' + str(self.callcount))
            print(str(self.aacount) + " removed because of illegal amino acid,")
            print(str(self.removeCcount) + " removed because of containing cystein,")
            print(str(self.lencount) + " removed because of illegal length of peptides")
        except:
            print("LOCAL FILTER REPORT ERROR")
            return True            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        