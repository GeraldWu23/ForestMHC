# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 14:14:04 2020

@author: wukak
"""

with open('d:/ForestMHC_origin/ForestMHC_master/forests/mono_HASM/correct_output.txt') as correctf:
    correct = correctf.readlines()

with open('d:/ForestMHC_origin/ForestMHC_master/forests/mono_HASM/my_output.txt') as myoutputf:
    myoutput = myoutputf.readlines()
    
print(correct == myoutput)