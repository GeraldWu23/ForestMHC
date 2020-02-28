# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 16:49:46 2020

@author: wukak
"""
import random


PATH = 'D:/ForestMHC_origin/ForestMHC_master/outterdata/zeynap_input.txt'
with open(PATH, 'r') as f:
    data = f.readlines()
    random.shuffle(data)

with open(PATH, 'w') as f:
    f.writelines(data)
    