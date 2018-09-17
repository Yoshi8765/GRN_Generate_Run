# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 10:46:56 2018

@author: Yoshi
"""

from RunModel import run_model
from GetModel import get_model

antStr= open('8gene_network.txt','r').read()
antStr,biotap_str=get_model(8)

noiseLevel = 0.15 # put in a percentage. 0.05 = 5%
r,res,resN = run_model(antStr,noiseLevel,exportData=[ 0,'M','y','y','y','y','y'],biotap=biotap_str)