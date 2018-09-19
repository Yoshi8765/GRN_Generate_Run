# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 10:46:56 2018

@author: Yoshi
"""

from RunModel import run_model
from GetModel import get_model

model_name='exp2'
numGenes = 8
seednum = 6372 #seed not working

# TODO: fix seed
# TODO: check to make sure all species in model move concentration

antStr,biotap_str=get_model(numGenes,model_name=model_name,seed=seednum,export=True)


#antStr= open('8gene_network.txt','r').read()
noiseLevel = 0.05 # put in a percentage. 0.05 = 5%
tmax=200 # The clean data will have tmax*5 datapoints
resolution = 10 # How much you want to divide from the clean data (ie 10 => tmax/10 datapoints)

r,res,resN = run_model(antStr,noiseLevel,exportData=[0,'M',True,True,True,True],inputData=[1,tmax,resolution],bioTap=biotap_str)

print('done!')