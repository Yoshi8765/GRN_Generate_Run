# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 10:46:56 2018

@author: Yoshi
"""

from RunModel import run_model

antStr = open('pathway_antimony.txt','r').read()
noiseLevel = 0.15 # put in a percentage. 0.05 = 5%
r,res,resN = run_model(antStr,noiseLevel,exportData=[ [2,3,4],'P','y','y','y','y'],showTimePlots='y',drawModel=['y','fdp'])