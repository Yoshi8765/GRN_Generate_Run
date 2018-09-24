# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 10:10:32 2018

@author: Kateka Seth
"""

from GetModel import get_model
from RunModel import run_model
import pandas as pd
#def get_model(num_genes, reg_probs = [0.2, 0.2, 0.2, 0.2, 0.2], model_name="pathway",init_params=[1.0/60,1,1.0/60,1,5.0/60,5,1.0/60],
#              param_std = 0.25, seed = 0, reachability=0.9, self_feedback_min = 0, max_builds = 1000, export = False):

data = get_model(8, seed=123, init_params=[10,1/60,1/60,1,5,3,1/100])
ant_str = data[0]
biotap = data[1]

#def run_model(antStr,noiseLevel,exportData=[ [0],'P',True,True,True],inputData=[1,40,1, [0],[0]],bioTap='',
#              savePath='\\model_output\\',showTimePlots=False,seed=0,drawModel=[False,'fdp'],runAttempts=5):

run_model(ant_str, noiseLevel=8, exportData= [[1,2,3,4,5,6,7,8], 'M', True,True,True], inputData=[1,200,10, [3], ["UP",20,5]])

data_table2 = pd.read_csv('Random_GRNs/pathway/Results_Clean.csv') # RNASeq_HiRes has timepoints = [0,200,20]
data_table2.set_index('time', inplace=True) 
data_table2.plot()
