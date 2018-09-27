# -*- coding: utf-8 -*-
"""
Created on Thu Jul 05 12:28:35 2018

@author: Yoshi
"""

import os
import tellurium as te
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(linewidth=160)

#%% Visualize the models used in DREAM7

r = te.loadAntimonyModel('C:\Users\Yoshi\Documents\GitHub\DREAM-work\model1_0\model_without_parameters\model1_antimony.txt')

r.draw()