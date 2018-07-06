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

#%% Test model importing two models together

ffl1 = te.readFromFile('FFL1_ant.txt')
ffl2 = te.readFromFile('FFL2_ant.txt')

#model1 = te.tellurium.loadAntimonyModel('FFL_ant.txt')
#model2 = te.tellurium.loadAntimonyModel('FFL2_ant.txt')

b = te.loada(ffl1 + ffl2 + '''
model main
    var species p_con
    model1: ffl1();
    model2: ffl2();
    model1.p_output is p_con;
    model2.p_input is p_con;
end
''')
b.draw()

r = te.loada(ffl1 + '''
             model test1
             model1: ffl1();
             end
             ''')
r.draw()