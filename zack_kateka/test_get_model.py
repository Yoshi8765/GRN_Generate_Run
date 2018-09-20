# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 11:42:03 2018

@author: Kateka
"""

import tellurium as te
import time
import numpy as np

from GetModel import get_model
from GetModel import convert_to_biotapestry
from Biotapestry import convert_biotapestry_to_antimony


start = time.time()

ant_strs = []

for _ in range(10):
    next_antstring, biotap = get_model(8, seed = 6432, model_name="debug")
    #print (next_antstring)
    print ("-"*40)
    #ant_strs.append(next_antstring)








#ant_str,biotap = get_model(8, init_params=[1,2,3,4,5,0.5,7], model_name = "debug", param_std = 0, export=True)
#ant_str2 = convert_biotapestry_to_antimony("debug_biotapestry.csv",num_genes=8, init_params=[1,2,3,4,5,0.5,7], model_name="debug")

#print(ant_str)
#print("-"*40)
#print(ant_str2)

#print(ant_str == ant_str2)
#ant_str = convert_biotapestry_to_antimony("Biotapestry/8gene_network.csv", 8, [.5,.9,.8,30,30,1,.5])
#print(ant_str)
#r = te.loada(ant_str)
#leaks=["time","P8","P5","P7","P4"]
#protein=["time", "INPUT","P1","P2","P3","P4","P5","P6","P7","P8"]
#r.simulate(0,20,100,leaks)
#r.plot(figsize=[7,7],ylim=[0,300],linewidth=2)       

#ant_str = get_model(3, reg_probs=[0.3, 0.3, 0.2, 0.1, 0.1], reachability=1, init_params = [1]*7, param_std = 0, model_name = "test") 
#ant_str2 = convert_biotapestry_to_antimony("test_biotapestry.csv",3,[1]*7 , model_name = "test")

#print (ant_str)

#print ("-"*40)
#print ("\n\n\n")
#print (ant_str2)
#print ("-"*40)
#print (ant_str2 == ant_str)

#r=te.loada(ant_str)
#r.reset()
#res = r.simulate(0,50,1000)
#r.plot()
#r.draw(layout='fdp')

print ("\ndone!")
print ("runtime: " + str(time.time() - start))


