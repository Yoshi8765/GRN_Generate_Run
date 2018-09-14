# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 11:42:03 2018

@author: Kateka
"""

import tellurium as te
import time

from GetModel import get_model
from GetModel import convert_to_biotapestry
from Biotapestry import convert_biotapestry_to_antimony


start = time.time()

ant_str = convert_biotapestry_to_antimony("Biotapestry/8gene_broken.csv", 8, [1]*7)

print (ant_str)
       

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


