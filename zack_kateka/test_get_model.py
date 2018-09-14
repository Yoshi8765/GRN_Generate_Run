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

ant_str = get_model(20, reg_probs=[0.3, 0.3, 0.2, 0.1, 0.1], reachability= 1) 


#r=te.loada(ant_str)
#r.reset()
#res = r.simulate(0,50,1000)
#r.plot()
#r.draw(layout='fdp')

print ("\ndone!")
print ("runtime: " + str(time.time() - start))


