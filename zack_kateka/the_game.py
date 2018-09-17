# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 12:30:55 2018

@author: Kateka Seth
"""

from change_biotapestry import remove_biotapestry
from change_biotapestry import add_biotapestry
from GetModel import get_model
from GetModel import convert_to_biotapestry
from Biotapestry import convert_biotapestry_to_antimony
import tellurium as te
import roadrunner
import antimony

ant_str = convert_biotapestry_to_antimony("Biotapestry/8gene_broken.csv", 8, [1]*7)
print(ant_str)
r=te.loada(ant_str)
leaks=["time","P8","P5","P7","P4"]
protein=["time", "INPUT","P1","P2","P3","P4","P5","P6","P7","P8"]
r.simulate(0,20,100,protein)
r.plot(figsize=[7,7],ylim=[0,2],linewidth=2) 

