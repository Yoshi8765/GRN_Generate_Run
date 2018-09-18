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

ant_str = convert_biotapestry_to_antimony("../Biotapestry/8gene_broken.csv", 8, [1]*7)
r=te.loada(ant_str)
leaks=["time","P8","P5","P7","P4"]
protein=["time", "INPUT","P1","P2","P3","P4","P5","P6","P7","P8"]
r.simulate(0,20,100,leaks)
r.plot(figsize=[7,7],ylim=[0,300],linewidth=2) 

csv_filename="../Biotapestry/8gene_broken.csv"
csv_newfile="../Biotapestry/8gene_test.csv"

connect1= list(range(1,8))
connect2= list(range(1,8))

#add
add_biotapestry([(1,8,1),(2,8,1)],csv_filename, csv_newfile)
#test
ant_str2 = convert_biotapestry_to_antimony("../Biotapestry/8gene_test.csv", 8, [1]*7)


#remove
#can just use the original file for now
