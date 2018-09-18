# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 12:30:55 2018

@author: Kateka Seth
"""
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from change_biotapestry import remove_biotapestry
from change_biotapestry import add_biotapestry
from GetModel import get_model
from GetModel import convert_to_biotapestry
from Biotapestry import convert_biotapestry_to_antimony
import tellurium as te
import roadrunner
import antimony

csv_filename="../Biotapestry/8gene_broken.csv"
csv_newfile="../Biotapestry/8gene_test.csv"
#init_params = ['d_p', 'd_m' , 'L' , 'Vm' , 'a_p' , 'H', 'K']
add_biotapestry([(6,8,1),(3,7,-1)],csv_filename, csv_newfile)

ant_str = convert_biotapestry_to_antimony("../Biotapestry/8gene_test.csv", 8, [1/60, 1, 1/60, 1, 5/60, 5, 1/60])
r=te.loada(ant_str)
leaks=["time","P8","P5","P7","P4"]
protein=["time", "INPUT","P1","P2","P3","P4","P5","P6","P7","P8"]
mRNA=["time","mRNA1","mRNA2","mRNA3","mRNA4","mRNA5","mRNA6","mRNA7","mRNA8"]
r.simulate(0,300,100,mRNA)
r.plot(figsize=[7,7],xlim=(0,300),linewidth=2) 



connect1= list(range(1,8))
connect2= list(range(1,8))

#add

#exec("[(%d,8,%d),(%d,8,1)]",i,j,k)
#add_biotapestry([(1,8,1),(2,8,1)],csv_filename, csv_newfile)
#test
#ant_str2 = convert_biotapestry_to_antimony("../Biotapestry/8gene_test.csv", 8, [1]*7)


#remove
#can just use the original file for now
