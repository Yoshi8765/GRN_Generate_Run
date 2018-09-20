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
#from data_analysis import interaction_estimate
import tellurium as te
import roadrunner
import antimony

# Round 2 of game (new network) analysis
csv_broken = "../playing_game/model_files/biotapestry_broken.csv"
csv_new = "../playing_game/model_files/test_added_connections.csv"

selections = ['time'] + ["mRNA" + str(i) for i in range(1,9)]

for i in range(1,9):
    add=[(i,2,-1)]
    print(i,j)
    add_biotapestry(add,csv_broken, csv_new)
    ant_str = convert_biotapestry_to_antimony(csv_new, 8, [0.01556653, 9.959682  , 0.1056418 , 6.66957033, 0.08160472,
                                                                                                4.25284957, 0.06687737])
                                                                 
    r=te.loada(ant_str)
    r.simulate(0,300,100,selections=selections)
    r.plot(figsize=[7,7],xlim=(0,300),linewidth=2,linestyle='-') 
    r.resetToOrigin()







#----------------------------------------------------

#csv_filename="../Biotapestry/8gene_broken.csv"
#csv_newfile="../Biotapestry/8gene_test.csv"
#csv_remove="../Biotapestry/8gene_remove.csv"
#init_params = ['d_p', 'd_m' , 'L' , 'Vm' , 'a_p' , 'H', 'K']
#add_biotapestry([(6,8,1),(3,5,-1),(6,7,-1),(7,7,-1),(2,4,-1),(6,4,-1)],csv_filename, csv_newfile)

#remove_biotapestry([(7,7),(3,1),(1,8),(1,7)],"../Biotapestry/8gene_network.csv",csv_remove)

#actual parameters=[1/60, 1, 1/60, 1, 5/60, 5, 1/60]
#small range=[8.4111E-05, 5.919E-01, 1.88693E-01, 4.035919E-01,5.38276E-02, 5.92, 2.9811E-02]
#large range=[3.13755E-05, 9.5372, 3.21424, 6.66695, 3.6342586E-02, 9.9475, 5.24815884E-02]

#ant_str = convert_biotapestry_to_antimony("../Biotapestry/8gene_test.csv", 8, [2.48117E-06, 5.907979, 1/60, 4.08598971,
#                                                                               4.6555711E-02, 9.965566, 1.69489E-02])
#r=te.loada(ant_str)
#leaks=["time","P8","P5","P7","P4"]
#protein=["time", "INPUT","P1","P2","P3","P4","P5","P6","P7","P8"]
#mRNA=["time","mRNA1","mRNA2","mRNA3","mRNA4","mRNA5","mRNA6","mRNA7","mRNA8"]
#r.simulate(0,300,100,mRNA)
#r.plot(figsize=[7,7],xlim=(0,300),linewidth=2) 

#csv_filename="../Biotapestry/8gene_test.csv"
#csv_newfile="../Biotapestry/8gene_test2.csv"

#for i in range(1,9):
#    for j in range(1,9):
#        add=[(i,4,-1),(j,4,-1)]
#        print(i,j)
#        r.resetToOrigin()
#        add_biotapestry(add,csv_filename,csv_newfile)
#        ant_str = convert_biotapestry_to_antimony("../Biotapestry/8gene_test2.csv", 8, [8.4111E-05, 5.919E-01, 1.88693E-01, 4.035919E-01,
#                                                                                        5.38276E-02, 5.92, 2.9811E-02])
#        r=te.loada(ant_str)
#        r.simulate(0,300,100,mRNA)
#        r.plot(figsize=[7,7],xlim=(0,300),linewidth=2,linestyle='-') 

