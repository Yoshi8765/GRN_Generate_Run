# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 09:44:50 2018

@author: Kateka Seth
"""

from change_biotapestry import remove_biotapestry
from change_biotapestry import add_biotapestry


csv_filename="Biotapestry/8gene_network.csv"
csv_filename2="Biotapestry/8gene_broken.csv"
csv_newfile="Biotapestry/8gene_broken3.csv"

#remove=[(1,4),(3,5),(1,7),(1,8),(2,4),(7,7)]
#remove_biotapestry(remove, csv_filename, csv_filename2)

# 1=positive, -1=negative
add = [(1,8,1),(3,8,1)]
add_biotapestry(add, csv_filename2, csv_newfile)



    
        