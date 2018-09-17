# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 09:44:50 2018

@author: Kateka Seth
"""

from change_biotapestry import remove_biotapestry

remove=[(1,4),(3,5),(1,7),(1,8),(2,4),(7,7)]
csv_filename="Biotapestry/8gene_network.csv"
csv_newfile="Biotapestry/8gene_broken.csv"

remove_biotapestry(remove, csv_filename, csv_newfile)