# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 13:48:20 2018

@author: Kateka Seth
"""

from change_biotapestry import add_biotapestry

"""
given:
student's connections
broken csv file
name of new student's file
true csv file
"""
add=[(6,8,1),(3,5,-1),(6,7,-1),(7,7,-1)]
csv_broken="Biotapestry/8gene_broken.csv"
csv_network="Biotapestry/8gene_network.csv"
csv_student="Biotapestry/8gene_student.csv"

#add in student's connections
add_biotapestry(add, csv_broken, csv_student)

network = []

#strip useless stuff, add to set
#f = open(csv_network)
#for line in f: 
#    line = line.replace("\"", "")
#    line = line.strip()
#    words = line.split(",")
#    for i in range(len(words)):
#        words[i] = words[i].strip()
#    if "#" not in words[0] and words[0] != "model" and words[0] != "":
#        network.append(words)
#print(network)

f = open(csv_student)
f_true = open(csv_network)
for line in f:
    line = line.replace("\"", "")
    line = line.strip()
    print(line)
