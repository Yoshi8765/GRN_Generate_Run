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
csv_filename="Biotapestry/8gene_broken.csv"
csv_newfile="Biotapestry/8gene_test.csv"

test_filename="Biotapestry/test_network.csv"

#add_biotapestry([(5,7,1), (5,7,1), (4,3,1), (3,1,1)],csv_filename2,test_filename)


def write_fencepost(f_new, words):
    f_new.write(words[0].strip())
    for i in range(1, len(words)):
        f_new.write(", " + words[i].strip())
    f_new.write("\n")
    
#init_params = ['d_p', 'd_m' , 'L' , 'Vm' , 'a_p' , 'H', 'K']
#add_biotapestry([(1,7,1)],csv_filename, csv_newfile)
add=[(1,7,1)]

f = open(csv_filename)
f_new = open(csv_newfile,'w')
new_nodes = set()


# If a nodeOnly is inside new_nodes, erase the nodeOnly line
for i in add:
    if i[0] == "INPUT" or i[1] == "INPUT":
        new_nodes.add("INPUT");
    else:
        new_nodes.add("Gene " + str(i[0]))
        new_nodes.add("Gene " + str(i[1]))
    
# erase nodeOnly
for line in f:
    line = line.replace("\"", "")
    words = line.split(",")
    if "Instance" not in words[0].strip() and words[0].strip() != "nodeOnly":
        write_fencepost(f_new, words)
    elif words[0].strip() == "nodeOnly":
        if words[3].strip() not in new_nodes:
            write_fencepost(f_new, words)

# write new nodes
for i in range(0, len(add)):
    f_new.write("general, root, ")
    for j in range(0,2):
        if add[i][j] == "INPUT":
            f_new.write("box, ")
        else:
            f_new.write("gene, ")
        f_new.write("Gene " + str(add[i][j]) + ", ")    
    if add[i][2] == 1:
        f_new.write("positive\n")
    else:
        f_new.write("negative\n")
f.close()
f_new.close()

