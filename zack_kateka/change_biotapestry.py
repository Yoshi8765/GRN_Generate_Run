# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 09:38:56 2018

@author: Kateka Seth
"""
def remove_biotapestry(remove, csv_filename, csv_newfile):
    f = open(csv_filename)
    f_new = open(csv_newfile,'w')
    #count_node_only(f)
    og_nodes= set()
    new_sources= set()
    new_targets= set()
    i=0
    for line in f:
        line = line.replace("\"", "")
        words = line.split(",")
        if "#" not in words[0] and words[0] != "general":
            f_new.write(words[0])
            for j in range(1, len(words)):
                f_new.write(", " + words[j])
        elif words[0] == "general":
            gene_source=words[3]
            gene_target=words[5]
            og_nodes.add(gene_source)
            og_nodes.add(gene_target)
            given_source="Gene " + str(remove[i][0])
            if remove[i][0] == "INPUT":
                given_source="INPUT"
            given_target="Gene " + str(remove[i][1])
            if (given_source == gene_source) and (given_target == gene_target):
                i += 1
            else:
                f_new.write(words[0])
                for j in range(1, len(words)):
                    f_new.write(", " + words[j])
                new_sources.add(gene_source)
                new_targets.add(gene_target)
    # print the node only
    node_only=og_nodes.difference(new_sources).difference(new_targets)
    for node in node_only:
        f_new.write("nodeOnly, root, gene, " + node + "\n")
    f.close()
    f_new.close()