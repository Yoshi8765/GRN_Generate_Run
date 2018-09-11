import math
import numpy as np
import random
import time


## need to include initial parameters

def get_model(num_genes, reg_probs = [0.2, 0.2, 0.2, 0.2, 0.2]):
    #TO DO: Implement parameter for initial parameter values (rate constants etc)

    if not len(reg_probs) == 5:
        raise ValueError("There are 5 gene types, but " + str(len(reg_probs)) + " probabilities were given")

    # Strategy: look through each gene. Based on its type, assign other proteins/genes
    # to act as its activators/repressors. Before assigning however, check if this process
    # will create an orphan (within each disjoint set, keep track of how many inputs are
    # left between all genes within the set; if the previous action connects two genes within
    # a set AND the # inputs remaining is only 1, do not connect the two)

    # watch out: if you add two elements from different sets, you add their remaining
    # inputs and then remove one input, but if they are in the same set already you only
    # remove one input (NO ADDING)

    gene_types = ["SA", "SR", "DA", "SA+SR","DR"]

    np.random.seed(0) # for reproducibility (remove after testing)

    # TESTING: sampling works
    random_reg_types = np.random.choice(gene_types, size=num_genes, p=np.asarray(reg_probs), replace=True)

    # TESTING: gene creation works
    all_genes = []
    for i, reg_type in enumerate(random_reg_types):
        protein_name = "P" + str(i+1)
        all_genes.append(Gene(protein_name, reg_type))

    gene_sets = DisjointSet()
    for gene in all_genes:
        gene_sets.make_set(gene.protein_name, -1*gene.remaining_connections)

    # TESTING: Disjoint set appears to work
    #print (gene_sets)
    #print (gene_sets.find_set("P3"))
    #gene_sets.union("P1", "P4")
    #print (gene_sets)
    #gene_sets.decrement_value("P1")
    #print(gene_sets)


    assign_connections(all_genes, gene_sets)

    # TESTING: assignment protocol
    for gene in all_genes:
        print (str(gene.reg_type) + "(" + str(gene.protein_name) + "): " + str(gene.in_connections))

    print "done!"




# Notes: a protein cannot be connected to same gene more than once
def assign_connections(all_genes, gene_sets):
    # randomly choose a gene to have INPUT
    input_gene = np.random.choice(all_genes)
    input_gene.add_in_connection(Gene("INPUT"))
    gene_sets.decrement_value(input_gene.protein_name)

    total_connections_left = sum([gene.remaining_connections for gene in all_genes])

    for gene in all_genes:
        while (gene.remaining_connections > 0):
            available_genes = list(set(all_genes) - set(gene.in_connections))
            gene_to_add = np.random.choice(available_genes)

            # check if new connection is valid
            name1 = gene.protein_name
            name2 = gene_to_add.protein_name

            set1 = gene_sets.find_set(name1)
            set2 = gene_sets.find_set(name2)
            if not set1 == set2:
                gene_sets.union(name1, name2)
                gene.add_in_connection(gene_to_add)
                total_connections_left -= 1
            # else genes are part of same set, but is there still enough connections remaining to form
            # another without creating an orphan?

            # Second condition handles edge case where final connection is trying to be made.
            # In this case, the situation is similar to when an orphan is being formed. However,
            # in this case the "orphan" is the entire, complete network
            elif gene_sets.get_size(set1) > 1 or total_connections_left == 1:
                gene_sets.decrement_value(name1)
                gene.add_in_connection(gene_to_add)
                total_connections_left -= 1



def convert_to_antimony(all_genes):
    # protein/mRNA start at zero concentration
    return "not implemented"



class Gene():

    def __init__(self, protein_name, reg_type = None):
        # name of protein created by this gene
        self.protein_name = protein_name

        # regulation type this gene undergoes
        self.reg_type = reg_type

        # list of genes connected to this one (genes whose proteins are activators/repressors)
        self.in_connections = []
        if self.protein_name == "INPUT":
            self.remaining_connections = 0
        elif self.reg_type == "SA" or self.reg_type == "SR":
            self.remaining_connections = 1
        else:
            self.remaining_connections = 2

    def add_in_connection(self, other):
        if self.remaining_connections <= 0:
            raise ValueError("Too many connections have been added to a gene of this type")
        self.in_connections += [other]
        self.remaining_connections -= 1

    def equals(self, other):
        return self.protein_name == other.protein_name

    # toString method mostly used for debugging
    def __repr__(self):
        return str(self.protein_name)
        #return (str(self.reg_type) + " (" + str(self.protein_name) + "): " + str(self.remaining_connections) + " connection(s) remaining")




# TO DO: provide link to CSE 373 page for explanation of disjoint set
class DisjointSet():

    def __init__(self):
        # stores values that point towards head of disjoint set, or stores a sentinal value representing
        # total # of remaining connections within that set if that element is the head of a disjoint set
        self.pointers = []

        # allows us to find index in array that stores value associated to given gene
        self.converter = {}

        # number of total elements in all disjoint sets
        self.size = 0

    def make_set(self, item, value):
        if (self.contains(item)):
            raise ValueError("This item is already present in the DisjointSet")

        self.pointers.append(value)
        self.converter[item] = self.size
        self.size += 1

    # returns total # of connections remaining in set whose head is at the given index set_head
    def get_size(self, head_index):
        return -1 * self.pointers[head_index]

    def find_set(self, protein_name):
        index = self.converter.get(protein_name)
        return self.find_helper(index)

    # returns index of representative element of disjoint set
    def find_helper(self, index):
        if (self.pointers[index] < 0):
            return index
        else:
            representative_index = self.find_helper(self.pointers[index])
            self.pointers[index] = representative_index
            return representative_index

    def contains(self, item):
        return item in self.converter.keys()

    # combines the two sets of item1 and item2; item1 becomes the new head
    def union(self, item1, item2):
        if not self.contains(item1) or not self.contains(item2):
            raise ValueError("One of those items given is not part of this disjoint set collection")
        head1 = self.find_set(item1)
        head2 = self.find_set(item2)
        if head1 == head2:
            raise ValueError("These two items are already part of the same set")

        self.pointers[head1] += self.pointers[head2] + 1 # +1 is due to 1 connection lost as new connection forms
        self.pointers[head2] = head1

    # used for the case where we want to form a connection between two genes that are already part of the same
    # network, but that still has some total connections remaining to avoid forming an orphan
    def decrement_value(self, item):
        # reduces total # of connections remaining in this set (uses +1 vs -1 as sentinal values are negative)
        head = self.find_set(item)
        self.pointers[head] += 1

    def __repr__(self):
        return str(self.pointers) + str(self.converter)



get_model(10)