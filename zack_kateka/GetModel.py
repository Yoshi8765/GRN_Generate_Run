import math
import numpy as np
import random
import time

# TO DO:
# **Decide how to fix issue with INPUT interfering with algorithm
#       could just count how many disjoint sets there are and re run algorithm if it happens to
#       fail due to this issue. Failure should be uncommon. I.E. add these lines at bottom of method
#       if (# disjointSets > 1) {
#           return get_model(num_genes, reg_probs, model_name, init_params)
#       }
# **Adjust sentinal values in DisjointSet to be -(total connections remaining + 1)
#   so that we avoid the issue of the sentinal values reaching 0 at the end. This will allow us
#   to easily count # of disjoint sets at end (simply count negative values in disjointSet array)


# init_params = ['d_p', 'd_m' , 'L' , 'Vm' , 'a_p' , 'K' , 'H']
def get_model(num_genes, reg_probs = [0.2, 0.2, 0.2, 0.2, 0.2], model_name="pathway", init_params=[0.5, 0.9, 0.8, 30, 30, 0.5, 1]):
    # cases/errors to handle:
    #     probs dont add to 1
    #     probs are all >= 0
    #     init_params has length not == 6
    #     num genes is > 1
    #     init_params are all greater than 0

    # Invalid parameter handling
    while num_genes < 2:
        num_genes = int(input("the number of genes cannot be negative, and zero genes or only one gene is not very fun. \n Please input a new integer: "))

    while not sum(reg_probs) == 1 or not len(reg_probs) == 5 or any(x < 0 for x in reg_probs):
        answer = input("The probabilities must sum to 1, and 5 positive probabilities must be given.\n Please enter 5 decimals separated only by spaces that add to 1: ")
        reg_probs = [int(s) for s in answer.split()]

    while not len(init_params) == 7 or any(x < 0 for x in init_params):
        answer = input("The initial parameter means must all be greater than zero, and there must be seven of them.\n Please enter 7 numbers separated only by spaces: ")
        init_params = [int(s) for s in answer.split()]

    while model_name == None:
        model_name = input("Please enter a model name")

    # Strategy: look through each gene. Based on its type, assign other proteins/genes
    # to act as its activators/repressors. Before assigning however, check if this process
    # will create an orphan (within each disjoint set, keep track of how many inputs are
    # left between all genes within the set; if the previous action connects two genes within
    # a set AND the # inputs remaining is only 1, do not connect the two)

    # watch out: if you add two elements from different sets, you add their remaining
    # inputs and then remove one input, but if they are in the same set already you only
    # remove one input (NO ADDING)

    gene_types = ["SA", "SR", "DA", "SA+SR","DR"]

    #np.random.seed(123) # for reproducibility (remove after testing)

    # TESTING: sampling works
    random_reg_types = np.random.choice(gene_types, size=num_genes, p=np.asarray(reg_probs), replace=True)

    # TESTING: gene creation works
    all_genes = []
    for i, reg_type in enumerate(random_reg_types):
        protein_name = "P" + str(i+1)
        all_genes.append(Gene(protein_name, reg_type))

    gene_sets = DisjointSet()
    for gene in all_genes:
        gene_sets.make_set(gene.protein_name, -1*(gene.remaining_connections + 1))

    # TESTING: Disjoint set appears to work
    #print (gene_sets)
    #print (gene_sets.find_set("P3"))
    #gene_sets.union("P1", "P4")
    #print (gene_sets)
    #gene_sets.decrement_value("P1")
    #print(gene_sets)


    assign_connections(all_genes, gene_sets)

    # TESTING: assignment protocol seems to work and generate no orphans
    # Notes:
    #for gene in all_genes:
    #    print (str(gene.protein_name) + "(" + str(gene.reg_type) + "): " + str(gene.in_connections))

    # Handles the case where INPUT cases algorithm to fail; this is only likely when the proportion
    # of double input genes (DA, DR, SA+SR) is low
    if (not gene_sets.get_set_count() == 1):
        ant_str = get_model(num_genes, reg_probs, model_name, init_params)
    else:
        ant_str = convert_to_antimony(all_genes, model_name, init_params)

    #print (ant_str)

    f = open(model_name + "_antimony.txt", 'w')
    f.write(ant_str)
    f.close()

    print ("done!")

    return ant_str




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
            elif gene_sets.get_total_connections(set1) > 1 or total_connections_left == 1:
                gene_sets.decrement_value(name1)
                gene.add_in_connection(gene_to_add)
                total_connections_left -= 1



def convert_to_antimony(all_genes, model_name, init_params):
    # add variation to parameters around their mean found in init_params
    std_dev_perc = 0.25

    ant_str = "\'\'\'\n"
    ant_str += "model *" + model_name + "()\n\n"
    ant_str += "\t// Compartments and Species:\n"

    species = "".join([gene.protein_name + ", " + "mRNA" + str(i+1) + ", " for i, gene in enumerate(all_genes)])
    ant_str += "\tspecies INPUT, " + species[:-2] + ";\n\n"

    rules = {}
    ant_str += "\t// Assignment Rules (production rates used in reactions):\n"
    for i,gene in enumerate(all_genes):
        type = gene.reg_type
        inputs = gene.in_connections
        P1 = inputs[0].protein_name
        P2 = ""
        if len(inputs) == 2:
            P2 = inputs[1].protein_name

        # TO DO: incorporate variable initial parameters

        #Ki_j corresponds to rate constant i for protein j
        K1 = "K1_" + str(i+1)
        K2 = "K2_" + str(i + 1)
        K3 = "K3_" + str(i + 1)
        H1 = "H" + str(i + 1)

        num = ""
        denom = ""
        if type == "SA":
            num = '(' + K1 + '*' + P1 + '^' + H1 + ')'
            denom = '(1 +' + K1 + '*' + P1 + '^' + H1 + ')'

        elif type == "SR":
            num = ' 1 '
            denom = '(1 +' + K1 + '*' + P1 + '^' + H1 + ')'

        elif type == "DA":
            num = '(' + K1 + '*' + P1 + '^' + H1 + ' + ' + K2 + '*' + P2 + '^' + H1 + ' + '+ K1 + '*' + K3 + '*' + P1 + '^' + H1 + '*' + P2 + '^' + H1 + ')'
            denom = '(1 + ' + K1 + '*' + P1 + '^' + H1 + ' + ' + K2 + '*' + P2 + '^' + H1 + ' + ' + K1 + '*' + K3 + '*' + P1 + '^' + H1 + '*' + P2 + '^' + H1 + ')'

        elif type == "DR":
            num = ' 1 '
            denom = '(1 + ' + K1 + '*' + P1 + '^' + H1 + ' + ' + K2 + '*' + P2 + '^' + H1 + ' + ' + K1 + '*' + K3 + '*' + P1 + '^' + H1 + '*' + P2 + '^' + H1 + ')'

        elif type == "SA+SR":
            num = '(' + K1 + '*' + P1 + '^' + H1 + ')'
            denom = '(1 + ' + K1 + '*' + P1 + '^' + H1 + ' + ' + K2 + '*' + P2 + '^' + H1 + ' + '+ K1 + '*' + K2 + '*' + P1 + '^' + H1 + '*' + P2 + '^' + H1 + ')'

        else:
            raise ValueError("gene type does not match any of the 5 standard varieties")

        expression = "Vm" + str(i+1) + '*(' + num  + '/' + denom + ')'
        rules["v" + str(i+1)] = expression
        ant_str += "\t// transcription" + str(i+1) + "("+str(gene.reg_type)+" : in connections = " + str(gene.in_connections)+ ")" + " uses production rate := " + expression +  ";\n"

    ant_str += "\n\t// Reactions:\n"
    for i in range(len(all_genes)):
        ant_str += "\ttranscription" + str(i+1) + ": => mRNA" + str(i+1) + " ; L" + str(i+1) + " + "+ rules["v" + str(i+1)] + " - d_mRNA"+ str(i+1) + " * mRNA" + str(i+1) + ";\n"
        ant_str += "\ttranslation" + str(i+1) + ": => P" + str(i+1) + " ; " + "a_protein"+ str(i+1) + " * mRNA" + str(i+1) + " - d_protein " + str(i+1) + " * protein" + str(i+1) + ";\n"



    ant_str += "\n\t// Species initializations:\n"
    ant_str += "\tINPUT = 1;\n"
    for i in range(len(all_genes)):
        ant_str += "\tmRNA" + str(i+1) + " = 0;\n"
        ant_str += "\tP" + str(i + 1) + " = 0;\n"


    ant_str += "\n\t// Variable initializations:\n"
    var_names = ["H", "Vm", "L", "d_mRNA", "a_protein", "d_protein", "K1_", "K2_", "K3_"]
    for i in range(len(all_genes)):
        for k, var in enumerate(var_names):
            mean = 0
            if k <= 5:
                mean = init_params[k]
            else:
                mean = init_params[6]
            std = mean * std_dev_perc
            value = np.random.normal(loc=mean, scale=std)
            ant_str += "\t" + var + str(i+1) + " = " + str(value) + ";\n"

    ant_str += "\n\t// Other declarations:\n"

    #variables = "".join(["v" + str(i+1) + ", " for i in range(len(all_genes))])
    #ant_str += "\tvar " + variables[:-2] + ";\n"
    ant_str += "\tconst "
    for i in range(len(all_genes)):
        for var in var_names:
            ant_str += var + str(i+1) + ", "
    ant_str = ant_str[:-2] + ";\n"


    ant_str += "\n\nend\n\'\'\'"
    return ant_str



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
        self.set_count = 0

    def make_set(self, item, value):
        if (self.contains(item)):
            raise ValueError("This item is already present in the DisjointSet")

        self.pointers.append(value)
        self.converter[item] = self.size
        self.size += 1
        self.set_count += 1

    # returns total # of connections remaining in set whose head is at the given index set_head
    def get_total_connections(self, head_index):
        return -1 * self.pointers[head_index]

    def get_set_count(self):
        return self.set_count

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
        self.set_count -= 1

    # used for the case where we want to form a connection between two genes that are already part of the same
    # network, but that still has some total connections remaining to avoid forming an orphan
    def decrement_value(self, item):
        # reduces total # of connections remaining in this set (uses +1 vs -1 as sentinal values are negative)
        head = self.find_set(item)
        self.pointers[head] += 1

    def __repr__(self):
        return str(self.pointers) + str(self.converter)



get_model(2)