import numpy as np
import time


# Generates and returns an antimony string for a random biological pathway involving num_genes genes.
# The pathway is fully connected, and contains no orphans. Also, saves a plain text file
# containing the antimony string to a plain text file.

# init_params = ['d_p', 'd_m' , 'L' , 'Vm' , 'a_p' , 'K' , 'H']
# reg_probs = [prob(SA), prob(SR), prob(DA), prob(SA+SR), prob(DR)]
def get_model(num_genes, reg_probs = [0.2, 0.2, 0.2, 0.2, 0.2], model_name="pathway", init_params=[0.5, 0.9, 0.8, 30, 30, 0.5, 1],seed = 0):
    """Docstring for help command: Demonstrate docstrings and does nothing really."""
    # Invalid parameter handling

    #catch for illegal file names
    forbiddenChar = ['<','>',':','"','/','\'','|','?','*']
    try:
        assert list(set(model_name).intersection(forbiddenChar))==[]
    except AssertionError:
        raise AssertionError('You used an illegal character in your model name!')

    if not type(num_genes) == int or num_genes < 2:
        raise ValueError("num_genes is invalid: it must be a single integer greater than 1")

    if not type(reg_probs) == list or any([not (type(x) == int or type(x) == float) for x in reg_probs]) or not sum(reg_probs) == 1 or not len(reg_probs) == 5 or any(x < 0 for x in reg_probs):
        raise ValueError("reg_probs is invalid: reg_probs must be a list of 5 positive decimals that sum to 1")

    if not type(init_params) == list or any([not (type(x) == int or type(x) == float) for x in init_params]) or not len(init_params) == 7 or any(x < 0 for x in init_params):
        raise ValueError("init_params is invalid: init_params must be a list of 7 positive numbers")

    if not type(model_name) == str or model_name == None or len(model_name) == 0:
        raise ValueError("model name is invalid: must be a valid string with no illegal characters")

    # Algorithm: look through each gene. Based on its type, assign other proteins/genes
    # to act as its activators/repressors. Before assigning however, check if this process
    # will create an orphan (within each disjoint set, keep track of how many inputs are
    # left between all genes within the set; if the previous action connects two genes within
    # a set AND the # inputs remaining is only 1, do not connect the two). The INPUT messes with
    # this algorithm because it acts as a gene with no in-connections. Thus, if it initially connects
    # to a single-regulation type (SA or SR) it may create an orphan, as there is no guarantee
    # this SA/SR will connect to other genes, and since it used up its only in connection with
    # INPUT, it might not be connected at all. This is likely only an issue if the proportion of
    # the double input genes (DA, DR, SA+SR) is significantly low

    gene_types = ["SA", "SR", "DA", "SA+SR","DR"]

    if seed!=0:
        np.random.seed(seed)
    try:
        assert seed >= 0 and seed < 0xffffffff
    except AssertionError:
        raise AssertionError("Please use a 32-bit unsigned integer or 0.")


    # chooses a random collection of gene types of the required size with the given probabilities
    random_reg_types = np.random.choice(gene_types, size=num_genes, p=np.asarray(reg_probs), replace=True)

    # assigns genes their names and regulation types
    all_genes = []
    for i, reg_type in enumerate(random_reg_types):
        protein_name = "P" + str(i+1)
        all_genes.append(Gene(protein_name, reg_type))

    # DisjointSets helps keeps track of which genes are connected, directly or indirectly.
    # Used to prevent the formation of "orphans"
    gene_sets = DisjointSets()
    for gene in all_genes:
        gene_sets.make_set(gene.protein_name, -1*(gene.remaining_connections + 1))

    assign_connections(all_genes, gene_sets)

    #for gene in all_genes:
    #    print (str(gene.protein_name) + "(" + str(gene.reg_type) + "): " + str(gene.in_connections))

    # Handles the case where INPUT causes algorithm to fail; this is only likely when the proportion
    # of double input genes (DA, DR, SA+SR) is low
    if (not gene_sets.get_set_count() == 1):
        ant_str = get_model(num_genes, reg_probs, model_name, init_params)
    else:
        ant_str = convert_to_antimony(all_genes, model_name, init_params)

    f = open(model_name + "_antimony.txt", 'w')
    f.write(ant_str)
    f.close()

    return ant_str



# Assigns all in connections for each gene (i.e. assigns proteins which act as regulators for each gene)
# params
#  ** all genes = a list of all the genes in the network
#  ** gene_sets = a collection of disjoint sets, keeping track of which genes are already connected

# Assumptions:
# ** a protein cannot be connected to same gene more than once
# ** connections should not be formed if they result in the formation of an orphan
def assign_connections(all_genes, gene_sets):

    # randomly choose a gene to have INPUT as one of its regulators
    input_gene = np.random.choice(all_genes)
    input_gene.add_in_connection(Gene("INPUT"))
    gene_sets.decrement_value(input_gene.protein_name)

    total_connections_left = sum([gene.remaining_connections for gene in all_genes])

    # go gene to gene and fill in any empty regulatory sites
    for gene in all_genes:

        # while the gene still has unoccupied/unassigned regulatory sites
        while (gene.remaining_connections > 0):

            # to avoid connecting same protein to both regulatory sites on a gene, do not consider
            # genes whose proteins are already acting as a regulator for this current gene
            available_genes = list(set(all_genes) - set(gene.in_connections))
            gene_to_add = np.random.choice(available_genes)

            # check if new connection is valid (does not result in orphan)
            name1 = gene.protein_name
            name2 = gene_to_add.protein_name

            set1 = gene_sets.find_set(name1)
            set2 = gene_sets.find_set(name2)

            # If the two genes are not already connected, connect them
            if not set1 == set2:
                gene_sets.union(name1, name2)
                gene.add_in_connection(gene_to_add)
                total_connections_left -= 1

            # else genes are part of same set (already connected), but is there
            # still enough connections remaining to form another without creating an orphan?

            # Second condition handles edge case where final connection is trying to be made.
            # In this case, the situation is similar to when an orphan is being formed. You are
            # connecting two genes in same network, and only one connection remains. However,
            # the "orphan" is the entire, complete network
            elif gene_sets.get_total_connections(set1) > 1 or total_connections_left == 1:
                gene_sets.decrement_value(name1)
                gene.add_in_connection(gene_to_add)
                total_connections_left -= 1


# converts the generated network to an antimony string
def convert_to_antimony(all_genes, model_name, init_params):
    # add variation to parameters around their mean found in init_params
    std_dev_perc = 0.25

    ant_str = ""
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
        ant_str += "\t//transcription" + str(i+1) + "\n"
        ant_str += "\tJ" + str(i+1) + ": => mRNA" + str(i+1) + " ; L" + str(i+1) + " + "+ rules["v" + str(i+1)] + " - d_mRNA" + str(i+1) + " * mRNA" + str(i+1) + ";\n"
        ant_str += "\t//translation" + str(i+1) + "\n"
        ant_str += "\tF" + str(i+1) + ": => P" + str(i+1) + " ; " + "a_protein" + str(i+1) + " * mRNA" + str(i+1) + " - d_protein" + str(i+1) + " * P" + str(i+1) + ";\n"



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


    ant_str += "\n\nend"
    return ant_str



# Keeps track of the gene regulatory type, the protein name associated to this gene,
# the other genes this gene has already been connected to, and the remaining connections to be made
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

    # connect another gene to this one (add another regulator)
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



# DisjointSets keeps a collection of disjoint sets. In this code, these sets represent collections of
# genes that are connected, directly or indirectly, to each other (when the network is considered as an
# undirected graph). Helps keep track of conditions important for avoiding the creation of orphans.

# The following link gives some more information on disjoint sets (pg 32-57)
# https://courses.cs.washington.edu/courses/cse373/18su/files/slides/lecture-19.pdf
class DisjointSets():

    def __init__(self):
        # stores values that point towards index of head of disjoint set, or stores a sentinel value representing
        # total # of remaining connections within that set if that element is the head of a disjoint set
        self.pointers = []

        # allows us to find index in array that stores value associated to given gene (keys are protein names)
        self.converter = {}

        # number of total elements in all disjoint sets
        self.size = 0

        # count of the number of disjoint sets in the collection
        self.set_count = 0

    # makes a new set of a single element, and stores the given value as a sentinel value
    def make_set(self, item, value):
        if (value >= 0):
            raise ValueError("The sentinel value must be less than zero")
        if (self.contains(item)):
            raise ValueError("This item is already present in the DisjointSet")

        self.pointers.append(value)
        self.converter[item] = self.size
        self.size += 1
        self.set_count += 1

    # returns total # of connections remaining in set whose head is at the given index set_head
    def get_total_connections(self, head_index):
        return -1 * self.pointers[head_index]

    # returns the # of sets in the collection
    def get_set_count(self):
        return self.set_count

    # returns the index (in self.pointers) of the representative element of the disjoint
    # set that the given protein_name is in
    def find_set(self, protein_name):
        index = self.converter.get(protein_name)
        return self.find_helper(index)

    def find_helper(self, index):
        if (self.pointers[index] < 0):
            return index
        else:
            representative_index = self.find_helper(self.pointers[index])
            self.pointers[index] = representative_index
            return representative_index

    # returns true if any disjoint set in the collection contains the given item
    def contains(self, item):
        return item in self.converter.keys()

    # combines the two sets of item1 and item2; item1 becomes the new head (arbitrary in this case)
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

    # used for the case where we want to form a connection between two genes that are already indirectly
    # connected. This will not create an orphan as long as that set has some remaining unassigned connections
    def decrement_value(self, item):
        # reduces total # of connections remaining in this set (uses +1 vs -1 as sentinal values are negative)
        head = self.find_set(item)
        self.pointers[head] += 1

    # toString mostly used for debugging purposes
    def __repr__(self):
        return str(self.pointers) + str(self.converter)



get_model(10)