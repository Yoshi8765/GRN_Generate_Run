
import tellurium as te
from GetModel import convert_to_antimony
from GetModel import Gene

# From Biotapestry: Export > Export to SBML
def convert_biotapestry_to_antimony(csv_filename, num_genes, init_param):
    f = open(csv_filename)
    all_genes = []
    
    for line in f:
        words = line.split(",")
        if words[0] == "general":


print(convert_biotapestry_to_antimony("pathway_biotapestry.csv",8,init_param = [1,1,1,1,1,1,1]))
