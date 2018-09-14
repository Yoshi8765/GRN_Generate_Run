
import tellurium as te
from GetModel import convert_to_antimony
from GetModel import Gene

# From Biotapestry: Export > Export to SBML
def convert_biotapestry_to_antimony(csv_filename, num_genes, init_param):
    f = open(csv_filename)
    all_connects = {}
    for i in range(num_genes):
        all_in_connects["P"+ str(i+1)] = []
   
   all_genes = []
    # name : [(P1, Activator), (P2, Repressor)] 
    for line in f:
        words = line.split(",")
        if words[0].rstrip() == "general":
            gene_source = words[3].rstrip()
            if gene_source == "INPUT":
                source_name = "INPUT"
            else:
                source__name = "P" + gene_source[5:].rstrip()

            gene_target = words[5].rstrip()
            target_name = "P" + gene_target[5:].rstrip()
            
            reg_type = words[7].rstrip()

            all_in_connects[target_name].append((source_name, reg_type))

    for protein_name in all_in_connects.keys():
        in_connects = all_in_connects[protein_name]
        if len(in_connects) == 1:
            if  in_connects[0][1] == "positive":
                next_gene = Gene(protein_name, "SA")
                Gene.add_in_connection
            else:
        else:

print(convert_biotapestry_to_antimony("pathway_biotapestry.csv",8,init_param = [1,1,1,1,1,1,1]))
