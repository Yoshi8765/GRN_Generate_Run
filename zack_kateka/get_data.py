from run_experiments import export_experiments

"""
This file will send out the next batch of experimental data. 
Simply ensure that num_genes is the number of genes in the network, and run the script.
"""
num_genes = 8

export_experiments(num_genes, sendEmail=True, updateMoney=True)
