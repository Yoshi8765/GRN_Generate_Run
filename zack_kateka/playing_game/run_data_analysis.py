
'''
@author: Zachary McNulty & Kateka Seth
'''


import pandas as pd
import tellurium as te
import matplotlib.pyplot as plt
import scipy
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Biotapestry import convert_biotapestry_to_antimony
from data_analysis import estimate_connections
from data_analysis import objective_func
from data_analysis import single_parameter_sweep

"""
Runs differential evolution in an attempt to roughly estimate parameter values using data
Uses "broken" model, so the estimate is based on an incorrect model, but filtering out
likely poor/incorrect species from this estimation can improve its performance. You can
select which species to use in the estimation by altering the "selections" variable
Timepoints is the a vector describing the timepoints present in the data. It is in the form
[start, stop, step_size]. This is needed so the simulated data from the partially complete model
lines up with the experimental data
"""
broken_model = "model_files/biotapestry_broken.csv"
mRNA_selections = ["mRNA" + str(i) for i in [1,2,3,4,5,6,7,8]]
p_selections = ["P" + str(i) for i in [1,2,3,4,5,6,7,8]]
selections = p_selections + mRNA_selections
timepoints = [0,200, 11]
num_genes = 8
param_ranges = [(0,10)] * num_genes

'''
Load in experimental data.
'''
# from first RNAseq test (high res timepoints = [0,200,20])
rna_data_table = pd.read_csv('model_files/RNASeq_HiRes.csv') # RNASeq_HiRes has timepoints = [0,200,20]
rna_data_table.set_index('time', inplace=True) 
rna_data_table[mRNA_selections].plot(style='.-')
rna_data = rna_data_table[mRNA_selections].values # converts data to 2D numpy array

# from Mass Spec (low res timepoints = [0,200,10])
ms_data_table = pd.read_csv('model_files/massSpec_lowRes.csv')
ms_data_table.set_index('time', inplace=True)
ms_data_table[p_selections].plot(style='.-')
protein_data = ms_data_table[p_selections].values # converts data to 2D numpy array

# aggregate mrna and protein data
rna = rna_data_table.iloc[::2][mRNA_selections] 
aggregated = pd.concat([ ms_data_table[p_selections], rna], axis = 1)
aggregated_data = aggregated.values
aggregated.plot()
plt.show()

init_params =  [0.037, 9.959682  , 0.1056418 , 6.66957033,0.0816 , 4.25284957, 0.06687737]

# Load in our current "best guess" for the model
ant_str = convert_biotapestry_to_antimony(broken_model, 8,init_params)
r = te.loada(ant_str)
r.simulate(timepoints[0],timepoints[1], 1000, selections = ['time'] + p_selections) 
r.plot()

plt.show()

'''
Run objective_func through differential evolution to estimate parameters ['d_protein', 'd_mRNA', 'L', 'Vm', 'a_protein', 'H', 'K']
'''
# other parameters for optimize.differential_evolution
#   popsize : increasing this will increase search radius; may lead to better solution but slows algorithm
#   mutation : scales mutation phase. Larger numbers increase search radius (improves solution), but slows convergence
#   recombination : higher numbers increase randomness. May lead to better solutions, but can increase instability

#single_param_ranges = [(0,10)] * num_genes # 8 genes
#ID = "Vm"
#opt_sol = scipy.optimize.differential_evolution(lambda x: single_parameter_sweep(x,ID, r, aggregated_data, timepoints, selections=selections), single_param_ranges, disp=True, popsize=40, mutation = (1,1.9))

#opt_sol = scipy.optimize.differential_evolution(lambda x: objective_func(x, r, protein_data, timepoints, selections=p_selections), param_ranges, disp=True, popsize=40, mutation = (1,1.9))
#print("\nParameter Estimation: [d_protein, d_mRNA, L, Vm, a_protein, H, K ] = " + str(opt_sol))


'''
Probes for possible connections; we can investigate the feasibility of these connections using further experimental data
'''
# def estimate_connections(gene, data, timepoints, csv_filename, csv_newfile, selections, params):  
#connection = estimate_connections([2,7], 8, data, timepoints, broken_model, "model_files/temp.csv", selections, init_params)
#print("Best connection " + str(connection))
