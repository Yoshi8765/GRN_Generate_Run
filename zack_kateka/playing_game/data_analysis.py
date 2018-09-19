import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Biotapestry import convert_biotapestry_to_antimony
from change_biotapestry import add_biotapestry

import scipy

"""
Runs differential evolution in an attempt to roughly estimate parameter values using data
Uses "broken" model, so the estimate is based on an incorrect model, but filtering out
liking poor/incorrect species from this estimation can improve its performance. You can
select which species to use in the estimation by altering the "selections" variable
Timepoints is the a vector describing the timepoints present in the data. It is in the form
[start, stop, step_size]. This is needed so the simulated data from the partially complete model
lines up with the experimental data
"""
broken_model = "model_files/biotapestry_broken.csv"
selections = ["mRNA" + str(i) for i in range(1,9)] 
timepoints = [0,200, 20]


'''
Load in experimental data.
'''
# from first RNA seqeunce test
data_table = pd.read_csv('model_files/RNASeq_HiRes.csv')
data_table.set_index('time', inplace=True) # need to wait for Yoshi to update .csv to include time
data_table[selections].plot(style='.-')


# converts dataframe into numpy 2D array
#data = data_table.as_matrix(columns=selections)-----(depreciated)
data = data_table[selections].values

# Load in our current "best guess" for the model
ant_str = convert_biotapestry_to_antimony(broken_model, 8, [0.01556653, 9.959682  , 0.1056418 , 6.66957033, 0.08160472, 4.25284957, 0.06687737])
r = te.loada(ant_str)
r.simulate(timepoints[0],timepoints[1], timepoints[2], selections = ['time'] + selections) 
r.plot()
#plt.show()

"""
param_ranges is a list of rough bounds for each parameter value in the form (min, max). For example,
the first entry (0,10) says we think d_protein is somewhere between 0 and 10.
"""

# smaller range
#param_ranges = [(0,0.25), (0,4), (0,0.25), (0,4), (0,0.5), (0,6), (0,0.25)]
# larger range
param_ranges = [(0,10), (0,10), (0,10), (0,10), (0,10), (0,10), (0,10)]


# gene interaction


"""
the one for gene interaction
gene: the source gene 
data: the results
timepoints: [start,stop, step] for r.simulate
"""
def interaction_estimate(gene, data, timepoints, csv_filename, csv_newfile, selections=None):  
    bestError = 0xFFFFFFFF
    connection = []
    # 0 for single source
    for i in range(0, 9):
        for j in range(0, 9):
            for k in (-1, 1):
                for m in (-1, 1):
                    if (i != j):
                        if (i != 0 and j != 0):
                            if (i == 0):
                                add = [(j, gene, k)]
                            elif (j == 0):
                                add = [(i, gene, k)]
                            else:
                                add = [(i, gene, k), (j, gene, m)]
                            add_biotapestry(add, csv_filename, csv_newfile)
                            ant_str = convert_biotapestry_to_antimony(csv_newfile, 8, 
                                                                       [1/60, 1, 1/60, 1, 5/60, 5, 1/60])
                            r = te.loada(ant_str)
                            if selections == None:
                                selections = r.getFloatingSpeciesIds() + r.getBoundarySpeciesIds()
                            start = timepoints[0]
                            stop = timepoints[1]
                            steps = timepoints[2]
                            result = r.simulate(start,stop,steps, selections=selections)
                            diff = data - result
                            error = np.sum(np.power(diff, 2))
                            if error < bestError:
                                bestError = error
                                connection = add
    return connection
    
    



# returns the total squared error between the data and the simulated data generated from
# "best guess" model. Used for parameter estimation
#
# paramset: vector of values for each parameter (length of 7)
# r : roadrunner instance storing model
# data: numpy 2D array storing concentrations of each species at each timepoint
# timepoints: [start, stop, steps] gives information on timepoints that our model (r) 
#             will simulate species concentrations at; needs to match timepoints provided in data
# selections: which species should be compared; note, you will have to pre-process "data"
#             so that it only contains the species you are interested in

def objective_func(paramset, r, data, timepoints, selections=None):
    if selections == None:
        selections = ['time'] + r.getFloatingSpeciesIds() + r.getBoundarySpeciesIds()
    ids = r.getGlobalParameterIds()
    r.resetToOrigin()
    for i, next_id in enumerate(ids):
        if i % 9 < 6:
            val = paramset[i%9]
        else:
            val = paramset[6]

        exec("r.%s = %f" % (next_id, val)) #runs the given formatted string as if it is a line of code

    start = timepoints[0]
    stop = timepoints[1]
    steps = timepoints[2]
    result = r.simulate(start,stop,steps, selections=selections)
    diff = data - result
    return  np.sum(np.power(diff, 2))



'''
Run objective_func through differential evolution to estimate parameters ['d_protein', 'd_mRNA', 'L', 'Vm', 'a_protein', 'H', 'K']
'''
# other parameters for optimize.differential_evolution
#   popsize : increasing this will increase search radius; may lead to better solution but slows algorithm
#   mutation : scales mutation phase. Larger numbers increase search radius (improves solution), but slows convergence
#   recombination : higher numbers increase randomness. May lead to better solutions, but can increase instability

#opt_sol = scipy.optimize.differential_evolution(lambda x: objective_func(x, r, data, timepoints, selections=selections), param_ranges, disp=True, popsize=40, mutation = (1,1.9))
#print("\nParameter Estimation: [d_protein, d_mRNA, L, Vm, a_protein, H, K ] = " + str(opt_sol))

'''
Tries to find possible gene connections. To test is a list of genes you want to find the best in-connections to, given the current model. Does NOT 
overwrite broken model, so can only test one gene at a time. Multiple iterations would be required to attempt to repair network
'''
to_test = [2,7]
for i in to_test:
    connection = interaction_estimate(i, data, timepoints, broken_model, "model_files/8gene_ie.csv", selections=selections)
    print("Best connection for " + str(i) + str(connection))
