import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
import os, sys
import itertools as it
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Biotapestry import convert_biotapestry_to_antimony
from change_biotapestry import add_biotapestry

import scipy

"""
Runs differential evolution in an attempt to roughly estimate parameter values using data
Uses "broken" model, so the estimate is based on an incorrect model, but filtering out
liking poor/incorrect species from this estimation can improve its performance. You can
select which species to use in the estimation by altering the "selections" variable

"""
selections = ["mRNA" + str(i) for i in [1,2,3,6]]

# from first RNA seqeunce test
#rna_data = pd.read_csv("Orders/rnaseq_Noisy_Result.csv")
#rna_data.set_index('time', inplace=True) # need to wait for Yoshi to update .csv to include time
#rna_data.iloc[::30].plot(style='.-')
#plt.yscale('log')


# Load in our current "best guess" for the model
ant_str = convert_biotapestry_to_antimony("8gene_broken.csv", 8, [0.05, 2, 0.05, 2, 0.1, 3, 0.02])
r = te.loada(ant_str)

# fake data: remove when real data is obtained
data_str = convert_biotapestry_to_antimony("8gene_network.csv", 8,  [1.0/60, 1, 1.0/60, 1,5.0/60, 5, 1.0/60])
data_model = te.loada(data_str)
data = data_model.simulate(0,200,40)#, selections=selections)
#data_model.plot()

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
gene: source gene you want to look at
data: the results
timepoints: [start,stop, step] for r.simulate
"""
# TODO: return the top 10 best
# TODO: try/catch for add_biotap overflow
def estimate_connections(gene, data, timepoints, csv_filename, csv_newfile):  
    permConnection = []
    permError = float("inf")
    perms = list(it.permutations(gene))

    for ii in range(len(perms)):
        choice = perms[ii]
        numAdded = -1
        connection = [] # stores the chosen gene connection
        for jj in range(len(choice)):
            gene = choice[jj]
            print(connection)
            numAdded = -1
            singleConnection = connection[:]
#            if numAdded == 1: 
#                singleConnection.pop() # prepare singleConnection
#            elif numAdded == 2:
#                singleConnection.pop()
#                singleConnection.pop()
            singleError = float("inf")
            print(singleConnection)
            for i in range(0, 9): # 0 = flag for single connection
                for j in range(0, 9):
                    for k in (-1, 1):
                        for m in (-1, 1):
                            if (i != j):
                                if (i == 0):
                                    add = [(j, gene, k)]
                                elif (j == 0):
                                    add = [(i, gene, k)]
                                else:
                                    add = [(i, gene, k), (j, gene, m)]
                                # add the new connection
                                singleConnection.extend(add)
                                #print(add)
                                #print(singleConnection)
                                
                                add_biotapestry(singleConnection, csv_filename, csv_newfile)
                                ant_str = convert_biotapestry_to_antimony(csv_newfile, 8, 
                                                  [1/60, 1, 1/60, 1, 5/60, 5, 1/60])
                                # simulate
                                r = te.loada(ant_str)
                                start = timepoints[0]
                                stop = timepoints[1]
                                steps = timepoints[2]
                                result = r.simulate(start, stop, steps)
                                diff = data - result
                                error = np.sum(np.power(diff, 2))
                                # test error for best error
                                if error < singleError:
                                    # removes old (worst) connection
                                    if numAdded == 1:  
                                        connection.pop()
                                    elif numAdded == 2: 
                                        connection.pop()
                                        connection.pop()
                                    # stores how many new connections are added
                                    if len(add) == 1:
                                        numAdded = 1
                                    else:
                                        numAdded = 2
                                    singleError = error
                                    print("connection before add " + str(connection))
                                    connection.extend(add)
                                    #print(add)
                                    print("connection " + str(connection))
                                # remove choice for next loop
                                #print("before pop" + str(singleConnection))
                                singleConnection.pop()
                                if len(add) == 2:
                                    singleConnection.pop()   
                                #print(add)
                                #print("after " + str(singleConnection))
        print(choice)
        print(connection)                   
        print(str(singleError) + " " + str(permError))
        print()
        if singleError < permError:              
            permConnection = connection
            permError  = singleError
    return permConnection



# returns the total squared error between the data and the model output
# for parameter estimation
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

        exec("r.%s = %f" % (next_id, val))

    start = timepoints[0]
    stop = timepoints[1]
    steps = timepoints[2]
    result = r.simulate(start,stop,steps, selections=selections)
    diff = data - result
    return  np.sum(np.power(diff, 2))




#opt_sol = scipy.optimize.differential_evolution(lambda x: objective_func(x, r, data, [0,200,40], selections=selections), param_ranges, disp=True, popsize=30)
#print("\nParameter Estimation: [d_protein, d_mRNA, L, Vm, a_protein, H, K ] = " + str(opt_sol))

connection = estimate_connections([7,5], data, [0,200,40], "../Biotapestry/8gene_broken.csv", "../Biotapestry/8gene_ie.csv")
print("Best connection " + str(connection))