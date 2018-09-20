'''
@author: Zachary McNulty & Kateka Seth
'''


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
import itertools as it
import scipy

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Biotapestry import convert_biotapestry_to_antimony
from change_biotapestry import add_biotapestry


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
selections = ["mRNA" + str(i) for i in range(1,9)] 
timepoints = [0,200, 40]


'''
Load in experimental data.
'''
# from first RNA seqeunce test
data_table = pd.read_csv('model_files/RNASeq_HiRes.csv') # RNASeq_HiRes has timepoints = [0,200,20]
data_table.set_index('time', inplace=True) 
#data_table[selections].plot(style='.-')


# converts dataframe into numpy 2D array
#data = data_table.as_matrix(columns=selections)-----(depreciated)
data = data_table[selections].values

# Load in our current "best guess" for the model
ant_str = convert_biotapestry_to_antimony(broken_model, 8, [0.01556653, 9.959682  , 0.1056418 , 6.66957033, 0.08160472, 4.25284957, 0.06687737])
r = te.loada(ant_str)
r.simulate(timepoints[0],timepoints[1], timepoints[2], selections = ['time'] + selections) 
#r.plot()
#plt.show()

"""
param_ranges is a list of rough bounds for each parameter value in the form (min, max). For example,
the first entry (0,10) says we think d_protein is somewhere between 0 and 10.
"""
# smaller range
#param_ranges = [(0,0.25), (0,4), (0,0.25), (0,4), (0,0.5), (0,6), (0,0.25)]
# larger range
param_ranges = [(0,10), (0,10), (0,10), (0,10), (0,10), (0,10), (0,10)]



#testing purposes
data_str = convert_biotapestry_to_antimony("../Biotapestry/8gene_network.csv", 8,  [1.0/60, 1, 1.0/60, 1,5.0/60, 5, 1.0/60])
data_model = te.loada(data_str)	
data = data_model.simulate(0,200,40, selections=selections)



"""
the one for gene interaction
gene: source gene you want to look at
data: the experimental data
timepoints: [start,stop, step] for r.simulate
DONT RUN 5 
"""
# TODO: plz optimize
def estimate_connections(gene, data, timepoints, csv_filename, csv_newfile, selections):  
    mapping = {}
    permConnections = []
    permError = [float("inf")]*10
    perms = []
    
    ant_str = convert_biotapestry_to_antimony(csv_filename, 8, 
                      [1/60, 1, 1/60, 1, 5/60, 5, 1/60])
    # simulate
    r = te.loada(ant_str)
    start = timepoints[0]
    stop = timepoints[1]
    steps = timepoints[2]
    result = r.simulate(start, stop, steps, selections=selections)
    diff = data - result
    error = np.sum(np.power(diff, 2))
    
    permError[0] = error
    mapping[str(error)] = [(-1,-1,-1)]
    
    
    for i in range(1, len(gene) + 1):
        perms.extend(list(it.permutations(gene, i)))
    length = len(perms)
    count =0
    #print(perms)
    for ii in range(len(perms)):
        choice = perms[ii]
        numAdded = -1
        connection = [] # stores the chosen gene connection
        for jj in range(len(choice)):
            gene = choice[jj]
            #print(connection)
            numAdded = -1
            singleConnection = connection[:]
            singleError = float("inf")
            #print(singleConnection)
            for i in range(0, 9): # 0 = flag for single connection
                for j in range(1, 9):
                    for k in (-1, 1):
                        for m in (-1, 1):   
                            if (i != j):
                                if (i == 0):
                                    add = [(j, gene, k)]
                                else:
                                    add = [(i, gene, k), (j, gene, m)]
                                # add the new connection
                                singleConnection.extend(add)
                                #print(add)
                                #print(singleConnection)
                                try:
                                    add_biotapestry(singleConnection, csv_filename, csv_newfile)
                                except ValueError:
                                    break
                                ant_str = convert_biotapestry_to_antimony(csv_newfile, 8, 
                                                  [1/60, 1, 1/60, 1, 5/60, 5, 1/60])
                                os.remove(csv_newfile)
                                # simulate
                                r = te.loada(ant_str)
                                result = r.simulate(timepoints[0], timepoints[1], 
                                                    timepoints[2], selections=selections)
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
                                    #print("connection before add " + str(connection))
                                    connection.extend(add)
                                    #print(add)
                                    #print("connection " + str(connection))
                                # remove choice for next loop
                                #print("before pop" + str(singleConnection))
                                singleConnection.pop()
                                if len(add) == 2:
                                    singleConnection.pop()   
                                #print(add)
                                #print("after " + str(singleConnection))
#        print(choice)
#        print(connection)                   
#        print(str(singleError) + " " + str(permError))
#        print()
        count += 1
        print("progress: " + str(count) + "/" + str(length))
        permError.sort()
        place = -1
        for kk in range(len(permError) - 1, -1, -1):
            if singleError <= permError[kk]:
                place = kk
        if place != -1:
            permError.insert(kk, singleError)
            mapping[str(singleError)] = connection
            
    permError.sort()
    print(permError)
    for i in range(len(permError)):
        if permError[i] != float("inf"):
            permConnections.append(mapping.get(str(permError[i])))
    return permConnections




# returns the total squared error between the data and the simulated data generated from
# "best guess" model. Used for parameter estimation with an optimization tool requiring an
# objective function
#
# paramset: vector of values for each parameter (length of 7)
# r : roadrunner instance storing model
# data: numpy 2D array storing concentrations of each species at each timepoint
# timepoints: [start, stop, steps] gives information on timepoints that our model (r) 
#             will simulate species concentrations at; needs to match timepoints provided in data
# selections: which species should be compared; note, you will have to pre-process "data"
#             so that it only contains the species you are interested in

def objective_func(paramset, r, data, timepoints, selections):
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
Probes for possible connections; we can investigate the feasibility of these connections using further experimental data
'''
connection = estimate_connections([2,8,7,5], data, timepoints, "../Biotapestry/8gene_broken.csv", "../Biotapestry/8gene_ie.csv", selections)
print("Best connection " + str(connection))
