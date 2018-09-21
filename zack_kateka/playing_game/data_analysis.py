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
    mapping[str(error)] = [(-1,-1,-1)] # flag for add nothing; current model is already best
    
    
#    for i in range(1, len(gene) + 1):
#        perms.extend(list(it.permutations(gene, i)))
#    print(perms)
    #place = -1
    perms = list(it.permutations(gene))
    length = len(perms)
    count = 0
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
            for i in (0, 2, 4, 6, 8): # 0 = flag for single connection
                for j in (1, 3, 5, 7):
                    for k in (-1, 1):
                        for m in (-1, 1, 0):   
                            if m == 0: # don't add connection just add it to permError
                                permError.sort()
                                #place = -1
                                for kk in range(len(permError) - 1, -1, -1):
                                    if singleError <= permError[kk]:
                                        place = kk
                                if place != -1:
                                    permError.insert(kk, singleError)
                                    mapping[str(singleError)] = connection
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


# param set i.e. [Vm1, Vm2, Vm3, ...]
# ID = string representing single parameter of interest (i.e. "Vm")
def single_parameter_sweep(paramset, ID,  r, data, timepoints, selections):
    r.resetToOrigin()
    for i, next_param in enumerate(paramset):
        exec("r.%s = %f" % (ID + str(i+1), next_param)) #runs the given formatted string as if it is a line of code

    result = simulate(timepoints[0], timepoints[1], timepoints[2], selections=selections)
    diff = data - result
    return np.sum(np.power(diff,2))
