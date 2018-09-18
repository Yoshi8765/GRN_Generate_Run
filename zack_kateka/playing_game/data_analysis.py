import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
import os, sys
import itertools as itl
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Biotapestry import convert_biotapestry_to_antimony

import scipy

# from first RNA seqeunce test
#rna_data = pd.read_csv("Orders/rnaseq_Noisy_Result.csv")
data_str = convert_biotapestry_to_antimony("8gene_network.csv", 8,  [1.0/60, 1, 1.0/60, 1,5.0/60, 5, 1.0/60])
#rna_data.set_index('time', inplace=True) # need to wait for Yoshi to update .csv to include time
ant_str = convert_biotapestry_to_antimony("8gene_broken.csv", 8, [0.05, 2, 0.05, 2, 0.1, 3, 0.02])

#rna_data.iloc[::30].plot(style='.-')
#plt.yscale('log')

selections =['time'] + ["mRNA" + str(i+1) for i in range(8)]
r = te.loada(ant_str)
#r.simulate(0,300,1000, selections=selections) #TO DO: update so time scale matches data
#r.plot()

# fake data: remove when real data is obtained
data_model = te.loada(data_str)
data = data_model.simulate(0,200,40, selections=selections)
#data_model.plot()

param_ranges = [(0,0.25), (0,4), (0,0.25), (0,4), (0,0.5), (0,6), (0,0.25)]
step_count = 6 
param_possibilities = [np.linspace(x[0],x[1],step_count) for x in param_ranges]

# Running a parameter scan to find minimum ideally
# parameters = init_params = ['d_proteinX', 'd_mRNAX' , 'LX' , 'VmX' , 'a_proteinX' , 'HX', 'K_X']

def next_paramset(x):
    for a in x[0]:
        for b in x[1]:
            for c in x[2]:
                for d in x[3]:
                    for e in x[4]:
                        for f in x[5]:
                            for g in x[6]:
                                yield (a,b,c,d,e,f,g)


def objective_func(paramset, r, data, timepoints, selections=None):
    print (paramset)
    if selections == None:
        selections = ['time'] + r.getFloatingSpeciesIds() + r.getBoundarySpeciesIds()
    ids = r.getGlobalParameterIds()
    print(ids)
    r.resetAll()
    r.resetToOrigin()
    for i, next_id in enumerate(ids):
        if i % 9 < 6:
            val = paramset[i%9]
        else:
            val = paramset[6]

        exec("r.%s = %d" % (next_id, val))
        print(next_id)

    print(r.getGlobalParameterValues())

    start = timepoints[0]
    stop = timepoints[1]
    steps = timepoints[2]
    result = r.simulate(start,stop,steps, selections=selections)
    diff = data - result
    return  np.sum(np.power(diff, 2))

opt_sol = scipy.optimize.differential_evolution(lambda x: objective_func(x, r, data, [0,200,40], selections=selections), param_ranges, disp=True)
print(opt_sol)











def find_param_estimate(param_possibilities, r, data, timepoints, selections=None):
    if selections == None:
        selections = ['time'] + r.getFloatingSpeciesIds() + r.getBoundarySpeciesIds()
    ids = r.getGlobalParameterIds()
    all_paramsets = next_paramset(param_possibilities)
    
    best_guess = [1,1,1,1,1,1,1]
    guess_error = 2**32 

    for paramset in all_paramsets:
        r.resetToOrigin()
        for i, next_id in enumerate(ids):
            if i % 9 < 6:
                val = paramset[i%9]
            else:
                val = paramset[6]
            exec("r.%s = %d" % (next_id, val))
        
        start = timepoints[0]
        stop = timepoints[1]
        steps = timepoints[2]
        result = r.simulate(start,stop,steps, selections=selections)

        diff = data - result
        error = np.sum(np.power(diff, 2))
        if error < guess_error:
            best_guess = paramset
            guess_error = error
            print (best_guess)
            print (guess_error)

    return best_guess
    
#param_est  = find_param_estimate(param_possibilities, r, data, [0,200,40], selections=selections)
#print (param_est)
