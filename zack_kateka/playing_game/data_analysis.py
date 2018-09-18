import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
import os, sys
import itertools as itl
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Biotapestry import convert_biotapestry_to_antimony

# from first RNA seqeunce test
rna_data = pd.read_csv("Orders/rnaseq_Noisy_Result.csv")
ant_str = convert_biotapestry_to_antimony("8gene_network.csv", 8,  [1.0/60,1,1.0/60,1,5.0/60,5,1.0/60])
rna_data.set_index('time', inplace=True) # need to wait for Yoshi to update .csv to include time

print(rna_data.head())

var_names = list(rna_data)

rna_data.iloc[::30].plot(style='.-')
plt.yscale('log')

selections =['time'] + ["mRNA" + str(i+1) for i in range(8)]
r = te.loada(ant_str)
r.simulate(0,300,1000, selections=selections) #TO DO: update so time scale matches data
r.plot()



#plt.show()


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

    ids = r.getGlobalParameterIds()
    paramsets = next_paramset(param_possibilities)
    for paramset  in paramsets:
        for i, next_id in enumerate(ids):
            if i % 9 < 6:
                val = paramset[i%9]
            else:
                val = paramset[6]
            exec("r.%s = %d" % (next_id, val))

        r.resetToOrigin()
        result = r.simulate(start,stop,steps, selections=selections)

param_ranges = [(0,5), (0,5), (0,5), (0,60), (0,60), (0,5), (0,5)]
step_count = 2
param_possibilities = [np.linspace(x[0],x[1],step_count) for x in param_ranges]
