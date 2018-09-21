
'''
@author: Zachary McNulty & Kateka Seth
'''


import pandas as pd
import tellurium as te
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Biotapestry import convert_biotapestry_to_antimony
from data_analysis import estimate_connections
from data_analysis import objective_func
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
temp = "model_files/temp.csv"
selections = ["mRNA" + str(i) for i in [8]] 
selections2 = ["P" + str(i) for i in (2,8)]
selections3 = ["P" + str(i) for i in [1,2,3,4,5,6,8]]
timepoints = [0, 200, 10]

'''
Load in experimental data.
'''
## from first RNAseq test
#data_table = pd.read_csv('model_files/RNASeq_HiRes.csv') # RNASeq_HiRes has timepoints = [0,200,20]
#data_table.set_index('time', inplace=True) 
#data_table[selections].plot(style='.-')

# upreg 7, monitor 2 and 8
data_table2 = pd.read_csv('model_files/up7_fl2_fl8.csv') # RNASeq_HiRes has timepoints = [0,200,20]
data_table2.set_index('time', inplace=True) 
data_table2[selections2].plot(style='.-')

# mass spec test
#data_table3 = pd.read_csv('model_files/massSpec_lowRes.csv')
#data_table3.set_index('time', inplace=True) 
#data_table3[selections3].plot(style='.-')

# converts dataframe into numpy 2D array
#data = data_table.as_matrix(columns=selections)-----(depreciated)
data = data_table2[selections2].values


init_params =  [0.01556653, 9.959682  , 0.1056418 , 6.66957033, 0.08160472, 4.25284957, 0.06687737]

add = [(7, 2, -1)]
add_biotapestry(add, broken_model, temp)
ant_str = convert_biotapestry_to_antimony(temp, 8, init_params)
print(ant_str)
r = te.loada(
"""
model *pathway()

        // Compartments and Species:
        species INPUT, P1, mRNA1, P2, mRNA2, P3, mRNA3, P4, mRNA4, P5, mRNA5, P6, mRNA6, P7, mRNA7, P8, mRNA8;

        // Assignment Rules (production rates used in reactions):
        // transcription1(SR : in connections = [SA+SR (P4): ]) uses production rate := Vm1*( 1 /(1 +K1_1*P4^H1));
        // transcription2(SR : in connections = [N/A (P7): ]) uses production rate := Vm2*( 1 /(1 +K1_2*P7^H2));
        // transcription3(DR : in connections = [SA (P5): , SA (P6): ]) uses production rate := Vm3*( 1 /(1 + K1_3*P5^H3 + K2_3*P6^H3 + K1_3*K3_3*P5^H3*P6^H3));
        // transcription4(SA+SR : in connections = [SR (P1): , N/A (P8): ]) uses production rate := Vm4*((K1_4*P1^H4)/(1 + K1_4*P1^H4 + K2_4*P8^H4 + K1_4*K2_4*P1^H4*P8^H4));
        // transcription5(SA : in connections = [N/A (P8): ]) uses production rate := Vm5*((K1_5*P8^H5)/(1 +K1_5*P8^H5));
        // transcription6(SA : in connections = [None (INPUT): ]) uses production rate := Vm6*((K1_6*INPUT^H6)/(1 +K1_6*INPUT^H6));
        // transcription7(N/A : in connections = []) uses production rate := Vm7*(/);
        // transcription8(N/A : in connections = []) uses production rate := Vm8*(/);

        const d_protein1, d_mRNA1, L1, Vm1, a_protein1, H1, K1_1, K2_1, K3_1, d_protein2, d_mRNA2, L2, Vm2, a_protein2, H2, K1_2, K2_2, K3_2, d_protein3, d_mRNA3, L3, Vm3, a_protein3, H3, K1_3, K2_3, K3_3, d_protein4, d_mRNA4, L4, Vm4, a_protein4, H4, K1_4, K2_4, K3_4, d_protein5, d_mRNA5, L5, Vm5, a_protein5, H5, K1_5, K2_5, K3_5, d_protein6, d_mRNA6, L6, Vm6, a_protein6, H6, K1_6, K2_6, K3_6, d_protein7, d_mRNA7, L7, Vm7, a_protein7, H7, K1_7, K2_7, K3_7, d_protein8, d_mRNA8, L8, Vm8, a_protein8, H8, K1_8, K2_8, K3_8;

        // Reactions:
        //transcription1
        J1: => mRNA1 ; L1 + Vm1*( 1 /(1 +K1_1*P4^H1)) - d_mRNA1 * mRNA1;
        //translation1
        F1: => P1 ; a_protein1 * mRNA1 - d_protein1 * P1;
        //transcription2
        J2: => mRNA2 ; L2 + Vm2*( 1 /(1 +K1_2*P7^H2)) - d_mRNA2 * mRNA2;
        //translation2
        F2: => P2 ; a_protein2 * mRNA2 - d_protein2 * P2;
        //transcription3
        J3: => mRNA3 ; L3 + Vm3*( 1 /(1 + K1_3*P5^H3 + K2_3*P6^H3 + K1_3*K3_3*P5^H3*P6^H3)) - d_mRNA3 * mRNA3;
        //translation3
        F3: => P3 ; a_protein3 * mRNA3 - d_protein3 * P3;
        //transcription4
        J4: => mRNA4 ; L4 + Vm4*((K1_4*P1^H4)/(1 + K1_4*P1^H4 + K2_4*P8^H4 + K1_4*K2_4*P1^H4*P8^H4)) - d_mRNA4 * mRNA4;
        //translation4
        F4: => P4 ; a_protein4 * mRNA4 - d_protein4 * P4;
        //transcription5
        J5: => mRNA5 ; L5 + Vm5*((K1_5*P8^H5)/(1 +K1_5*P8^H5)) - d_mRNA5 * mRNA5;
        //translation5
        F5: => P5 ; a_protein5 * mRNA5 - d_protein5 * P5;
        //transcription6
        J6: => mRNA6 ; L6 + Vm6*((K1_6*INPUT^H6)/(1 +K1_6*INPUT^H6)) - d_mRNA6 * mRNA6;
        //translation6
        F6: => P6 ; a_protein6 * mRNA6 - d_protein6 * P6;
        //transcription7
        J7: => mRNA7 ; L7 + 0 - d_mRNA7 * mRNA7;
        //translation7
        F7: => P7 ; a_protein7 * mRNA7 - d_protein7 * P7;
        //transcription8
        J8: => mRNA8 ; L8 + 0 - d_mRNA8 * mRNA8;
        //translation8
        F8: => P8 ; a_protein8 * mRNA8 - d_protein8 * P8;

        // Species initializations:
        INPUT = 1;
        mRNA1 = 0;
        P1 = 0;
        mRNA2 = 0;
        P2 = 0;
        mRNA3 = 0;
        P3 = 0;
        mRNA4 = 0;
        P4 = 0;
        mRNA5 = 0;
        P5 = 0;
        mRNA6 = 0;
        P6 = 0;
        mRNA7 = 0;
        P7 = 0;
        mRNA8 = 0;
        P8 = 0;
        d_protein1 = 0.01556653;
        d_mRNA1 = 9.959682;
        L1 = 0.1056418;
        Vm1 = 6.66957033;
        a_protein1 = 0.08160472;
        H1 = 5.0;
        K1_1 = 0.06687737;
        K2_1 = 0.06687737;
        K3_1 = 0.06687737;
        d_protein2 = 0.01556653;
        d_mRNA2 = 9.959682;
        L2 = 0.1056418;
        Vm2 = 6.66957033;
        a_protein2 = 0.08160472;
        H2 = 5.0;
        K1_2 = 0.06687737;
        K2_2 = 0.06687737;
        K3_2 = 0.06687737;
        d_protein3 = 0.01556653;
        d_mRNA3 = 9.959682;
        L3 = 0.1056418;
        Vm3 = 6.66957033;
        a_protein3 = 0.08160472;
        H3 = 5.0;
        K1_3 = 0.06687737;
        K2_3 = 0.06687737;
        K3_3 = 0.06687737;
        d_protein4 = 0.01556653;
        d_mRNA4 = 9.959682;
        L4 = 0.1056418;
        Vm4 = 6.66957033;
        a_protein4 = 0.08160472;
        H4 = 5.0;
        K1_4 = 0.06687737;
        K2_4 = 0.06687737;
        K3_4 = 0.06687737;
        d_protein5 = 0.01556653;
        d_mRNA5 = 9.959682;
        L5 = 0.1056418;
        Vm5 = 6.66957033;
        a_protein5 = 0.08160472;
        H5 = 5.0;
        K1_5 = 0.06687737;
        K2_5 = 0.06687737;
        K3_5 = 0.06687737;
        d_protein6 = 0.01556653;
        d_mRNA6 = 9.959682;
        L6 = 0.1056418;
        Vm6 = 6.66957033;
        a_protein6 = 0.08160472;
        H6 = 5.0;
        K1_6 = 0.06687737;
        K2_6 = 0.06687737;
        K3_6 = 0.06687737;
        d_protein7 = 0.01556653;
        d_mRNA7 = 9.959682;
        L7 = 0.1056418;
        Vm7 = 6.66957033;
        a_protein7 = 0.08160472;
        H7 = 5.0;
        K1_7 = 0.06687737;
        K2_7 = 0.06687737;
        K3_7 = 0.06687737;
        d_protein8 = 0.01556653;
        d_mRNA8 = 9.959682;
        L8 = 0.1056418;
        Vm8 = 6.66957033;
        a_protein8 = 0.08160472;
        H8 = 5.0;
        K1_8 = 0.06687737;
        K2_8 = 0.06687737;
        K3_8 = 0.06687737;


end

"""
)
r.Vm2 = 12.5
r.Vm7 = 10
r.Vm1 = 8
r.simulate(timepoints[0],timepoints[1], timepoints[2], selections = ['time'] + selections3) 
r.plot()
# Load in our current "best guess" for the model
#for i in range(1,9):
#    for j in range(1,9):
#        for k in (-1, 1):
#            for m in (-1, 1):
#                add = [(6,8,1),(8,8,1),(7,2,-1),(i,5,k),(j,5,m)]
#                print(add)
#                add_biotapestry(add, broken_model, temp)
#                ant_str = convert_biotapestry_to_antimony(temp, 8, init_params)
#                r = te.loada(ant_str)
#                r.Vm2 = 12.5
#                r.Vm7 = 10
#                r.Vm1 = 8
#                r.simulate(timepoints[0],timepoints[1], timepoints[2], selections = ['time'] + selections3) 
#                r.plot()
#for i in range(1,9):
#    for k in (-1, 1):
#        
#        add = [(5, 7, -1), (7, 7, -1), (8, 2, -1), (2, 8, -1), (3, 8, 1),(7, 6, -1),(i,1,k)]
#        print(add)
#        add_biotapestry(add, broken_model, temp)
#        ant_str = convert_biotapestry_to_antimony(temp, 8, init_params)
#        r = te.loada(ant_str)
#        r.Vm2 = 12.5
#        r.Vm7 = 10
#        r.Vm1 = 8
#        r.simulate(timepoints[0],timepoints[1], timepoints[2], selections = ['time'] + selections3) 
#        r.plot()
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
# def estimate_connections(gene, data, timepoints, csv_filename, csv_newfile, selections, params):  
#connection = estimate_connections([2,7,8], 8, data, timepoints, broken_model, "model_files/temp2.csv", selections3, init_params)
#print("Best connection " + str(connection))
