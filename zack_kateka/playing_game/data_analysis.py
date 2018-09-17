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
ant_str = """
model *pathway()

        // Compartments and Species:
        species INPUT, P1, mRNA1, P2, mRNA2, P3, mRNA3, P4, mRNA4, P5, mRNA5, P6, mRNA6, P7, mRNA7, P8, mRNA8;

        // Assignment Rules (production rates used in reactions):
        // transcription1(DA : in connections = [P3, P4]) uses production rate := Vm1*((K1_1*P3^H1 + K2_1*P4^H1 + K1_1*K3_1*P3^H1*P4^H1)/(1 + K1_1*P3^H1 + K2_1*P4^H1 + K1_1*K3_1*P3^H1*P4^H1));
        // transcription2(SR : in connections = [P8]) uses production rate := Vm2*( 1 /(1 +K1_2*P3^H2));
        // transcription3(SR : in connections = [INPUT]) uses production rate := Vm3*( 1 /(1 +K1_3*P3^H3));
        // transcription4(N/A : in connections = []) uses production rate := Vm4*(/);
        // transcription5(N/A : in connections = []) uses production rate := Vm5*(/);
        // transcription6(DA : in connections = [P3, P6]) uses production rate := Vm6*((K1_6*P3^H6 + K2_6*P6^H6 + K1_6*K3_6*P3^H6*P6^H6)/(1 + K1_6*P3^H6 + K2_6*P6^H6 + K1_6*K3_6*P3^H6*P6^H6));
        // transcription7(N/A : in connections = []) uses production rate := Vm7*(/);
        // transcription8(N/A : in connections = []) uses production rate := Vm8*(/);

        // Reactions:
        //transcription1
        J1: => mRNA1 ; L1 + Vm1*((K1_1*P3^H1 + K2_1*P4^H1 + K1_1*K3_1*P3^H1*P4^H1)/(1 + K1_1*P3^H1 + K2_1*P4^H1 + K1_1*K3_1*P3^H1*P4^H1)) - d_mRNA1 * mRNA1;
        //translation1
        F1: => P1 ; a_protein1 * mRNA1 - d_protein1 * P1;
        //transcription2
        J2: => mRNA2 ; L2 + Vm2*( 1 /(1 +K1_2*P3^H2)) - d_mRNA2 * mRNA2;
        //translation2
        F2: => P2 ; a_protein2 * mRNA2 - d_protein2 * P2;
        //transcription3
        J3: => mRNA3 ; L3 + Vm3*( 1 /(1 +K1_3*P3^H3)) - d_mRNA3 * mRNA3;
        //translation3
        F3: => P3 ; a_protein3 * mRNA3 - d_protein3 * P3;
        //transcription4
        J4: => mRNA4 ; L4 + 0 - d_mRNA4 * mRNA4;
        //translation4
        F4: => P4 ; a_protein4 * mRNA4 - d_protein4 * P4;
        //transcription5
        J5: => mRNA5 ; L5 + 0 - d_mRNA5 * mRNA5;
        //translation5
        F5: => P5 ; a_protein5 * mRNA5 - d_protein5 * P5;
        //transcription6
        J6: => mRNA6 ; L6 + Vm6*((K1_6*P3^H6 + K2_6*P6^H6 + K1_6*K3_6*P3^H6*P6^H6)/(1 + K1_6*P3^H6 + K2_6*P6^H6 + K1_6*K3_6*P3^H6*P6^H6)) - d_mRNA6 * mRNA6;
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

        // Variable initializations:
        d_protein1 = 0.5;
        d_mRNA1 = 0.9;
        L1 = 0.8;
        Vm1 = 30.0;
        a_protein1 = 30.0;
        H1 = 1.0;
        K1_1 = 0.5;
        K2_1 = 0.5;
        K3_1 = 0.5;
        d_protein2 = 0.5;
        d_mRNA2 = 0.9;
        L2 = 0.8;
        Vm2 = 30.0;
        a_protein2 = 30.0;
        H2 = 1.0;
        K1_2 = 0.5;
        K2_2 = 0.5;
        K3_2 = 0.5;
        d_protein3 = 0.5;
        d_mRNA3 = 0.9;
        L3 = 0.8;
        Vm3 = 30.0;
        a_protein3 = 30.0;
        H3 = 1.0;
        K1_3 = 0.5;
        K2_3 = 0.5;
        K3_3 = 0.5;
        d_protein4 = 0.5;
        d_mRNA4 = 0.9;
        L4 = 0.8;
        Vm4 = 30.0;
        a_protein4 = 30.0;
        H4 = 1.0;
        K1_4 = 0.5;
        K2_4 = 0.5;
        K3_4 = 0.5;
        d_protein5 = 0.5;
        d_mRNA5 = 0.9;
        L5 = 0.8;
        Vm5 = 30.0;
        a_protein5 = 30.0;
        H5 = 1.0;
        K1_5 = 0.5;
        K2_5 = 0.5;
        K3_5 = 0.5;
        d_protein6 = 0.5;
        d_mRNA6 = 0.9;
        L6 = 0.8;
        Vm6 = 30.0;
        a_protein6 = 30.0;
        H6 = 1.0;
        K1_6 = 0.5;
        K2_6 = 0.5;
        K3_6 = 0.5;
        d_protein7 = 0.5;
        d_mRNA7 = 0.9;
        L7 = 0.8;
        Vm7 = 30.0;
        a_protein7 = 30.0;
        H7 = 1.0;
        K1_7 = 0.5;
        K2_7 = 0.5;
        K3_7 = 0.5;
        d_protein8 = 0.5;
        d_mRNA8 = 0.9;
        L8 = 0.8;
        Vm8 = 30.0;
        a_protein8 = 30.0;
        H8 = 1.0;
        K1_8 = 0.5;
        K2_8 = 0.5;
        K3_8 = 0.5;

        // Other declarations:
        const d_protein1, d_mRNA1, L1, Vm1, a_protein1, H1, K1_1, K2_1, K3_1, d_protein2, d_mRNA2, L2, Vm2, a_protein2, H2, K1_2, K2_2, K3_2, d_protein3, d_mRNA3, L3, Vm3, a_protein3, H3, K1_3, K2_3, K3_3, d_protein4, d_mRNA4, L4, Vm4, a_protein4, H4, K1_4, K2_4, K3_4, d_protein5, d_mRNA5, L5, Vm5, a_protein5, H5, K1_5, K2_5, K3_5, d_protein6, d_mRNA6, L6, Vm6, a_protein6, H6, K1_6, K2_6, K3_6, d_protein7, d_mRNA7, L7, Vm7, a_protein7, H7, K1_7, K2_7, K3_7, d_protein8, d_mRNA8, L8, Vm8, a_protein8, H8, K1_8, K2_8, K3_8;
        end
"""
#convert_biotapestry_to_antimony("8gene_broken.csv",8, [1]*7)
#rna_data.set_index("time") # need to wait for Yoshi to update .csv to include time


var_names = list(rna_data)

rna_data.iloc[::30].plot(style='.-')
plt.yscale('log')


selections =['time'] + ["mRNA" + str(i+1) for i in range(8)]
r = te.loada(ant_str)
r.simulate(0,50,100, selections=selections) #TO DO: update so time scale matches data
r.plot()

plt.show()


# Running a parameter scan to find minimum ideally
# parameters = init_params = ['d_proteinX', 'd_mRNAX' , 'LX' , 'VmX' , 'a_proteinX' , 'HX', 'K_X']

# (min, max) for each parameter value (choose reasonable values that we would expect from biology of system)
#param_ranges = [(0,5), (0,5), (0,5), (0,100), (0,100), (0,10), (0,5), (0,5),(0,5)]
#step_count = 10

#param_possibilities = [np.linspace(x[0],x[1],step_count) for x in param_ranges]
#x = param_possibilities
#choices = [(a,b,c,d,e,f,g,h,i) for a in x[0] for b in x[1] for c in x[2] for d in x[3] for e in x[4] for f in x[5] for g in x[6] for h in x[7] for i in x[8]]
#for choice in choices:
#    print (choice)
