# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 12:30:55 2018

@author: Kateka Seth
"""
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from change_biotapestry import remove_biotapestry
from change_biotapestry import add_biotapestry
from GetModel import get_model
from GetModel import convert_to_biotapestry
from Biotapestry import convert_biotapestry_to_antimony
#from data_analysis import interaction_estimate
import tellurium as te
import roadrunner
import antimony

# Round 2 of game (new network) analysis
csv_broken = "../playing_game/model_files/biotapestry_broken.csv"
csv_new = "../playing_game/model_files/test_added_connections.csv"

selections = ['time'] + ["mRNA" + str(i) for i in range(1,9)] + ['P2'] + ['P8']

#for i in range(1,9):
#    add=[(i,2,-1)]
#    print(i,j)
#    add_biotapestry(add,csv_broken, csv_new)
#    ant_str = convert_biotapestry_to_antimony(csv_new, 8, [0.01556653, 9.959682  , 0.1056418 , 6.66957033, 0.08160472,
#                                                                                                4.25284957, 0.06687737])
#                                                                 
#    r=te.loada(ant_str)
#    r.simulate(0,300,100,selections=selections)
#    r.plot(figsize=[7,7],xlim=(0,300),linewidth=2,linestyle='-') 
#    r.resetToOrigin()


# round 2 parameters: [0.40242652, 8.39507097, 0.0136378 , 9.99897755, 0.98195867,
#                      8.1766647 , 0.5442412]
add_biotapestry([(7,2,-1)], csv_broken, csv_new)
ant_str = convert_biotapestry_to_antimony(csv_new, 8, [0.01556653, 9.959682  , 0.1056418 , 6.66957033, 0.08160472,
                                                                                            4.25284957, 0.06687737])
print(ant_str)


r=te.loada(
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
        Vm1 = 8;
        a_protein1 = 0.08160472;
        H1 = 5.0;
        K1_1 = 0.06687737;
        K2_1 = 0.06687737;
        K3_1 = 0.06687737;
        d_protein2 = 0.01556653;
        d_mRNA2 = 9.959682;
        L2 = 0.1056418;
        Vm2 = 12.5;
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
        Vm7 = 10;
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
""")
    
r=te.loada(ant_str)
r.simulate(0,220,100,selections=selections)
r.plot(figsize=[7,7],xlim=(0,220),ylim=(0,8),linewidth=2,linestyle='-')



#----------------------------------------------------

#csv_filename="../Biotapestry/8gene_broken.csv"
#csv_newfile="../Biotapestry/8gene_test.csv"
#csv_remove="../Biotapestry/8gene_remove.csv"
#init_params = ['d_p', 'd_m' , 'L' , 'Vm' , 'a_p' , 'H', 'K']
#add_biotapestry([(6,8,1),(3,5,-1),(6,7,-1),(7,7,-1),(2,4,-1),(6,4,-1)],csv_filename, csv_newfile)

#remove_biotapestry([(7,7),(3,1),(1,8),(1,7)],"../Biotapestry/8gene_network.csv",csv_remove)

#actual parameters=[1/60, 1, 1/60, 1, 5/60, 5, 1/60]
#small range=[8.4111E-05, 5.919E-01, 1.88693E-01, 4.035919E-01,5.38276E-02, 5.92, 2.9811E-02]
#large range=[3.13755E-05, 9.5372, 3.21424, 6.66695, 3.6342586E-02, 9.9475, 5.24815884E-02]

#ant_str = convert_biotapestry_to_antimony("../Biotapestry/8gene_test.csv", 8, [2.48117E-06, 5.907979, 1/60, 4.08598971,
#                                                                               4.6555711E-02, 9.965566, 1.69489E-02])
#r=te.loada(ant_str)
#leaks=["time","P8","P5","P7","P4"]
#protein=["time", "INPUT","P1","P2","P3","P4","P5","P6","P7","P8"]
#mRNA=["time","mRNA1","mRNA2","mRNA3","mRNA4","mRNA5","mRNA6","mRNA7","mRNA8"]
#r.simulate(0,300,100,mRNA)
#r.plot(figsize=[7,7],xlim=(0,300),linewidth=2) 

#csv_filename="../Biotapestry/8gene_test.csv"
#csv_newfile="../Biotapestry/8gene_test2.csv"

#for i in range(1,9):
#    for j in range(1,9):
#        add=[(i,4,-1),(j,4,-1)]
#        print(i,j)
#        r.resetToOrigin()
#        add_biotapestry(add,csv_filename,csv_newfile)
#        ant_str = convert_biotapestry_to_antimony("../Biotapestry/8gene_test2.csv", 8, [8.4111E-05, 5.919E-01, 1.88693E-01, 4.035919E-01,
#                                                                                        5.38276E-02, 5.92, 2.9811E-02])
#        r=te.loada(ant_str)
#        r.simulate(0,300,100,mRNA)
#        r.plot(figsize=[7,7],xlim=(0,300),linewidth=2,linestyle='-') 

