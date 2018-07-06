# -*- coding: utf-8 -*-
"""
Created on Thu Jul 05 12:28:35 2018

@author: Yoshi
"""

import os
import tellurium as te
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(linewidth=160)

#%% Feed-Forward Loop Motif

print('\n')
print('*')*70
print('\n')

r = te.loada('''
//create a FFL: P1->P2->P3 ; P1->P4->P3
//assume all activation TF 
//assume time units are in minutes

model ffl()

#todo: needs to be more user-friendly to read (ex: add species declaration like below)
#//species declarations
#gene g1, g2 g3, g4;
#species p1, p2, p3, p4, m1, m2, m3, m4;

//reactions

r1: => m1 ; a_m1 * g1 + L1
r2: m1 => ; d_m1*m1
r3: => p1 ; a_p1*m1 
r4: p1 => ; d_p1*p1

r5: => m2 ; v1 * g2 + L2
r6: m2 => ; d_m2*m2
r7: => p2 ; a_p2*m2
r8: p2 => ; d_p2*p2

//r98: => m3 ; v2*g3
//r99: => m3 ; v3*g3

r9:  => m4 ; v1 * g4 + L3
r10: m4 => ; d_m4*m4
r11: => p4 ; a_p4*m4
r12: p4 => ; d_p4*p4

r13: => m3 ; v2 * g3 + L4
r14: m3 => ; d_m3*m3
r15: => p3 ; a_p3*m3
r16: p3 => ; d_p3*p3

//reaction rates
v1 := (p1^h1)/(K1^h1 + p1^h1)
v2 := Vm2 * (p2^h2)/(K2^h2 + p2^h2) * (Vm3 * (p4^h3)/(K3^h3 + p4^h3))

//species and parameters
g1 = 1 ; g2 = 1  ; g3 = 1  ; g4 = 1  ; //DNA count
m1 = 0 ; m2 = 0  ; m3 = 0  ; m4 = 0  ; //mRNA count
p1 = 0 ; p2 = 0  ; p3 = 0  ; p4 = 0  ; //protein count

//a_m, d_m from bionumbers DNA binding model

L1 = .1     ; L2 = .1     ; L3 = .1     ; L4 = .1    ; //leak rate
a_m1 = 20   ; a_m2 = 20   ; a_m3 = 20   ; a_m4 = 20  ; //transc rate
a_p1 = 1    ; a_p2 = 1    ; a_p3 = 1    ; a_p4 = 1   ; //transl rate
d_m1 = .6   ; d_m2 = .6   ; d_m3 = .6   ; d_m4 = .6  ; //mRNA deg rate
d_p1 = 0.09 ; d_p2 = 0.09 ; d_p3 = 0.09 ; d_p4 = 0.09; //prot deg rate

Vm1 = 20 ; Vm2 = 20 ; Vm3 = 20 ; //max transc rate
h1 = 2   ; h2 = 2   ; h3 = 2   ; //hill coeff
K1 = .2  ; K2 = .2  ; K3 = .2  ; //dissociation constant
''')

r.reset()
res = r.simulate(0,50,1000)
r.plot()