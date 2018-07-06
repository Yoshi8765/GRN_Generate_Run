# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 12:16:49 2018

@author: Alan
"""

import tellurium as te
import roadrunner
import matplotlib.pyplot as plt
import numpy as np
#   Activation    a3 + (P2/k3)/(1+P2/k3)^H3
#   Repression    a4 + 1/(1+ (P1/k4)^H4)
r=te.loada("""
    R0:      => M1   ; L1 + TM1 - dm1*M1
    R1:      => M2   ; L2 + TM2 * 1/(1+(P1/k2)^H2) - dm2*M2
    R2:      => M3   ; L3 + TM3 * (P2/k3)^H3/(1+P2/k3)^H3 - dm3*M3
    R3:      => M4   ; L4 + TM4 * 1/(1+ (P1/k4)^H4)* (P3/k4)^H4/(1+(P3/k4)^H4) - dm4*M4
    R8:      => P1   ; Tr1*M1 - dp1*P1
    R9:      => P2   ; Tr2*M2 - dp2*P2
    R10:     => P3   ; Tr3*M3 - dp3*P3
    R11:     => P4   ; Tr4*M4 - dp4*P4

     M1 = 0;      M2 = 0;      M3 = 0;      M4 = 0;
     P1 = 0;      P2 = 0;      P3 = 0;      P4 = 0;
     L1 = .01;    L2 = .01;    L3 = .01;    L4 = .01;
     k1 = .65;    k2 = .65;    k3 = .65;    k4 = .65;
    dm1 = .5;    dm2 = .5;    dm3 = .5;    dm4 = .5;
    dp1 = .5;    dp2 = .5;    dp3 = .5;    dp4 = .5;
    TM1 = 10;     TM2 = 10;    TM3 = 10;    TM4 = 10;
    Tr1 = .5;    Tr2 = .5;    Tr3 = .5;    Tr4 = .5;
    H1  =  1;    H2  =  1;    H3  =  1;    H4  =  2;
    
""")
result = r.simulate(0, 50, 200,)
plt.figure()
plt.grid(color='k', linestyle='-', linewidth=1)
plt.ylim(0,np.max(result[:,2:8]))
plt.yticks(np.arange(0,np.max(result[:,2:7]),np.max(result[:,2:8])/20))
#M1 , = plt.plot (result[:,0],result[:,1], label = 'M1')
#M2 , = plt.plot (result[:,0],result[:,3], label = 'M2')
#M3 , = plt.plot (result[:,0],result[:,5], label = 'M3')
#M4 , = plt.plot (result[:,0],result[:,6], label = 'M4')
P1 , = plt.plot (result[:,0],result[:,4], label = 'P1')
P2 , = plt.plot (result[:,0],result[:,2], label = 'P2')
P3 , = plt.plot (result[:,0],result[:,7], label = 'P3')
P4 , = plt.plot (result[:,0],result[:,8], label = 'P4')
plt.legend([P1, P2, P3, P4], ['P1', 'P2', 'P3', 'P4'])

#plt.legend([M1, M2, M3, M4, P1, P2, P3, P4], ['M1', 'M2', 'M3', 'M4', 'P1', 'P2', 'P3', 'P4'])
