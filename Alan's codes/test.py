# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 17:07:21 2018

@author: Yoshi
"""
import tellurium as te

r= te.loada('''
         model Random_GRN()
Rm1:=> M1;L1 + TRa1*((1+K2*P2 + K3*P3)/(1 +K2*P2 + K3*P3 + Kd3*P2*P3)) - M1*Dm1;
#Rm2:=> M2;L2 + TRa2*((K2*K3*P2*P3)/(1 +K2*P2 + K3*P3 + K2*K3*P2*P3)) - M2*Dm2;
Rm3:=> M3;L3 + TRa3*(K3*P3/(1 +K3*P3 + K2*P2 + Kd2*P3*P2)) - M3*Dm3;
Rp1:=> P1;TRb1*M1 - P1*Dp1;
Rp2:=> P2;TRb2*M2 - P2*Dp2;
Rp3:=> P3;TRb3*M3 - P3*Dp3;
P1 = 0.731;
P2 = 0.355;
P3 = 0.196;
Dm1 = 0.384;
Dm2 = 0.629;
Dm3 = 1.277;
Kd1 = 0.275;
Kd2 = 0.575;
Kd3 = 0.284;
TRa1 = 29.797;
TRa2 = 30.13;
TRa3 = 30.33;
H1 = 0.766;
H2 = 1.12;
H3 = 0.972;
K1 = 0.506;
K2 = 0.988;
K3 = 0.458;
TRb1 = 30.147;
TRb2 = 30.165;
TRb3 = 30.296;
M1 = 0;
M1 = 0;
M2 = 0;
M2 = 0;
M3 = 0;
M3 = 0;
L1 = 1.11;
L2 = 0.848;
L3 = 0.956;
Dp1 = 0.846;
Dp2 = 0.459;
Dp3 = 0.259;
end
         '''
         )

r.simulate(0,100,200)
r.plot()
r.steadyState()