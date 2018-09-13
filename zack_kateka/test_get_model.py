# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 11:42:03 2018

@author: Kateka
"""

import tellurium as te

from GetModel import get_model

#antimony_str = get_model(5)
#r = te.loada(antimony_str)

#s2 = get_model(-1)

#s3 = get_model(5, reg_probs=["hi"])

s4 = get_model(20, reg_probs=[0.2, 0.2, 0.2, 0.2, 0.2], reachability= 0.9)

r = te.loada('''
model *pathway()

	// Compartments and Species:
	species INPUT, P1, mRNA1, P2, mRNA2, P3, mRNA3, P4, mRNA4, P5, mRNA5;

	// Assignment Rules (production rates used in reactions):
	// transcription1(SA : in connections = [P2]) uses production rate := Vm1*((K1_1*P2^H1)/(1 +K1_1*P2^H1));
	// transcription2(DA : in connections = [P4, P1]) uses production rate := Vm2*((K1_2*P4^H2 + K2_2*P1^H2 + K1_2*K3_2*P4^H2*P1^H2)/(1 + K1_2*P4^H2 + K2_2*P1^H2 + K1_2*K3_2*P4^H2*P1^H2));
	// transcription3(SA : in connections = [P2]) uses production rate := Vm3*((K1_3*P2^H3)/(1 +K1_3*P2^H3));
	// transcription4(DA : in connections = [P2, P5]) uses production rate := Vm4*((K1_4*P2^H4 + K2_4*P5^H4 + K1_4*K3_4*P2^H4*P5^H4)/(1 + K1_4*P2^H4 + K2_4*P5^H4 + K1_4*K3_4*P2^H4*P5^H4));
	// transcription5(DR : in connections = [INPUT, P2]) uses production rate := Vm5*( 1 /(1 + K1_5*INPUT^H5 + K2_5*P2^H5 + K1_5*K3_5*INPUT^H5*P2^H5));

	// Reactions:
	//transcription1
	J1: => mRNA1 ; L1 + Vm1*((K1_1*P2^H1)/(1 +K1_1*P2^H1)) - d_mRNA1 * mRNA1;
	//translation1
	F1: => P1 ; a_protein1 * mRNA1 - d_protein1 * P1;
	//transcription2
	J2: => mRNA2 ; L2 + Vm2*((K1_2*P4^H2 + K2_2*P1^H2 + K1_2*K3_2*P4^H2*P1^H2)/(1 + K1_2*P4^H2 + K2_2*P1^H2 + K1_2*K3_2*P4^H2*P1^H2)) - d_mRNA2 * mRNA2;
	//translation2
	F2: => P2 ; a_protein2 * mRNA2 - d_protein2 * P2;
	//transcription3
	J3: => mRNA3 ; L3 + Vm3*((K1_3*P2^H3)/(1 +K1_3*P2^H3)) - d_mRNA3 * mRNA3;
	//translation3
	F3: => P3 ; a_protein3 * mRNA3 - d_protein3 * P3;
	//transcription4
	J4: => mRNA4 ; L4 + Vm4*((K1_4*P2^H4 + K2_4*P5^H4 + K1_4*K3_4*P2^H4*P5^H4)/(1 + K1_4*P2^H4 + K2_4*P5^H4 + K1_4*K3_4*P2^H4*P5^H4)) - d_mRNA4 * mRNA4;
	//translation4
	F4: => P4 ; a_protein4 * mRNA4 - d_protein4 * P4;
	//transcription5
	J5: => mRNA5 ; L5 + Vm5*( 1 /(1 + K1_5*INPUT^H5 + K2_5*P2^H5 + K1_5*K3_5*INPUT^H5*P2^H5)) - d_mRNA5 * mRNA5;
	//translation5
	F5: => P5 ; a_protein5 * mRNA5 - d_protein5 * P5;

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

	// Variable initializations:
	H1 = 0.632931770488;
	Vm1 = 1.04197851992;
	L1 = 0.97064617655;
	d_mRNA1 = 35.5346820902;
	a_protein1 = 43.1925549623;
	d_protein1 = 0.461411387941;
	K1_1 = 0.933461981124;
	K2_1 = 0.933070835883;
	K3_1 = 1.27462918078;
	H2 = 0.58805147899;
	Vm2 = 1.03697354515;
	L2 = 0.626151010343;
	d_mRNA2 = 26.688763913;
	a_protein2 = 38.9758843306;
	d_protein2 = 0.404148297226;
	K1_2 = 1.02371192618;
	K2_2 = 0.857748488476;
	K3_2 = 1.28864115241;
	H3 = 0.516935894237;
	Vm3 = 1.07436893564;
	L3 = 0.940519097317;
	d_mRNA3 = 32.1991680568;
	a_protein3 = 15.0019925094;
	d_protein3 = 0.534198845158;
	K1_3 = 1.2699519038;
	K2_3 = 0.563280067641;
	K3_3 = 0.9685320833;
	H4 = 0.478703739015;
	Vm4 = 1.14803739216;
	L4 = 0.465281311194;
	d_mRNA4 = 32.5851868639;
	a_protein4 = 21.2206544058;
	d_protein4 = 0.429940380266;
	K1_4 = 0.960724445272;
	K2_4 = 1.11865834228;
	K3_4 = 0.981420910407;
	H5 = 0.559780284189;
	Vm5 = 1.24295704444;
	L5 = 0.876857235995;
	d_mRNA5 = 24.8165832972;
	a_protein5 = 29.2478262299;
	d_protein5 = 0.426231248139;
	K1_5 = 0.966073415871;
	K2_5 = 0.844702779634;
	K3_5 = 0.763326923443;

	// Other declarations:
	const H1, Vm1, L1, d_mRNA1, a_protein1, d_protein1, K1_1, K2_1, K3_1, H2, Vm2, L2, d_mRNA2, a_protein2, d_protein2, K1_2, K2_2, K3_2, H3, Vm3, L3, d_mRNA3, a_protein3, d_protein3, K1_3, K2_3, K3_3, H4, Vm4, L4, d_mRNA4, a_protein4, d_protein4, K1_4, K2_4, K3_4, H5, Vm5, L5, d_mRNA5, a_protein5, d_protein5, K1_5, K2_5, K3_5;


end



''')
r.simulate(0,50,1000)
r.plot()

#s5 = get_model(5)


#%% Running the model to see if it works
#r=te.loada(s5)
#r.reset()
#plt.close("all")
#res = r.simulate(0,50,1000)
#r.plot()
#r.draw(layout='fdp')

#print ("\n\n\ndone!")



