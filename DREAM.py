# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 15:00:49 2018

@author: Yoshi
"""

import os
import tellurium as te
import matplotlib as plt
import structural as st
import numpy as np

np.set_printoptions(linewidth=160)

#%% Model 1
r.reset()

r=te.loadSBMLModel('model1_0\model_without_parameters\model1.sbml')
ant = te.sbmlToAntimony('model1_0\model_without_parameters\model1.sbml')
te.saveToFile('model1_0\model_without_parameters\model1_antimony.txt',ant)
a = r.simulate(0,10,100)

time = a[:,0]
output = a[:,1::]
print 'Boundary Species:'
print r.model.getBoundarySpeciesIds()
print 'Floating Species:'
print r.model.getFloatingSpeciesIds()
print 'Global Parameters:'
print r.model.getGlobalParameterIds()


ls = st.LibStructural()
ls.loadSBMLFromString(r.getSBML())
print(ls.getSummary())
print(ls.getTestDetails())
print ls.getStoichiometryMatrix()
print ls.getElementaryModes()

#%% Trying to make a model generator

print('\n')
print('*')*70
print('\n')

r.reset()

r = te.loada('''
//create a FFL: A->B->C ; A->D->C
//assume all activation TF 
$g1 -> m1 ; a_m1*g1 - d_m1
m1 -> p1 ; a_p1*m1 - d_p1*p1

$g2 -> m2 ; a_m2*g2 - d_m2
m2 -> p2 ; a_p2*m2 - d_p2*p2

$g3 -> m3 ; a_m3*g3 - d_m3
m3 -> p3 ; a_p3*m3 - d_p3*p3

$g4 -> m4 ; a_m4*g4 - d_m4
m4 -> p4 ; a_p4*m4 - d_p4*p4

a_m1=1 ; g1=1
a_p1=1 ; m1=1
d_m1=.1 ; p1=1
d_p1=.1 ; g2=1
a_m2=1 ; m2=1
a_p2=1 ; p2=1
d_m2=.1 ; g3=1
d_p2=.1 ; m3=1
a_m3=1 ; p3=1
a_p3=1 ; g4=1
d_m3=.1 ; m4=1
d_p3=.1 ; p4=1
a_m4=1 ;
a_p4=1 ;
d_m4=.1 ;
d_p4=.1 ;
''')

res = r.simulate(0,100,1000)
r.plot()