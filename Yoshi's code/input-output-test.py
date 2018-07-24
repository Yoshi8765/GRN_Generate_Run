# -*- coding: utf-8 -*-
"""
Created on Thu Jul 05 12:28:35 2018
Made in Python 2.7

@author: Yoshi
"""

#todo:
#change to correct CRN scheme
#randomize occurance of same motif based on literature
#user-input
#steady-state check

import os
import glob
import tellurium as te
import numpy as np
import antimony
import string
import re

np.set_printoptions(linewidth=160)
#%% Combining multiple motif Antimony models together using submodules

### Importing in motifs from a motifs folder in working directory.
print('Looking for motifs folder in '+ str(os.getcwd()) + " ...\n")
motifList = glob.glob(".\\motifs\\test\\*.txt")

numReads = len(motifList)
print('There are '+ str(numReads) + ' motifs.\n')

reads = {}
for i,ascii in enumerate(string.ascii_uppercase[:numReads],0):
    reads[ascii] = te.readFromFile(motifList[i])

### Initialize
models = list([reads.values()[0]])
species = []
currStr = reads.values()[0]
splitStr = re.split('(model)( [\*]?)(.*?\(*\))*\\n\\n',currStr)
modules = list([reads.keys()[0] + " : " + splitStr[3] + "; \n"])
modelList = list([splitStr[3]])
inputs = []
outputs = []

### Construct the combined model Antimony Code
for i in np.arange(1,numReads):
    currStr = reads.values()[i]
    splitStr = re.split('(model)( [\*]?)(.*?\(*\))*\\n\\n',currStr)
    if splitStr[3] not in modelList:
        models.append(reads.values()[i])
    modelList.append(splitStr[3])
    #Note: the modules have to call the model name within the .txt file for each motif, not the name of the .txt file
    modules.append(reads.keys()[i] + " : " + splitStr[3] + "; \n") #i
    species.append("var species p_c"+str(i)+"; \n") #i-1
    
    inputs.append(reads.keys()[i-1] + ".p_input" + " is p_c" + str(i) + "; \n") #i-1
    outputs.append(reads.keys()[i] + ".p_output" + " is p_c" + str(i) + "; \n") #i-1

#combine lists to create a final 'combined' Antimony model 
combined = "".join("".join(models) + 'model combined()\n'+ "".join(modules) + "".join(species) + "".join(inputs) + "".join(outputs) + 'end')

#flatten the combined model by converting it to sbml and then converting back to Antimony
antimony.clearPreviousLoads()
antimony.loadAntimonyString(combined)
sbmlStr = antimony.getSBMLString('combined')
antimony.loadSBMLString(sbmlStr)
flatComb = antimony.getAntimonyString('combined')

r = te.loada(flatComb)             # If it fails here, the combined model is not formatted correctly
r.exportToAntimony('combined.txt') # Save to combined.txt in working directory

### Optional Outputs
#r.draw(layout='dot')


### What I want

#combined = add(all models I have) + '''
#model combined
#    #repeat for n-1 modules
#    #!!module names have to be identical to model names from imported models!!
#    A : model1();
#    B : model2();
#    #repeat for all models I have
#    #specify a global input node
#    var species p_c;
#    A.p_input is p_c;
#    B.p_output is p_c;
#    #repeat for n-1 modules
# end 
# '''

# Code to fix the same motif problem with libsbml (not preferred, as it renames submodel IDs)
#import libsbml
#    if splitStr[3] not in modules:
#        antimony.clearPreviousLoads()
#        antimony.loadAntimonyFile(motifList[i])
#        currSbmlstr = antimony.getSBMLString(antimony.getMainModuleName())
#        currSbmldoc = libsbml.readSBMLFromString(currSbmlstr)
#        model = currSbmldoc.getModel()
#        model.getId()
#        model.setId('Test')
#        model = libsbml.writeSBMLToString(currSbmldoc)
#        currStr = te.sbmlToAntimony(model)