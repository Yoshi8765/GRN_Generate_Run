# -*- coding: utf-8 -*-
"""
Created on Thu Jul 05 12:28:35 2018
Made in Python 2.7

@author: Yoshi
"""

#TODO:randomize occurance of same motif based on literature
#---------fan-out
#---------config file
# add autoregulation and long-range connections, and cascades
#optional insertion of an oscillator or bistable switch
# - oscillator (have both types!)
# -- use ML objective functions to optimize parameters until you get an oscillator
# -- check that it is sustained using eigenvalues
# -- veronica is doing this

import os
import glob
import tellurium as te
import numpy as np
import antimony
import string
import re
import importlib
import itertools
import collections
import inspect
import time
#import traceback
import sys

#### SUPRESS TRACEBACK ####
sys.tracebacklimit = 0

np.set_printoptions(linewidth=160)

#%% Functions used

def motifsFromSource():
    print('Looking for \'motifs\' folder in '+ str(os.getcwd()) + " ...\n")
    motifList = glob.glob('./motifs/*.txt')
    if motifList == []:
        print('\nCannot find scripts in motifs folder! Defaulting to second option...\n')
        time.sleep(2)
        motifsFromScratch()
    return(motifList)

def motifsFromScratch():
    motifListDir = glob.glob("fxn_motif*.py") #finds all py files that start with fxn_motif
    if motifListDir == []:
        raise ImportError('Cannot find motif creation scripts in current working directory!')    
    print('===Motif List===')
    for i in range(len(motifListDir)):
        print(str(i)+'.'),
        print(motifListDir[i]+'\n'),
    print('================')
    motifInc(0,motifListDir)
    return()

def motifInc(counter,motifListDir):
    if counter == 3:
        raise Exception('Too many invalid inputs. Exiting...')
    counter += 1
    for retry in range(3):
        try:
            chosen = input('Which motifs do you want to include? [Comma separated list]\n')
            if isinstance(chosen, (int, long)): #turn an input with only one motif into a digit list
                chosen = map(int, str(chosen))
            if all(elem in range(len(motifListDir)) for elem in chosen):
                for i in chosen:
                    print('Generating motif from:' + str(motifListDir[i])+'\n')
                    currmotifDir = str(motifListDir[i])    
                    fxn = importlib.import_module(currmotifDir[:-3])
                    functions = inspect.getmembers(fxn,inspect.isfunction)
                    try:
                        test = functions[2] #EAFP test to see if the tuple `function` has more than one element
                        raise ValueError('The script ' + str(motifListDir[i]) + ' has more than one function (Contains: '+str(test)+'). Recheck script format.')
                    except:
                        method = getattr(fxn,functions[0][0])
                        try:
                            model = method() 
                            motifList.append(model)
                        except RuntimeError as error: 
                            raise
                        except:
                            print('Something went wrong!')
                            raise
                return()
            else:
                print('Entry includes indexes that don\'t exist. Try again.')
        except AttributeError: #Error caught from fxn.*
            raise AttributeError('model function doesn\'t exist. Please check model ' + str(motifListDir[i])+'.')
        except RuntimeError: #Error caught from fxn.*
            print('ERROR: Please check fxn ' + str(motifListDir[i]) + '.\n')
            raise RuntimeError('From fxn: ' + str(error))
        except:
            print('Invalid entries in input. Please enter index number of motifs.')
    else:
        raise ValueError('Too many invalid Entries.')
    return()


def chooseDupMotif(motifListDir):
    #Print motifs in use.
    print('===Included Motifs===')
    for i in range(len(motifList)):
        print(str(i)+'.'),
        print(motifList[i]+'\n'),
    print('=====================')
    #Takes care of duplicating motifs
    for retry in range(3): # This for loop allows for retries in case of invalid user entries.
        try:
            selMotifs = input('Which motifs do you want to use multiple times? (Comma-separated entries, or \'none\')\n')
            if selMotifs == 'none':
                break
            if isinstance(selMotifs, (int, long)):
                selMotifs = map(int, str(selMotifs))
            if all(elem in range(len(motifList))  for elem in selMotifs):
                dupMotif(selMotifs,motifList)
                return()
            else:
                print('Entry includes indexes that don\'t exist. Try again.')
        except:
            print('Invalid entries in input. Please enter index number of motifs.')    
    else:
        raise ValueError('Too many invalid Entries.')
    return()

def dupMotif(selMotifs,motifList):
    if set(selMotifs).issubset(range(len(motifList))):
        for j in range(len(selMotifs)):
            for retry1 in range(3):
                print('Current Motif:' + str(motifList[selMotifs[j]]))
                for retry2 in range(3):
                    try:
                        repNum = input('How many times do you want to use this motif?\n')
                        if repNum == int(repNum) and repNum>0:
                            for k in range(repNum-1):
                                motifList.append(motifList[selMotifs[j]])
                            print('Motifs copied.\n')
                        else:
                            print('Please enter a positive integer.')
                        break
                    except:
                        print('Invalid entries in input.')
                        repNum = []
                    break
                else:
                    raise ValueError('Too many invalid Entries.')
                    return()
                break
            else:
                raise ValueError('Too many invalid Entries.')
                return()
    else:
        print('Invalid entries in input. Please enter index of motifs.')
    return() #TODO: return the number of times each motif was copied.

def combineMotifs():
    #Summary of motifs included
    numReads = len(motifList)
    print('-------------------------------')
    print('There are '+ str(numReads) + ' motifs.')
    #print(motifCopiedNum) ###See above TODO
    print('-------------------------------')
    
    #Label all models with an Alphabetically ordered key
    reads = collections.OrderedDict({})
    for i,ascii in zip(list(range(numReads)), product_gen(string.ascii_uppercase)):
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
        species.append("var species p_c"+str(i)+"; \n") #i
        
        inputs.append(reads.keys()[i-1] + ".p_i" + " is p_c" + str(i) + "; \n") #i-1
        outputs.append(reads.keys()[i] + ".p_o" + " is p_c" + str(i) + "; \n") #i
    
    #combine lists to create a final 'combined' Antimony model 
    combined = "".join("".join(models) + 'model combined()\n'+ "".join(modules) + "".join(species) + "".join(inputs) + "".join(outputs) + 'end')
    return(combined)

def flattenMotif(combined):
    #flatten the combined model by converting it to sbml and then converting back to Antimony
    antimony.clearPreviousLoads()
    code = antimony.loadAntimonyString(combined)
    if code <= 0:
        textfile = open('combined.txt', 'w')
        textfile.write(combined)
        textfile.close()
        raise AssertionError('combined model is not flattenable. Are you using the right species names? ' + 'Combined model saved as: ' + str(os.getcwd()) + '\\combined.txt')
    sbmlStr = antimony.getSBMLString('combined')
    antimony.loadSBMLString(sbmlStr)
    flatComb = antimony.getAntimonyString('combined')
    
####TODO
#Delete extraneous mRNA -> Protein reactions at P_c connections
#look for p_c#, regex
#delete second equation with occurance of p_c#curr, regex

#todo: remove extraneous species + parameters in the removed equation

####TODO
    return(flatComb)

def testCombMotif(flatComb,motifList):
    try:
        combMotif = te.loada(flatComb)    # If it fails here, the combined model is not formatted correctly
        SS = combMotif.steadyState()
    except:
        raise AssertionError('Failed to create a combined model! Check code.')
    try:
        SS = combMotif.steadyState()
        print('Model reaches Steady-State!')
    except:
        print('Model Fails to reach Steady-state. Retrying with new parameter values.')
        global SSFAILED
        SSFAILED = True
        motifList=[]
        motifListDir = motifsFromScratch()
        motifList = chooseDupMotif(motifListDir)
        combined = combineMotifs()
        flatComb = flattenMotif(combined)
        testCombMotif(flatComb,motifList)
    if SS == 0:
        print('Steady-state = 0. Retrying with new parameter values.')
        global SSFAILED
        SSFAILED = True
        motifList=[]
        motifListDir = motifsFromScratch()
        chooseDupMotif(motifListDir)
        combined = combineMotifs()
        flatComb = flattenMotif(combined)
        testCombMotif(flatComb,motifList)
    return(combMotif,SSFAILED)

def product_gen(n):
    for r in itertools.count(1):
        for i in itertools.product(n, repeat=r):
            yield "".join(i)

def drawOutput(combMotif):
    ans = raw_input('Do you want to draw the model (Requires working PyGraphViz)? y/n\n')
    if ans == 'y':
        combMotif.draw(layout='dot')
        return()
    elif ans == 'n':
        return()
    else:
        print('Invalid entry in input.')
        drawOutput(combMotif)
        return()
    return()

####################################################################################
#%% Main Method: Combining multiple motif Antimony models together using submodules#
####################################################################################

global SSFAILED
SSFAILED = False   #Tag for when the model fails to run.

for retry in range(3): # This for loop allows for retries in case of invalid user entries.
    motifSource = raw_input('Do you want to import motifs from a folder?: y/n \n')
    motifList = []

# Option 1:
### Importing in motifs from a motifs folder in working directory.
    if motifSource == 'y':
        motifList = motifsFromSource()
        break

# Option 2:
### Generate motifs on the spot
    if motifSource == 'n':
        motifListDir = motifsFromScratch()
        break
    else:
        print('Invalid Entry.')
        
else:
    raise ValueError('Too many invalid Entries.')

# Create duplicates of motifs if necessary
chooseDupMotif(motifListDir)

# Create Combined Motif
combined = combineMotifs()

#Flatten the combined motif
flatComb = flattenMotif(combined)

#Test if Combined Motif Works
combMotif,SSFAILED = testCombMotif(flatComb,motifList)

### Optional Outputs
drawOutput(combMotif)

#Export Model
print('Combined model created.')
combMotif.exportToAntimony('combined.txt') # Save to combined.txt in working directory
print('Combined model saved as: ' + str(os.getcwd()) + '\\combined.txt')

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