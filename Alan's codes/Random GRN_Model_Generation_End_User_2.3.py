# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 14:52:12 2018

@author: Alan
"""
"""
2.3 Update:
    Cleaned up functions.
    Cleaned up verbose output.
    Renamed variables to match documentation and across function calls
    Removed extraneous code (making unnecessary variables)
    Added more error handling
    Made error outputs more visible
    Change user input to be forgiving to all entries 3 times
    Add verbose error mode
"""

# TODO:
#Figure out what to do about Perturbation method
#Figure out what to do about GetEvents method
#implement hill coefficient -> HILL COEFF NOT WORKING, BREAKS MODEL GENERATION
# TODO: There are still orphans being made.

import tellurium as te
import numpy as np
import math
import matplotlib.pyplot as plt
import os
import sys
import traceback
from collections import Counter
import logging
from time import sleep

######VERBOSE ERRORS SWITCH######
VERBOSE_ERRORS = True

#####Some useful variables to change when needed#####
GRAPH_LABEL_FONTSIZE = 5
GRAPH_TITLE_FONTSIZE = 8

NUCLEOTIDE_SOURCE = '' #'$N'
AMINO_ACID_SOURCE = '' #'$A'

    #Functions

    #Handles Errors
def ErrorPrinting(Error,customHeader='Error: '):
    if VERBOSE_ERRORS == True:
        errorStr = str(Error)
        print('!'*len(errorStr))
        print(customHeader + errorStr)
        print('!'*len(errorStr))
        sleep(.5)
        return()
    else:
        sleep(.5)
        return()

    #Creating a model that statically makes sense.
    #This whole model (and its submodels) are used to generate an Antimony model composed of various logic gates
def GetModel(tries,numGenes,regProb, InitParams, Name):
    print 'Model simulation tries = ' + str(tries)
    #Set Lists and Empty Dictionaries
    Interactions = ["SingleAct", "SingleRep", "and", "or", "Nand", "Nor", "Counter"]
    GRN = {'Rate': [], 'mRNA':[], 'Prot': [],'PDeg':[],'mDeg': [],'Leak': [],'Vm':[], 'a_p':[], 'DualDissoc':[], 'Dissoc':[], 'HillCoeff': []}
    GRNInt = {'TF':[],'TFs':[]}
    #Perturbation = {'TRa1':[10]}
    StringList=[]
    modelCreationTries = 0
    while modelCreationTries < 3:
        try:
            ModelInitNames(GRN)
            print 'ModelInitNames done.'
            AssignRegulation(numGenes,GRN,GRNInt,Interactions)
            print 'AssignRegulation done.'
            DetermineInteractors(GRN,GRNInt,regProb)
            print 'DetermineInteractors done.'
            ChooseProtein(GRNInt)
            print 'ChooseProtein done.'
            Rates = GetReactionRates(GRN,GRNInt)
            print 'GetReactionRates done.'
            GetReactions(GRN,StringList,Rates)
            print 'GetReactions done.'
            ValueGeneration(InitParams,StringList)
            #GetEvents(Perturbation,StringList)
            print 'ValueGeneration done.'
            antStr = GetModelString(StringList)
            print('Loading model into Antimony...\n')
            #print StringList
        except Exception as error:
            print "Something went wrong when constructing Model (GetModel method)!"
            ErrorPrinting(error)
            return(0,0,0)
        try:
            model = te.loadAntimonyModel(antStr)
            print('Success!')
            break
        except:
            modelCreationTries = modelCreationTries + 1
            print "Failed to load model code into Antimony."
            if antStr == '':
                return(0,StringList.join('This is a pre-loaded string into Antimony that failed.\n'),GRN)
            return(0,antStr,GRNInt)
    return (model, antStr,GRNInt)

#%%
    #Create species names in model
def ModelInitNames(GRN):
    for k in range(numGenes):
        num = str(k+1)
        GRN['Leak'].append("L" + num)
        GRN['Rate'].append('R' + num +' := ')
        GRN['mDeg'].append('d_m' + num)
        GRN['PDeg'].append('d_p' + num)
        GRN['mRNA'].append('M' + num)
        GRN['Prot'].append('P' + num)
        GRN['Vm'].append('Vm' + num)
        GRN['a_p'].append('a_p' + num)
        GRN['HillCoeff'].append('H' + num)
        GRN['Dissoc'].append('K' + num)
        GRN['DualDissoc'].append('Kd' + num)
    return

    #identifies activators among N genes in model. where N = # of numGenes
    #The remainder are set as repressors.
    # ISSUE: Currently any gene is set to either being an activator or a repressor, for any regulation.
def AssignRegulation(numGenes,GRN,GRNInt,Interactions):
    # Creating a list of activator proteins, numbering about half the number of genes.
    nums =[int(math.ceil(float(numGenes)/2)), int(math.floor(float(numGenes)/2))]
    RanAct = np.random.choice(GRN['Prot'], nums[np.random.randint(2)])
    #print "Randomly generated Activators: " + str(RanAct)
    GRNInt.update({'activators': RanAct})
    scannedActivators = []
    for i in (range(len(GRNInt['activators']))):
        ActivatorCur =GRNInt['activators'][i]
        if ActivatorCur in scannedActivators:
            i -= 1
            continue
        scannedActivators.append(ActivatorCur)
        i += 1
    GRNInt.update({'activators': scannedActivators})
    #print "Activators: " + str(GRNInt['activators'])

    #All proteins not turned into activators will be set as a repressor
    repressors = list(GRN['Prot'])
    for s in GRNInt['activators']:
        if s in repressors:
            repressors.remove(s)
        GRNInt.update({'repressors': repressors})
    #print "Repressors: " + str(GRNInt['repressors'])

    #Creating a list of regulation types, not assigned to proteins yet.
    ProbabilityMatrix = []
    ProbSing = regProb[0]/2
    ProbDub = regProb[1]/4
    ProbCounter = regProb[2]
    ProbabilityMatrix = [ProbSing, ProbSing, ProbDub, ProbDub, ProbDub, ProbDub, ProbCounter]
    pickedInteractions = np.random.choice(Interactions,numGenes,True,ProbabilityMatrix)
    interactionCounts = Counter(pickedInteractions)

    if interactionCounts['SingleRep'] > len(GRNInt['repressors']):
        numChangeInteractions = interactionCounts['SingleRep'] - len(GRNInt['repressors'])
        #print 'number of SingleRep behavior we have above capacity.: ' + str(numChangeInteractions)
        counter = 0
        singleRepCounter = 0
        for x in pickedInteractions: #when pickedInteraction == singlerep
            if numChangeInteractions == singleRepCounter:
                break
            if x == 'SingleRep':
                singleRepCounter += 1
            while x == 'SingleRep':
                x = np.random.choice(Interactions,1,True,ProbabilityMatrix)
                pickedInteractions.put(counter,x)
            counter += 1

    if interactionCounts['SingleAct'] > len(GRNInt['activators']):
        numChangeInteractions = interactionCounts['SingleAct'] - len(GRNInt['activators'])
        #print 'number of SingleAct behavior we have above capacity.: ' + str(numChangeInteractions)
        counter = 0
        singleActCounter = 0
        for x in pickedInteractions: #when pickedInteraction == singlerep
            if numChangeInteractions == singleActCounter:
                break
            if x == 'SingleAct':
                singleActCounter += 1
            while x == 'SingleAct':
                x = np.random.choice(Interactions,1,True,ProbabilityMatrix)
                pickedInteractions.put(counter,x)
            counter += 1

    GRNInt.update({'interactions': pickedInteractions})
    #print "Regulation types: " + str(GRNInt['interactions'])
    return


    #Generates a N set of interactions
    #Correct number and type are randomly assigned.
    #Assignments prevent dual regulation.
def DetermineInteractors(GRN,GRNInt,regProb):
    for Inter in range(len(GRNInt['interactions'])):
        reg = GRNInt['interactions'][Inter]
        RegType = 'repressors'
        if 'Single' in reg :
            if 'Act' in reg:
                RegType = 'activators'
            pickedProt = np.random.randint(len(GRNInt[RegType]))
            currTF = [reg, [GRNInt[RegType][pickedProt]] ]
            testcounter = 0
            while currTF in GRNInt['TF']:
                testcounter += 1
                pickedProt = np.random.randint(len(GRNInt[RegType]))
                currTF = [reg, [GRNInt[RegType][pickedProt]] ]
            GRNInt['TF'].append(currTF)
        else:
            Tps={'activators':[], 'repressors':[]}
            for m in Tps.keys():
                n=1
                Tps[m].append(np.random.choice(GRNInt['activators'],2))
                while Tps[m][0][n-1] == Tps[m][0][n]:
                    Tps[m][0][n] =np.random.choice(GRNInt['activators'])
            TypeA = 'activators'
            TypeR = 'repressors'
            if reg == "and" or reg == "or":
                GRNInt['TF'].append([reg, [Tps[TypeA][0][0],Tps[TypeA][0][1]]])
            elif reg == "Nand" or reg == "Nor":
                GRNInt['TF'].append([reg, [Tps[TypeR][0][0],Tps[TypeR][0][1]]])
            else:
                GRNInt['TF'].append([reg, [Tps[TypeA][0][0],Tps[TypeR][0][1]]] )
    return

    # Assigns network regulation between genes
def ChooseProtein(GRNInt):
    for i in range(numGenes):
        TF=[]
        string = GRNInt['TF'][i][1]
        for n in string:
            TF.append(int(n.replace("P","")))
        GRNInt['TFs'].append(TF)
    return

    # Generate mRNA production rates according to regulatory interaction type
def GetReactionRates(GRN,GRNInt):
    ReactionRates=[]
    WorkingString=''
    for i in range(numGenes):
        Reg = GRNInt['TF'][i][0]
        WorkingString += 'Rm_' + str(i+1) + ' := '
        TF1 = GRNInt['TFs'][i][0] - 1
        P1 = GRN['Prot'][TF1]
        K1 =  GRN['Dissoc'][TF1]
#        HCoeff = GRN['HillCoeff'][i]
        num = ''
        denom = ''
        #frac = ''
        if 'Single' in Reg:
            #frac = '(' + P1 + '*' + K1 + ')'
#            if 'Act' in Reg:
#                hillEq = '(' + frac + '/(1 + ' + frac + ')**' + HCoeff + ')'
#            else:
#                hillEq = '(1 /(1 + ' + frac + ')**' + HCoeff + ')'
            WorkingString = GRN['Leak'][i] + ' + ' + GRN['Vm'][i] + '- ' + GRN['mRNA'][i] + '*' + GRN['mDeg'][i] + ';\n'
        elif len(GRNInt['TFs'][i]) >1 :
            TF2 = GRNInt['TFs'][i][1] - 1

            P2 = GRN['Prot'][TF2]
            K2 =  GRN['Dissoc'][TF2]
            K3 = GRN['DualDissoc'][TF2]
            if 'Counter' in Reg:
                num = K1 + '*' + P1
                denom = '(1 +' + K1 + '*' + P1 + ' + ' + K2 + '*' + P2 + ' + ' + K3 + '*' + P1 + '*' + P2 + ')'#K3 is a place holder.. need to generate dissoc associated with dual reg
                eq = '(' + num  + '/' + denom + ')'

            elif  not('N' in Reg) and("and" in Reg or 'or' in Reg):
                if 'and' in Reg:
                    num = '(' + K1 + '*' + K2 + '*' + P1 + '*' + P2  + ')'
                    denom = '(1 +' + K1 + '*' + P1 + ' + ' + K2 + '*' + P2 + ' + ' + K1 + '*' + K2 + '*' + P1 + '*' + P2 + ')'
                    eq = '(' + num  + '/' + denom + ')'
                else:
                    num = '(' + K1 + '*' + P1 + ' + ' + K2 + '*' + P2 + ')'
                    denom = '(1 +' + K1 + '*' + P1 + ' + ' + K2 + '*' + P2  + ')'
                    eq = '(' + num  + '/' + denom + ')'
            else:
                denom = '(1 +' + K1 + '*' + P1 + ' + ' + K2 + '*' + P2 + ' + ' + K3+ '*' + P1 + '*' + P2 + ')'
                if 'and' in Reg:
                    num = '(1+' + K1 + '*' + P1 + ' + ' + K2 + '*' + P2 + ')'
                    eq = '(' + num  + '/' + denom + ')'
                else:
                    eq =  '(1' + '/' +  denom + ')'
            WorkingString = GRN['Leak'][i] + ' + ' + GRN['Vm'][i] + '*' + eq + ' - ' + GRN['mRNA'][i] + '*' + GRN['mDeg'][i]+ ';\n'
        ReactionRates.append(WorkingString)
        #print 'ReactionRates ' + str(i) + ' :'+ str(ReactionRates)
        #print 'Reg ' + str(i) + ' :'+ str(Reg)
    print '\n'
    return ReactionRates

    # Generates chemical network reaction model for the model
def GetReactions(GRN,StringList,Rates):
    ReactionM =''
    ReactionP =''
    ReactionSumString= 'model Random_GRN()\n'
    for i in np.arange(1, numGenes+1):
        M = GRN['mRNA'][i-1]
        P = GRN['Prot'][i-1]
        ReactionM += 'Rm' + str(i) + ':' + NUCLEOTIDE_SOURCE + '=> ' + M + ';' + Rates[i-1]
        ReactionP += 'Rp' + str(i) + ':' + AMINO_ACID_SOURCE + '=> ' + P + ';' + GRN['a_p'][i-1] + '*' + M + ' - ' + P + '*' + GRN['PDeg'][i-1] + ';\n'
    ReactionSumString += ReactionM +ReactionP
    StringList.append(ReactionSumString)

    #Generates a random number close to a given value for parameters and initial conditions
def ValueGeneration(InitParams,StringList):
    Species=''
    Values = {'M':[InitParams[0]], 'P':[InitParams[1]],'d_p':[InitParams[2]],'d_m':[InitParams[3]],'L':[InitParams[4]],'Vm':[InitParams[5]], 'a_p':[InitParams[6]],'Kd':[InitParams[7]], 'K':[InitParams[8]], 'H': [InitParams[9]]}
    Keys = Values.keys()
    for i in Keys:
        for n in np.arange(1, numGenes+1):
            val = 0
            if i == 'M' or i=='4P':
                Values[i].append(val)
                Species += i + str(n) + ' = ' + str(val) + ';\n'
            else:
                while val <= 0:
                    val = round(np.random.normal(Values[i][0],.25),3)
            Values[i].append(val)
            Species += i + str(n) + ' = ' + str(val) + ';\n'
    StringList.append(Species)
    return

#%%Events-- attempt to perturb a user defined gene's transcription at the associated time
def GetEvents(Perturbation,StringList):
    EventString = ''
#        for m in range(numGenes):
#            EventString += 'Em' + str(m+1) + ': at(' + GRN['mRNA'][m] + '<0): ' + GRN['mRNA'][m] + '=0;\n'
#            EventString += 'Ep' + str(m+1) + ': at(' + GRN['Prot'][m] + '<0): ' + GRN['Prot'][m] + '=0;\n'    Old Catch. Previous issue.
    counter = 0
    for k in Perturbation.keys():
        counter +=1
        EventString += 'Ed' + str(counter) + ': at(time >'  + str(Perturbation[k][0]) +'): ' + k + '= ' + str(39) +';\n'
    StringList.append(EventString)
    return

def GetModelString(StringList):
    antStr = ''
    for i in StringList:
        antStr += i
    antStr += 'end'
    return antStr

#%%

    #This function tries to run the constructed model
def RunModel(runAttempts,tries,model,antStr,tmax,regProb,Percent,modelName,DataOut,seed, filePath):
    if tries <= runAttempts:
        try:
            plt.close('all')
            model.reset()
            tStep = int(math.ceil(tmax)*5)
            result = model.simulate(0,tmax, tStep)
            tries = tries + 1
        except:
             tries = tries + 1
             print "Failed to solve Model"
             model,antStr,GRNInt = GetModel(tries,numGenes,regProb, InitParams, modelName)
             RunModel(runAttempts,tries,model,antStr,tmax,regProb,NLevel,modelName,DataOut,seed, filePath)
        else:
            NoisyResult = np.zeros([len(result[:,0]), len(result[0,:])])
            for k in range(len(result[:,0])):
                for i in range(len(result[0])):
                    if i == 0:
                        NoisyResult[k,i] = result[k,i]
                    else:
                        CurrVal = -1
                        while CurrVal < 0:
                            CurrVal = result[k,i] + np.random.normal(0,Percent*result[k,i])
                        NoisyResult[k,i] = CurrVal
            return(DataOut, model, antStr,NoisyResult, result)
    else:
        raise ValueError('Could not create working model more than ' + str(runAttempts+1) + ' times!')
    return()

    #This function makes two time-course plots: One based strictly on the model (non-noisy) and one with noise artificially added in (noisy)
def makePlots(tmax,result,NoisyResult):
    global GRAPH_LABEL_FONTSIZE
    global GRAPH_TITLE_FONTSIZE
    jet = plt.get_cmap('jet')
    pickColorSpace = [int(np.floor(i)) for i in np.linspace(0,255,numGenes)]
    colors = jet(pickColorSpace)
    colors = jet(np.linspace(0,1,numGenes))
    vars = {'P':[],'M':[]}
    for e in vars.keys():
        plt.figure(1)
        for k in np.arange(1,numGenes+1):
            vars[e].append(e + str(k))
            tStep = int(math.ceil(tmax)/20)
            tRange = np.arange(0,(int(math.ceil(tmax)))+ tStep,tStep)
            if e=='P':
                plt.subplot(2,1,1)
                plt.title('Protein Count Vs. Time', fontsize=GRAPH_TITLE_FONTSIZE)
                plt.grid(color='k', linestyle='-', linewidth=.4)
                plt.ylabel('count',fontsize=GRAPH_LABEL_FONTSIZE)
                plt.yticks(fontsize=GRAPH_LABEL_FONTSIZE)
                plt.plot (result[:,0],result[:,(k-1)+1], label = vars[e][k-1],color=colors[k-1])
                plt.yscale('log')
                if tmax > 5000:
                    plt.xscale('log')
                    plt.xticks(fontsize=GRAPH_LABEL_FONTSIZE)
                else:
                    plt.xticks(tRange, fontsize=GRAPH_LABEL_FONTSIZE)
                plt.xlim(0,tmax)
                plt.legend(loc='upper right',ncol=3, bbox_to_anchor=(1.13, 1.035), fontsize=GRAPH_LABEL_FONTSIZE)
            if e=='M':
               plt.subplot(2,1,2)
               plt.title('mRNA Count Vs. Time', fontsize=GRAPH_TITLE_FONTSIZE)
               plt.grid(color='k', linestyle='-', linewidth=.4)
               plt.xlabel('time(s)',fontsize=GRAPH_LABEL_FONTSIZE)
               plt.ylabel('count',fontsize=GRAPH_LABEL_FONTSIZE)
               plt.yticks(fontsize=GRAPH_LABEL_FONTSIZE)
               plt.plot (result[:,0],result[:,(k-1)+numGenes+1], label = vars[e][k-1],color=colors[k-1])
               plt.yscale('log')
               if tmax > 5000:
                    plt.xscale('log')
                    plt.xticks(fontsize=GRAPH_LABEL_FONTSIZE)
               else:
                   plt.xticks(tRange, fontsize=GRAPH_LABEL_FONTSIZE)
               plt.xlim(0,tmax)
               plt.legend(vars[e],loc='upper right',ncol=3, bbox_to_anchor=(1.13, 1.035), fontsize=GRAPH_LABEL_FONTSIZE)
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    manager.window.showMinimized()
    plt.savefig(modelPath + 'Simulation_Plot.png', dpi=400)

    vars = {'P':[],'M':[]}
    for e in vars.keys():
        plt.figure(2)
        for k in np.arange(1,numGenes+1):
            vars[e].append(e + str(k))
            if e=='P':
                plt.subplot(2,1,1)
                plt.title('Protein Count Vs. Time (Noisy)', fontsize=GRAPH_TITLE_FONTSIZE)
                plt.grid(color='k', linestyle='-', linewidth=.4)
                plt.ylabel('count',fontsize=GRAPH_LABEL_FONTSIZE)
                plt.xlim(0,tmax)
                plt.yticks(fontsize=GRAPH_LABEL_FONTSIZE)
                plt.plot (NoisyResult[:,0],NoisyResult[:,(k-1)+1], label = vars[e][k-1],color=colors[k-1])
                plt.legend(vars[e])
                plt.yscale('log')
                if tmax > 5000:
                    plt.xscale('log')
                    plt.xticks(fontsize=6)
                else:
                   plt.xticks(tRange, fontsize=GRAPH_LABEL_FONTSIZE)
                plt.legend(loc='upper right',ncol=3, bbox_to_anchor=(1.13, 1.035), fontsize=GRAPH_LABEL_FONTSIZE)
            if e=='M':
                plt.subplot(2,1,2)
                plt.title('mRNA Count Vs. Time (Noisy)', fontsize=GRAPH_TITLE_FONTSIZE)
                plt.grid(color='k', linestyle='-', linewidth=.4)
                plt.xlim(0,tmax)
                plt.yticks(fontsize=GRAPH_LABEL_FONTSIZE)
                plt.xlabel('time(s)',fontsize=GRAPH_LABEL_FONTSIZE)
                plt.ylabel('count',fontsize=GRAPH_LABEL_FONTSIZE)
                plt.plot (NoisyResult[:,0],NoisyResult[:,(k-1)+numGenes+1], label = vars[e][k-1],color=colors[k-1])
                plt.legend(vars[e])
                plt.yscale('log')
                if tmax > 5000:
                    plt.xscale('log')
                    plt.xticks(fontsize=GRAPH_LABEL_FONTSIZE)
                else:
                    plt.xticks(tRange, fontsize=GRAPH_LABEL_FONTSIZE)
                plt.legend(loc='upper right',ncol=3, bbox_to_anchor=(1.13, 1.035), fontsize=GRAPH_LABEL_FONTSIZE)
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    manager.window.showMinimized()
    plt.savefig(modelPath + 'Noisy Simulation_Plot.png', dpi=400)
    print('\nPlots created!\n')
    return()

    #This function saves the data generated by the program into csv, txt, and xml files
def Output(data, Name,model,seed,result,NoisyResult, fileName, antStr):
    for i in range(len(data)):
        if i == 0 and data[i] == 'y':
            writecsvFile(fileName + Name + '_Results.csv',model,result)
        elif data[i] == 'y':
            writecsvFile(fileName + Name + '_Noisy_Result.csv',model,NoisyResult)
    if seed != 0:
        fh = open(fileName + Name + '_Seed.txt', 'wb')
        fh.write('Random Seed = ' + str(seed))
    if seed == 0:
        fh = open(fileName + Name + '_Seed.txt', 'wb')
        fh.write('Random Seed not chosen.')
    sbmlStr = model.getSBML()
    te.saveToFile (fileName + Name + '_Model.xml', sbmlStr)
    fh = open(fileName + Name + '_Antimony.txt', 'wb')
    fh.write(str(antStr))
#    plt.close('all')
    print('\nData Saved!\n')

#%%
# fileName == entire path
def writecsvFile (fileName, model, data):
    names = model.getFloatingSpeciesIds()
    header_string = names[0]
    for name in names[1-len(names):]:
        header_string = header_string + ',' + name

    fh = open(fileName, "wb")
    fh.write(header_string)
    fh.write('\n')
    for x in xrange(data.shape[0]):
        for y in xrange(data.shape[1]):
            val = data[x][y]
            s1 = "{:6.5f}".format(val)
            fh.write (s1)
            if y != data.shape[1] - 1:
                fh.write (',')
        fh.write('\n')
    fh.close()

# drawing the model with pyGraphViz
def drawOutput(model):
    ans = raw_input('Do you want to draw the model (Requires working PyGraphViz)? y/n: ')
    if ans == 'y':
        model.draw(layout='fdp')
        return()
    elif ans == 'n':
        return()
    else:
        print('Invalid entry in input.')
    return()

# Global Exception Handling
def global_exception_handling(exctype, value, TB):
    if exctype == WindowsError:
       raise
    elif exctype == ValueError:
       print("There was an issue with your last entry, check your formatting.")
       ErrorPrinting(value)
    elif Exception:
       print ("Something bad happened. Crashlog made.")
       logging.error(traceback.format_exc())
    else:
       sys.__excepthook__(exctype, value, TB)

######################################
#%%%%%%%%%% MAIN METHOD %%%%%%%%%%%%%#
######################################

####User input based program####

#for retry in range(3):
#    try:
#        initPath = os.getcwd() +'\\'
#        print('Your current directory is ' + initPath)
#        folderoption = str(raw_input('Would you like to make your model folder somewhere else? y/n: '))
#        if folderoption == 'y':
#            initPath = str(raw_input('Please write the absolute directory where you would like to  save your model folder.: '))
#            if initPath == '':
#                retry += 1
#                raise AssertionError('You didn\'t write in a directory name!\n')
#            if initPath[-1] != '\\':
#                initPath += '\\'
#            print('\nYour new directory to save in is ' + initPath)
#            break
#        elif folderoption == 'n':
#            break
#        else:
#            print('Your response must be a y or n.\n')
#            retry += 1
#            continue
#    except AssertionError:
#        print('Folder name is empty, invalid or disallowed. Try again.\n')
#        retry += 1
#    except ValueError as error:
#        print("There was an issue with your last entry, check your formatting.\n")
#        ErrorPrinting(error)
#        retry += 1
#else:
#    raise AssertionError('You have too many invalid inputs!')
#
#for retry in range(3):
#    try:
#        folderName = ""
#        folderName = str(raw_input('Please choose a folder name in which you will save generated data(Just press enter to use your current directory): '))
#        forbiddenChar = ['<','>',':','"','/','\'','|','?','*']
#        #Illegal characters in foldernames
#        assert list(set(folderName).intersection(forbiddenChar))==[]
#        if folderName == '':
#            folderPath = initPath + 'GRN generator output\\'
#            if os.path.exists(folderPath) == True:
#                break
#        elif folderName[-1]== ' ' or folderName[-1] == '.':
#            #file and folder names cannot end in space or period
#                raise AssertionError
#        else:
#            folderPath = initPath + folderName
#        if os.path.exists(folderPath) == False:
#            os.mkdir(folderPath)
#            print('\nFolder created: ' + folderPath)
#            break
#            #If the folder exists, don't do anything.
#        break
#    except WindowsError as error:
#        ErrorPrinting(error,'Folder Error: ')
#        raise WindowsError
#    except AssertionError:
#        print('Folder name is invalid or disallowed. Try again.\n')
#        retry += 1
#    except ValueError as error:
#        print("There was an issue with your last entry, check your formatting.\n")
#        ErrorPrinting(error)
#        retry += 1
#else:
#    raise AssertionError('You have too many invalid inputs!')
#
#for retry in range(3):
#    try:
#        modelName = str(raw_input("Please enter the name of your network(str): "))
#        assert list(set(modelName).intersection(forbiddenChar))==[]
#        if modelName[-1]== ' ' or modelName[-1] == '.':
#            #file and folder names cannot end in space or period
#            raise AssertionError
#        modelPath = folderPath + modelName + '\\'
#        break
#    except ValueError as error:
#        print("There was an issue with your last entry, check your formatting.\n")
#        ErrorPrinting(error)
#        retry += 1
#    except AssertionError:
#        print('Network name includes a "forbidden" character. Try again.\n')
#        retry += 1
#else:
#    raise AssertionError('You have too many invalid inputs!')
#
#def newmodelpath(count,modelPath):
#    response = str(raw_input('You have a model output with this name. Would you like to replace the files?(y/n): '))
#    if count > 3:
#        raise AssertionError
#    if response == 'n':
#        modelName = str(raw_input("Please enter a different name for your network(str): "))
#        modelPath = folderPath + modelName + '\\'
#        if os.path.exists(modelPath) == True:
#            count += 2
#            count, modelPath = newmodelpath(count,modelPath)
#            return(count,modelPath)
#        else:
#            os.mkdir(modelPath)
#            return(count,modelPath)
#    elif response == 'y':
#        return(count,modelPath)
#    else:
#        print('Your response must be a y or n.\n')
#        count += 2
#        count, modelPath = newmodelpath(count,modelPath)
#        return(count,modelPath)
#
#for retry in range(3):
#    try:
#        if os.path.exists(modelPath) == False:
#            os.mkdir(modelPath)
#            break
#        if os.path.exists(modelPath) == True:
#            inner_retry = 0
#            inner_retry,modelPath = newmodelpath(inner_retry,modelPath)
#            retry += 1
#            break
#    except AssertionError:
#        print("You keep picking existing names.")
#        retry += 2
#    except ValueError as error:
#        print("There was an issue with your last entry, check your formatting.\n")
#        ErrorPrinting(error)
#        retry += 1
#    except WindowsError as error:
#        ErrorPrinting(error,'Folder Error: ')
#        raise WindowsError
#else:
#    raise AssertionError('You have too many invalid inputs!')
#
#for retry in range(3):
#    try:
#        numGenes = int(raw_input("Please enter the number of genes in your network(2<int<=30): "))
#        assert numGenes > 2 and numGenes <= 30
#        if numGenes <= 30 and numGenes >= 16:
#                    print('\nWarning: Due to the number of genes ('+str(numGenes)+'), the legend may go out of bounds.')
#        break
#    except AssertionError:
#        if numGenes <= 2:
#            print('The number of genes must larger than 2.\n')
#        if numGenes > 30:
#            print('The number of genes must be less than or equal to 30.\n')
#        retry += 1
#    except ValueError as error:
#        print("Please enter integers only.\n")
#        ErrorPrinting(error)
#        retry += 1
#else:
#    raise AssertionError('You have too many invalid inputs!')
#
#for retry in range(3):
#    try:
#        regProb=[]
#        regProb.append(float(raw_input("Please enter the probability of Single regulation.(float): ")))
#        regProb.append(float(raw_input("Please enter the probability of Double regulation.(float): ")))
#        CountProb = round(1-np.sum(regProb),3)
#        regProb.append(CountProb)
#        assert all(i>=0 for i in regProb) == True and round(np.sum(regProb),3) == 1.0
#        print('The probabilty of Counter regulation will be ' + str(CountProb))
#        break
#    except AssertionError:
#        print('The sum of the regulation Probabilities must equal 1 and cannot be negative.\n')
#        retry += 1
#    except ValueError as error:
#        print("Please enter numbers for your regulation Probabilities.\n")
#        ErrorPrinting(error)
#        retry += 1
#else:
#    raise AssertionError('You have too many invalid inputs!')
#
#for retry in range(3):
#    try:
#        tmax = float(raw_input("Please enter how long you'd like to simulate the network for. Recommended: 30-50. (float): "))
#        assert tmax > 0
#        if tmax > 50:
#            confirm = str(raw_input('You probably only need to run it for about 70 time points. Are you sure you want to run it for '+ str(tmax) + ' time points? (y/n): '))
#            if confirm == 'y':
#                break
#            elif confirm == 'n':
#                retry += 1
#                continue
#            else:
#                print('Your response must be a y or n.\n')
#                retry += 1
#                continue
#        break
#    except AssertionError:
#        print("Your simulation time must be positive.\n")
#        retry += 1
#    except ValueError as error:
#        print("Please enter a number.\n")
#        ErrorPrinting(error)
#        retry += 1
#else:
#    raise AssertionError('You have too many invalid inputs!')
#
##    Steps = int(raw_input("Please enter number of steps to generate for your simulation output (int): "))
##    if Steps <= 0 or Steps >=50:
##        Steps = int(raw_input("You entered a number of steps that is outside of the allowable range, please enter an integer between 0 and 50 (Exclusive): "))
##        if Steps <= 0 or Steps >=50:
##            int('i')
#
#for retry in range(3):
#    try:
#        NLevel = float(raw_input("Please enter the percentage of the noise level (float): "))
#        assert NLevel > 0 and NLevel <100
#        break
#    except AssertionError:
#        print("Please enter a number between 0 and 100.\n")
#        retry += 1
#    except ValueError as error:
#        print("Please enter a number.\n")
#        ErrorPrinting(error)
#        retry += 1
#else:
#    raise AssertionError('You have too many invalid inputs!')
#
#for retry in range(3):
#    try:
#        DataOut = []
#        if NLevel < 0.01:
#            DataOut.append ('n')
#        else:
#            DataOut.append(str(raw_input("Would you like simulated data exported(y/n): ")))
#            assert DataOut[0] == 'n' or DataOut[0] =='y'
#        DataOut.append(str(raw_input("Would you like noisy simulated data exported(y/n): ")))
#        assert DataOut[1] == 'n' or DataOut[1] =='y'
#        break
#    except AssertionError:
#        print("Your response must be a y or n.\n")
#        retry += 1
#    except ValueError as error:
#        print("There was an issue with your last entry, check your formatting.")
#        ErrorPrinting(error)
#        retry += 1
#else:
#    raise AssertionError('You have too many invalid inputs!')
#
#for retry in range(3):
#    try:
#        seed = int (raw_input ('Please enter an optional random seed (enter 0 for no seed): '))
#        assert seed >= 0 and seed < 0xffffffff
#        break
#    except AssertionError:
#        print("Your response was outside of the allowable range, please enter a 32-bit unsigned integer or 0.\n")
#        retry += 1
#    except ValueError as error:
#        print("Please enter a 32-bit unsigned integer.")
#        ErrorPrinting(error)
#        retry += 1
#else:
#    raise AssertionError('You have too many invalid inputs!')
#
##### END USER INPUT ####
#
#if seed != 0:
#    np.random.seed (seed)
#NLevel /= 100
#InitParams = [0,0,0.5,0.9,0.8,30,30,0.2,0.5,1]
#
#model,antStr,GRNInt = GetModel(tries,numGenes,regProb, InitParams, modelName)
#
#tries = 0
#DataOut, model, antStr, NoisyResult, result = RunModel(tries,model,antStr,tmax,regProb,NLevel,modelName,DataOut,seed, modelPath)
#Output(DataOut, modelName, model, seed, result, NoisyResult, modelPath, antStr)
#makePlots(tmax,result,NoisyResult)
#drawOutput(model)

###########################################
#############MAIN METHOD ENDS##############
###########################################

sys.excepthook = global_exception_handling

#########################################
#%%% Automatic block ####################
#########################################

#Running code without asking user for input (put in all the following lines into the console after runnning rest of script once). Or, comment out the main method and run.

regProb=[.3,.4,.3]
NLevel=float(10)
modelName = 'test'
print 'Name of model: ' + str(modelName)
numGenes=8
InitParams = [0,0,0.5,0.9,0.8,30,30,0.2,0.5,1] # Initial values the parameters should be generated around as a median on a normal distribution
tmax = 30
seed = 3
if seed != 0:
    np.random.seed (seed)
NLevel /= 100
runAttempts = 2
filePath = 'C:\\Users\\Yoshi\\Documents\\GitHub\\DREAM-work\\Alan\'s codes'


initPath = os.getcwd() +'\\'
folderPath = initPath + 'GRN generator output\\'
modelPath = folderPath + modelName + '\\'
if os.path.exists(modelPath) == False:
    if os.path.exists(folderPath) == False:
        os.mkdir(folderPath)
        print('\nFolder created: ' + folderPath)
    os.mkdir(modelPath)
    print('\nFolder created: ' + modelPath)
DataOut = ['y','y']
print 'finish dataout'
tries = 1
model,antStr,GRNInt = GetModel(tries,numGenes,regProb, InitParams, modelName)
print 'finish GetModel'
try:
    DataOut, model, antStr, NoisyResult, result = RunModel(runAttempts-1,tries,model,antStr,tmax,regProb,NLevel,modelName,DataOut,seed, modelPath)
except ValueError as error:
    ErrorPrinting(error)
    print('Writing: ' + modelPath + modelName + '_Antimony.txt for debugging.')
    fh = open(modelPath + modelName + '_Antimony.txt', 'wb')
    fh.write(str(antStr))
else:
    makePlots(tmax,result,NoisyResult)
    Output(DataOut, modelName, model, seed, result, NoisyResult, modelPath, antStr)

#model.draw(layout='fdp')