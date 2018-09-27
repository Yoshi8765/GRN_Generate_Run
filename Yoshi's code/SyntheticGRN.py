# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 14:52:12 2018
SyntheticGRN: Creates an Antimony model of a GRN
@author: Alan
"""
import tellurium as te
import numpy as np
import math
import matplotlib.pyplot as plt
#import os

NUCLEOTIDE_SOURCE = '' #'$N'
AMINO_ACID_SOURCE = '' #'$A'

# The ModelInitNames subroutine creates variable names associated with the proper number
# of numGenes. Once generated, names are appended to the corresponding dictionary
# key of the GRN dictionary. This sub routine also constructs a probability
# matrix with a probability for any given regulatory behavior.
# This subroutine returns the constructed probability matrix.
def ModelInitNames(numGenes,GRN,regProb):
    for k in range(numGenes):
        num = str(k+1)
        GRN['Leak'].append("L" + num)
        GRN['Rate'].append('R' + num +' := ')
        GRN['mDeg'].append('Dm' + num)
        GRN['PDeg'].append('Dp' + num)
        GRN['mRNA'].append('M' + num)
        GRN['Prot'].append('P' + num)
        GRN['TRa'].append('TRa' + num)
        GRN['TRb'].append('TRb' + num)
        GRN['HillCoeff'].append('H' + num)
        GRN['Dissoc'].append('K' + num)
        GRN['DualDissoc'].append('Kd' + num)
    ProbabilityMatrix = []
    ProbSingle = float(regProb[0])
    ProbDouble = float(regProb[1])
    ProbCounter = float(regProb[2])
    Psing = ProbSingle/2
    Pdub = ProbDouble/4
    ProbabilityMatrix = [Psing, Psing, Pdub, Pdub, Pdub, Pdub, ProbCounter]
    return ProbabilityMatrix

# The Assign Proteins subroutine randomly determines the number of numGenes that will function as activators and ensures that these proteins only occur in
# the set of activators once. The remainder of the proteins are asigned as repressors.
def AssignProteins(numGenes,GRN,GRNInt):
    nums =[int(math.ceil(float(numGenes)/2)), int(math.floor(float(numGenes)/2))]
    RanAct = np.random.choice(GRN['Prot'], nums[np.random.randint(2)])
    GRNInt.update({'activators': RanAct})
    index1 = np.arange(len(GRNInt['activators']))
    index2 = np.sort(-index1)
    for m in range(len(index1)):
        ActivatorCur =GRNInt['activators'][m]
        Activators1 = GRNInt['activators'][index2[m]:]
        Activators2 =  GRNInt['activators'][:index1[m]]
        if m > len(index1)-2:
            Activators1 = []
        if ActivatorCur in Activators1 or ActivatorCur in Activators2:
            AssignProteins(numGenes,GRN,GRNInt)
    Proteins = []
    for h in range(len(GRN['Prot'])):
        Proteins.append(GRN['Prot'][h])
    for s in GRNInt['activators']:
        if s in Proteins:
            Proteins.remove(s)
        GRNInt.update({'repressors': Proteins})
    return

# The DetermineInteractors subroutine randomly generates a set of regulatory interactions equal to the number of numGenes in the network. The correct number
# and type of regulatory protein are randomly assigned to any reaction, in a fashion that ensures no protein performs dual regulation.
def DetermineInteractors(ProbabilityMatrix,Interactions,numGenes,GRNInt):
    GRNInt.update({'interactions': np.random.choice(Interactions,numGenes,True,ProbabilityMatrix)})
    for Inter in range(len(GRNInt['interactions'])):
        reg = GRNInt['interactions'][Inter]
        RegType = 'repressors'
        TypeA = 'activators'
        TypeR = 'repressors'
        if 'single' in reg :
            if 'act' in reg:
                RegType = 'activators'
            GRNInt['TF'].append([reg, [GRNInt[RegType][np.random.randint(len(GRNInt[RegType]))]]])
        else:
            Tps={'activators':[], 'repressors':[]}
            for m in Tps.keys():
#                    n=1
                Tps[m] =np.random.choice(GRNInt[m],2)
                attempts = 0
                while Tps[m][0] == Tps[m][1] and attempts<=10:
                    Tps[m][1] = np.random.choice(GRNInt[m])
                    attempts+=1
            if reg == "and" or reg == "or":
                GRNInt['TF'].append([reg, [Tps[TypeA][0],Tps[TypeA][1]]])
            elif reg == "Nand" or reg == "Nor":
                GRNInt['TF'].append([reg, [Tps[TypeR][0],Tps[TypeR][1]]])
            else:
                GRNInt['TF'].append([reg, [Tps[TypeA][0],Tps[TypeR][1]]] )
    return

# Identify proteins regulating numGenes
def ChooseProtein(numGenes,GRNInt):
    for i in range(numGenes):
        TF=[]
        string = GRNInt['TF'][i][1]
        for n in string:
            TF.append(int(n.replace("P","")))
        GRNInt['TFs'].append(TF)
    return

# Generate Rxn Rates
def GetReactionRates(numGenes,GRNInt,GRN,InducedGenes):
    ReactionRates=[]
    WorkingString=''
    for i in range(numGenes):
        Reg = GRNInt['TF'][i][0]
        WorkingString += 'Rm_' + str(i+1) + ' := '
        TF1 = GRNInt['TFs'][i][0] - 1
        P1 = GRN['Prot'][TF1]
        K1 =  GRN['Dissoc'][TF1]
        H1 = GRN['HillCoeff'][TF1]
        num = ''
        denom = ''
        frac = ''
        if i in InducedGenes:
            WorkingString = GRN['Leak'][i] + '+' + GRN['TRa'][i] +  '-' + GRN['mRNA'][i] + '*' + GRN['mDeg'][i]
        if not(i in InducedGenes) and'single' in Reg:
            frac = '(' + P1 + '^' + H1 + '*' + K1 + ')'
            if 'act' in Reg:
                Hill = '(' + frac + '/(1 + ' + frac + '))'
            else:
                Hill = '(1 /(1 + ' + frac + '))'
            WorkingString = GRN['Leak'][i] + '+' + GRN['TRa'][i] + '*' + Hill + '-' + GRN['mRNA'][i] + '*' + GRN['mDeg'][i]
        elif not(i in InducedGenes) and len(GRNInt['TFs'][i]) >1 :
            TF2 = GRNInt['TFs'][i][1] - 1
            P2 = GRN['Prot'][TF2]
            K2 =  GRN['Dissoc'][TF2]
            K3 = GRN['DualDissoc'][TF2]
            H2 = GRN['HillCoeff'][TF2]
            if 'Counter' in Reg:
                num = '(' + K1 + '*' + P1+ '^' + H1 + ')'
                denom = '(1+' + K1 + '*' + P1 + '^' + H1 + '+' + K2 + '*' + P2 + '^' + H2 + '+' + K3 + '*' + P1 + '^' + H1 + '*' + P2 + '^' + H2 + ')'
                eq = '(' + num  + '/' + denom + ')'
            elif  not('N' in Reg) and("and" in Reg or 'or' in Reg):
                if 'and' in Reg:
                    num = '(' + K1 + '*' + K2 + '*' + P1 + '^' + H1 + '*' + P2 + '^' + H2  + ')'
                    denom = '(1+' + K1 + '*' + P1 + '^' + H1 + '+' + K2 + '*' + P2 + '^' + H2 + '+' + K1 + '*' + K2 + '*' + P1 + '^' + H1 + '*' + P2 + '^' + H2 + ')'
                    eq = '(' + num  + '/' + denom + ')'
                else:
                    num = '(' + K1 + '*' + P1 + '^' + H1 + '+' + K2 + '*' + P2 + '^' + H2 + ')'
                    denom = '(1+' + K1 + '*' + P1 + '^' + H1 + '+' + K2 + '*' + P2  + '^' + H2 + ')'
                    eq = '(' + num  + '/' + denom + ')'
            else:
                denom = '(1+' + K1 + '*' + P1 + '^' + H1 + '+' + K2 + '*' + P2 + '^' + H2 + '+' + K3+ '*' + P1 + '^' + H1 + '*' + P2 + '^' + H2 + ')'
                if 'and' in Reg:
                    num = '(1+' + K1 + '*' + P1 + '^' + H1 + '+' + K2 + '*' + P2 + '^' + H2 + ')'
                    eq = '(' + num  + '/' + denom + ')'
                else:
                    eq =  '(1' + '/' +  denom + ')'
            WorkingString = GRN['Leak'][i] + '+' + GRN['TRa'][i] + '*' + eq + '-' + GRN['mRNA'][i] + '*' + GRN['mDeg'][i]
        ReactionRates.append(WorkingString)
    return ReactionRates

# Reaction Model Generation
def GetReactions(numGenes,GRN,Rates,StringList):
    ReactionM =''
    ReactionP =''
    ReactionSumString= 'model Random_GRN()\n\n'
    for i in np.arange(1, numGenes+1):
        M = GRN['mRNA'][i-1]
        P = GRN['Prot'][i-1]
        ReactionM += 'Rm' + str(i) + ':' + NUCLEOTIDE_SOURCE + '=> ' + M + ';' + Rates[i-1] + '\n'
        ReactionP += 'Rp' + str(i) + ':' + AMINO_ACID_SOURCE + '=> ' + P + ';' + GRN['TRb'][i-1] + '*' + M + '-' + P + '*' + GRN['PDeg'][i-1] + '\n'
    ReactionSumString += ReactionM + '\n' + ReactionP + '\n'
    StringList.append(ReactionSumString)

# Comment?
def ValueGeneration(paramRange,numGenes,InitVals,StringList):
    Species=''
    Values = {'M':[], 'P':[],'Dp':[],'Dm':[],'L':[],'TRa':[], 'TRb':[],'Kd':[], 'K':[], 'H': []}
    Keys = Values.keys()
    if len(paramRange) != 0:
        for i in Keys:
            for n in np.arange(1, numGenes+1):
                val = 0
                if i == 'M' or i=='P':
                    Values[i].append(val)
                elif i == 'H':
                    while val < 1 or val > 8:
                        val = round(np.random.uniform(paramRange[2],paramRange[3]),0)
                        Values[i].append(val)
                else:
                    count = 0
                    while val <= 0 and count <= 10: #bug: infinite loop here, fixed for now with count.
                        # TODO: fix this
                        val = round(np.random.uniform(paramRange[0], paramRange[1]),3)
                        Values[i].append(val)
                        count+=1
                Species += i + str(n) + ' = ' + str(val) + '; \n'
            Species+='\n'
    else:
        Values = {'M':[InitVals[0]], 'P':[InitVals[1]],'Dp':[InitVals[2]],'Dm':[InitVals[3]],'L':[InitVals[4]],'TRa':[InitVals[5]], 'TRb':[InitVals[6]],'Kd':[InitVals[7]], 'K':[InitVals[8]], 'H': [InitVals[9]]}
        for i in Keys:
            for n in np.arange(1, numGenes+1):
                val = 0
                if i == 'M' or i=='P':
                    Values[i].append(val)
                elif i == 'H':
                    while val < 1 or val > 8:
                       val = round(np.random.normal(Values[i][0],3),0)
                       Values[i].append(val)
                else:
                    while val <= 0:
                        val = round(np.random.normal(Values[i][0],.25),3)
                        Values[i].append(val)
                Species += i + str(n) + ' = ' + str(val) + '; \n'
            Species+='\n'
    StringList.append(Species)
    return

def GetModelString(StringList):
    antStr = ''
    for i in StringList:
        antStr += i
    antStr += 'end'
    return antStr

# Reporting a summary of what regulations are included in the model.
def GetRegulations(numGenes,GRNInt,InducedGenes):
    Regulations='Regulatory behavior: \n\n'
    for i in range(numGenes):
        Reg = GRNInt['TF'][i][0]
        if i in InducedGenes:
            Regulations+= 'Gene ' + str(i+1) + ' was induced.\n'
        else:
            for k in range(len(GRNInt['TFs'][i])):
                if k == 0:
                    Regulations += 'Gene ' + str(i+1) + ' is regulated by P' + str(GRNInt['TFs'][i][k])
                elif k == 1:
                    Regulations += ' and P' + str(GRNInt['TFs'][i][k])
            Regulations += ' via ' + Reg + ' regulation. \n'
    return Regulations


# The RandomGRN prepare random genetic regulatory networks. The routine takes 3 arguments
# 1. an integer equal to the number of numGenes in the network. 2. a list of 3 floats that are equal
# indicate the probability of single regulation, dual repression(activator or repression), and counter
# regulation(both activation and repression) respectively. 3. .
def RandomGRN(numGenes,regProb,paramRange,InducedGenes):
    InitVals = [0,0,0.5,0.9,0.8,30,30,0.2,0.5,1]
    Interactions = ["single activation", "single repression", "and", "or", "Nand", "Nor", "Counter"]
    GRN = {'Rate': [], 'mRNA':[], 'Prot': [],'PDeg':[],'mDeg': [],'Leak': [],'TRa':[], 'TRb':[], 'DualDissoc':[], 'Dissoc':[], 'HillCoeff': []}
    GRNInt = {'TF':[],'TFs':[]}
    StringList=[]

    ProbabilityMatrix = ModelInitNames(numGenes,GRN,regProb)
    AssignProteins(numGenes,GRN,GRNInt)
    DetermineInteractors(ProbabilityMatrix,Interactions,numGenes,GRNInt)
    ChooseProtein(numGenes,GRNInt)
    Rates = GetReactionRates(numGenes,GRNInt,GRN,InducedGenes)
    GetReactions(numGenes,GRN,Rates,StringList)
    ValueGeneration(paramRange,numGenes,InitVals,StringList)
    antStr = GetModelString(StringList)
    Regulations = GetRegulations(numGenes,GRNInt,InducedGenes)
    return((Regulations,antStr))

# Comment?
def Output(modelName, folderPath, rSeed, r, antStr, Regulations):
    if rSeed != 0:
        fh = open(folderPath + modelName + ' Seed.txt', 'wb')
        fh.write('Random Seed = ' + str(rSeed))
    sbmlStr = r.getSBML()
    te.saveToFile (folderPath + modelName + ' Model.xml', sbmlStr)
    fh = open(folderPath + modelName + ' Antimony.txt', 'wb')
    fh.write(str(antStr))

    fh = open(folderPath + modelName + ' Network Regulations.txt', 'wb')
    fh.write(Regulations)
    plt.close('all')

# -----------------------------------------------------------
# This is the function a user should call if creating a model
# Comment?
def GetModel(numGenes, modelName, folderPath, regProb, rSeed = 0, paramRange=[], InducedGenes=[],MAX_TRIES=100):
    if len(InducedGenes)>0:
        for i in range(len(InducedGenes)):
            InducedGenes[i]-=1
    if rSeed != 0:
        np.random.seed (rSeed)
#    Percent /= 100
    tries = 0
    antStr = ''
    while tries < MAX_TRIES:
       try:
           Regulations,antStr=RandomGRN(numGenes,regProb, paramRange,InducedGenes)
#          print "AntString = ", antStr
           r = te.loadAntimonyModel(antStr)
           r.steadyState()
           print('Success! A random GRN was generated after ' + str(tries + 1) + ' tries')
           Output(modelName, folderPath, rSeed, r, antStr, Regulations)
           break
       except:
           tries = tries + 1
           print "Failed to find steady state, generating a new GRN"
           Output(modelName, folderPath, rSeed, r, antStr, Regulations)
           return antStr
    return antStr
