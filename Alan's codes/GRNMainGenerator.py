# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 09:49:46 2018

@author: Alan
"""
import pandas as pd
import csv
import numpy as np
import os
import matplotlib.pyplot as plt
import SyntheticGRN as Syn
import SimulateGRN as Sim
#User response limit variable initialization

MAX_GENES = 95
MIN_GENES = 3
MAX_TIME = 1E5
MIN_TIME = 10  
MAX_STEPS = 5000
MIN_STEPS = 0
MAX_NOISE = 100
MIN_NOISE = 0
#InitParams = [0,0,0.5,0.9,0.8,30,30,0.2,0.5,1]
ParamUpper = 1000
ParamLower = 0.1
HillUpper = 8
HillLower = 1
# User input
#try:
#Generate = str(raw_input('Would you like to generate a new GRN))
#    Simulate = str(raw_input('Would you like to simulate an Antimony Model))
#    initPath = os.getcwd()
#    print('Your current directory is ' + initPath)
#    folderName = str(raw_input('Please choose a folder for saving your generated data(press enter for default): '))
#    forbiddenChar = ['<','>',':','"','/','\'','|','?','*']
#    if any(i in folderName for i in forbiddenChar):
#          folderName = str(raw_input('You seem to have included a "forbidden" character in your folder name, please try another, or simply press enter for the default: '))      
#          if any(i in folderName for i in forbiddenChar):
#              int('i')
#    folderPath = initPath + '\\' + folderName + '\\' 
#    if folderName == '':
#        folderPath = initPath + '\\' + 'GRN generator output\\'
#        if os.path.exists(folderPath) == False:
#            os.mkdir(folderPath)
#    else:
#        if os.path.exists(folderPath) == False:
#            os.mkdir(folderPath)
#    modelName = str(raw_input("Please enter the name of your network(str): ")) 
#    if any(i in modelName for i in forbiddenChar):
#          modelName = str(raw_input('You seem to have included a "forbidden" character in your folder name, please try another, or simply press enter for the default: '))      
#          if any(i in modelName for i in forbiddenChar):
#              int('i')
#    filePath = folderPath + modelName + '\\'
#    if os.path.exists(filePath) == True:
#        response = str(raw_input('You have a model ouput with this name. Would you like to replace the files?(y/n): '))
#        if response != 'n' and response!='y':
#            response = str(raw_input('Your response must be a y or n. Please try again.(y/n): '))
#            if response != 'n' and response!='y':
#                int('i')
#        if response == 'n':
#            modelName = str(raw_input("Please enter a different name for your network(str): "))
#            filePath = folderPath + modelName + '\\'
#            os.mkdir(filePath)
#    elif os.path.exists(filePath) == False: 
#        os.mkdir(filePath)
#    print('Your files will be at this path: ' + filePath)
#    genes = int(raw_input('Please enter the number of genes in your network(' + str(MIN_GENES) + '<int<' + str(MAX_GENES) + '): '))
#    if genes < MIN_GENES:
#        print('You seem to have entered two, or less genes.\n')
#        genes = int(raw_input("Please enter the number of genes in your network(int): "))
#        if genes < MIN_GENES:
#            print('The number of genes must larger than two.\n')
#            int('i')
#    if genes > MAX_GENES:
#        print('You seem to have entered a larger number of genes than the program can handle.\n')
#        genes = int(raw_input('Please enter the number of genes in your network(' + str(MIN_GENES) + '<int<' + str(MAX_GENES) + '): '))
#        if genes > MAX_GENES:
#            int('i')
#    InitProb=[]
#    InitProb.append(float(raw_input("Please enter the probability of Single regulation.(float): ")))
#    InitProb.append(float(raw_input("Please enter the probability of Double regulation.(float): ")))
#    CountProb = round(1-np.sum(InitProb),3)
#    print('The probabilty of Counter regulation will be ' + str(CountProb) + '\n' )
#    InitProb.append(CountProb)
#    tries = 0
#    while any(i<0 for i in InitProb) == True or round(np.sum(InitProb),3) !=1.0:
#        if tries==2:
#            int('i')
#        else:   
#            print("Sorry, the sum of the probabilities must equal 1 or at least one is negative.")
#            InitProb[0] = float(raw_input("Please enter the probability (<=1) of having single regulated gene (float): "))
#            InitProb[1] = float(raw_input("Please enter the probability (<=1) of double regulation.(float): "))
#            CountProb = round(1-np.sum(InitProb[0:2]),3)
#            InitProb[2] = CountProb
#            print('\n Your new probabilty of counter regulation is ' + str(CountProb) + '\n' )
#            tries = tries + 1
#    tMax = float(raw_input('Please enter how long would you like to simulate the network for(' + str(MIN_TIME) +'<=time<=' + str(MAX_TIME) +'): '))
#    if tMax <MIN_TIME:
#        tMax = float(raw_input('Your simulation time must be larger than or equal to '+ str(MIN_TIME) + ', please enter an allowable value.'))
#    elif tMax > MAX_TIME:
#        tMax = float(raw_input('Your simulation time must be smaller than or equal to '+ str(MAX_TIME) + ', please enter an allowable value.'))
#    Steps = int(raw_input('Please enter number of steps to generate for your simulation output (' + str(MIN_STEPS) + '<steps<' + str(MAX_STEPS) + '): '))
#    if Steps <= MIN_STEPS or Steps >=MAX_STEPS:
#        Steps = int(raw_input('You entered a number of steps that is outside of the allowable range, please enter an integer between ' + str(MIN_STEPS) + ' and ' + str(MAX_STEPS) + ': '))
#        if Steps <= MIN_STEPS or Steps >=MAX_STEPS:
#            int('i')
#    
#    ParamRange=[]
#    response = str(raw_input('Would you like the parameters randomized?(y/n): '))
#    if response!= 'y' and response!= 'n':
#        response == str(raw_input('Your response must be a y or n. Please try again.(y/n)'))
#        if response!= 'y' and response!= 'n':
#            int('i')
#    elif response == 'y':
#        ParamRange.append(float(raw_input('Please enter the parameter lower limit.(limit>' + str(ParamLower) + '): ')))
#        ParamRange.append(float(raw_input('Please enter the parameter upper limit.(limit<' + str(ParamUpper) + '): ')))
#        ParamRange.append(float(raw_input('Please enter the Hill coefficient lower limit.(limit>=' + str(HillLower) + '): ')))
#        ParamRange.append(float(raw_input('Please enter the Hill coefficient upper limit.(limit<' + str(HillUpper) + '): ')))
#        if (ParamRange[1] - ParamRange[0]) < ParamLower or ParamRange[0] < ParamLower or ParamRange[1] > ParamUpper:
#            print('The parameters values chosen are out of the allowable range')
#            ParamRange[0] = (float(raw_input('Please enter the parameter lower limit.(limit>' + str(ParamLower) + '): ')))
#            ParamRange[1] = (float(raw_input('Please enter the parameter upper limit.(limit<' + str(ParamUpper) + '): ')))
#            if (ParamRange[1] - ParamRange[0]) < ParamLower or ParamRange[0] < ParamLower or ParamRange[1] > ParamUpper:
#                int('i')
#        elif (ParamRange[3] - ParamRange[2]) < HillLower or ParamRange[2] < HillLower or ParamRange[3] > HillUpper:
#            print('The Hill coefficient values chosen are out of the allowable range')
#            ParamRange[2] = (float(raw_input('Please enter the Hill coefficient lower limit.(limit>' + str(HillLower) + '): ')))
#            ParamRange[3] = (float(raw_input('Please enter the Hill coefficient upper limit.(limit<' + str(HillUpper) + '): ')))
#            if(ParamRange[3] - ParamRange[2]) < HillLower or ParamRange[2] < HillLower or ParamRange[3] > HillUpper:
#                int('i')
#    NLevel = float(raw_input("Please enter the percentage of the noise level (float): "))
#    if NLevel < MIN_NOISE or NLevel >MAX_NOISE:
#        NLevel = float(raw_input('You entered a number of steps that is outside of the allowable range, please enter a number between ' + str(MIN_NOISE) + ' and ' + str(MAX_NOISE) + ': '))
#        if NLevel < MIN_NOISE or NLevel >MAX_NOISE:
#            int('i')
#    DataExport = []
#    DataExport.append(str(raw_input("Would you like simulated data exported(y/n): ")))
#    if DataExport[0] != 'n' and DataExport[0]!='y':
#            DataExport[0] = str(raw_input('Your response must be a y or n. Please try again.(y/n)'))
#            if DataExport[0] != 'n' and DataExport[0]!='y':
#                int('i')
#    if NLevel < 0.01:
#        DataExport.append ('n')
#    else:
#        DataExport.append(str(raw_input("Would you like noisy simulated data exported(y/n): ")))
#        if DataExport[0] != 'n' and DataExport[0]!='y':
#            DataExport[0] = str(raw_input('Your response must be a y or n. Please try again.(y/n)'))
#            if DataExport[0] != 'n' and DataExport[0]!='y':
#                int('i')    
#    rSeed = int (raw_input ('Please enter an optional random seed (enter 0 for no seed): '))
#    if rSeed <0 or rSeed > 2**32 - 1:
#        rSeed = int (raw_input ('Your response was outside of the allowable range, please enter an integer between 0 and 2**32 - 1: '))
#        if rSeed <0 or rSeed > 2**32 - 1:
#            int('i')
    
#antStr = Syn.GetModel(genes, modelName, filePath, NLevel, rSeed, InitProb, InitParams, ParamRange , MAX_TRIES)
#    Sim.RunModel(genes,tMax,InitProb,NLevel,modelName,DataExport,rSeed, Steps, filePath, ModelString = antStr)
#
#except WindowsError:
#    print('folder error') 
#        
#except ValueError:
#    print("There was an issue with your last entry, check your formatting.")
#except:
#    print ("Somethng bad happenned")

#Largest network complete 80 genes
genes = 10
tMax = 500
Steps = 500
#AntimonyString = Syn.GetModel(60, 'Example', 'c:\\GRN generator output\\Example\\',  [1,0,0] , rSeed = 3456, ParamRange = [.01,100,1,8], MAX_TRIES = 10000,  InducedGenes = [3]) 
#AntimonyString = Syn.GetModel(10, 'Example', 'c:\\GRN generator output\\Example\\',  [0,1,0] , rSeed = 3456, ParamRange = [.01,100,1,8], MAX_TRIES = 10) 
AntimonyString = Syn.GetModel(50, 'Example', 'c:\\GRN generator output\\Example\\',  [1,0,0], MAX_TRIES=1000) 
#Sim.RunModel(20,10,[.4,.4,.2], 10.0,'foo0',['y','y'],3456,500,'c:\\GRN generator output\\foo0\\', ExternalModel = 'c:\\GRN generator output\\foo0\\foo0 Antimony.txt')
#Sim.RunModel(genes,'foo00','c:\\GRN generator output\\foo0\\',3456, [.4,.4,.2], 1.0, tMax,Steps,['y','y'], ModelString = antStr)

#Sim.RunModel(genes,'Example','c:\\GRN generator output\\Example\\',3456, 1.0, 500, 500, ['y','y'], ModelString = AntimonyString)
