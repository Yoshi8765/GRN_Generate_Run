# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 09:48:51 2018

@author: Alan
"""
import tellurium as te
import numpy as np
import math
import matplotlib.pyplot as plt
import os
#import SyntheticGRN as Syn


def writecsvFile (folderPath, r, data):
    names = r.getFloatingSpeciesIds()
    header_string = 'time,' + names[0] 
    for name in names[1-len(names):]:
        header_string = header_string + ',' + name 
    fh = open(folderPath, "wb")
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
    
    
def Output(DataExport, modelName,r,result,Percent, NoisyResult, folderPath):
    for i in range(len(DataExport)):
        if i == 0 and DataExport[i] == 'y':
            writecsvFile(folderPath + modelName + ' Results.csv',r,result)
        elif DataExport[i] == 'y' :
            writecsvFile(folderPath + modelName + ' Noisy Result.csv',r,NoisyResult)
    
# ---------------------------------------------------------------------------------
# This is the function a user should call  
            
# RunModel takes the number of numberOfGenes, the duration of the simulation, a list of three probabilities for regulation types, the noise level, 
# the name of the simulations, the user defined responses that reflect which data will be output, a seed for reproducibile simulations,
# the number of steps for the simulation, and a string associated with the directory the main output directory for all GRN generations.
#
# Given the correct inputs, this sub routine will call and sub routine GetModel, which will result in the creation of a random GRN(if no seed input)
# The model is simulated using user defined inputs(if the model fails a new GRN is generated up to 10 times). A noisy dataset is created by 
# adding noise to the simulated dataset. Figures containing Protein and mRNA Vs time plots are created and saved for regular and noisy simulations.
# Finally all the Outputs sub-routine is called which w
def RunModel(numberOfGenes, modelName,  filePath, rSeed, NLevel, tMax, Steps, DataExport, ModelString ='', ExternalModel = ''):
    if ModelString =='' and os.path.isfile(ExternalModel):
        AntImport = open(ExternalModel, 'r')
        ModelString = AntImport.read()
    if rSeed != 0:
        np.random.seed (rSeed)
    NLevel /= 100
#    Values = {'M':[InitVals[0]], 'P':[InitVals[1]],'Dp':[InitVals[2]],'Dm':[InitVals[3]],'L':[InitVals[4]],'TRa':[InitVals[5]], 'TRb':[InitVals[6]],'Kd':[InitVals[7]], 'K':[InitVals[8]], 'H': [InitVals[9]]}
#    InitParams = [0,0,0.5,0.9,0.8,30,30,0.2,0.5,1]
    r = te.loadAntimonyModel(ModelString)
#    Regulations,r,antStr = Syn.GetModel(numberOfGenes, InitProb, InitParams) 
    tries = 0
    while tries < 10:
        try:
            plt.close('all')
            r.reset()
            result = r.simulate(0,tMax, Steps)
            NoisyResult = np.zeros([len(result[:,0]), len(result[0,:])])        
            if numberOfGenes<11:
                vars = {'P':[],'M':[]}
                for e in vars.keys():
                    plt.figure(1)    
                    for k in np.arange(1,numberOfGenes+1):
                        vars[e].append(e + str(k))
                        tStep = int(math.ceil(tMax)/20)
                        if tMax < 20:
                            tStep = int(math.ceil(tMax)/(tMax/2))
                        tRange = np.arange(0,(int(math.ceil(tMax)))+ tStep,tStep)
                        if e=='P':
                            plt.subplot(2,1,1)
                            plt.title('Protein Vs. Time', fontsize=8)
                            plt.grid(color='k', linestyle='-', linewidth=.4)
                            plt.ylabel('count',fontsize=6)
                            plt.yticks(fontsize=6)
                            if tMax > 999: 
                                plt.loglog(result[:,0],result[:,(k-1)+1], label = vars[e][k-1])
                                plt.xticks(fontsize=6)
                            else:
                                plt.semilogy(result[:,0],result[:,(k-1)+1], label = vars[e][k-1])
                                plt.xticks(tRange, fontsize=6) 
                            plt.xlim(0,tMax)
                            plt.legend(loc='upper right', bbox_to_anchor=(1.13, 1.035), fontsize=6)
                        else:
                           plt.subplot(2,1,2) 
                           plt.title('mRNA Vs. Time', fontsize=8)
                           plt.grid(color='k', linestyle='-', linewidth=.4)
                           plt.xlabel('time(s)',fontsize=6)
                           plt.ylabel('count',fontsize=6)
                           plt.yticks(fontsize=6)
                           if tMax > 999: 
                               plt.loglog(result[:,0],result[:,(k-1)+numberOfGenes+1], label = vars[e][k-1])
                               plt.xticks(fontsize=6)
                           else:
                               plt.semilogy(result[:,0],result[:,(k-1)+numberOfGenes+1], label = vars[e][k-1])
                               plt.xticks(tRange, fontsize=6) 
                           plt.xlim(0,tMax)
                           plt.legend(vars[e],loc='upper right', bbox_to_anchor=(1.13, 1.035), fontsize=6)
                manager = plt.get_current_fig_manager()
                manager.window.showMaximized()
                plt.savefig(filePath + modelName + ' Simulation Plot.png', dpi=400)             
                plt.close('all')
            if NLevel != 0:
                for k in range(len(result[:,0])):
                    for i in range(len(result[0])):
                        if i == 0:
                            NoisyResult[k,i] = result[k,i]
                        else:    
                            CurrVal = -1
                            while CurrVal < 0:
                                CurrVal = result[k,i] + np.random.normal(0,NLevel*result[k,i])
                            NoisyResult[k,i] = CurrVal
                if numberOfGenes<11:
                    for e in vars.keys():
                        plt.figure(2)    
                        for k in np.arange(1,numberOfGenes+1):
                            vars[e].append(e + str(k))
                            if e=='P':
                                plt.subplot(2,1,1)
                                plt.title('Protein Vs. Time (Noisy)', fontsize=8)
                                plt.grid(color='k', linestyle='-', linewidth=.4)
                                plt.ylabel('count',fontsize=6)
                                plt.yticks(fontsize=6)
                                plt.legend(vars[e])
                                if tMax > 999: 
                                    plt.loglog(NoisyResult[:,0],NoisyResult[:,(k-1)+1], label = vars[e][k-1])
                                    plt.xticks(fontsize=6)
                                else:
                                   plt.semilogy(NoisyResult[:,0],NoisyResult[:,(k-1)+1], label = vars[e][k-1]) 
                                   plt.xticks(tRange, fontsize=6) 
                                plt.xlim(0,tMax)
                                plt.legend(loc='upper right', bbox_to_anchor=(1.13, 1.035), fontsize=6)
                            else:
                                plt.subplot(2,1,2)
                                plt.title('mRNA Vs. Time (Noisy)', fontsize=8)
                                plt.grid(color='k', linestyle='-', linewidth=.4)
                                
                                plt.yticks(fontsize=6)
                                plt.xlabel('time(s)',fontsize=6)
                                plt.ylabel('count',fontsize=6)
                                plt.legend(vars[e])
                                if tMax > 999:
                                    plt.loglog(NoisyResult[:,0],NoisyResult[:,(k-1)+numberOfGenes+1], label = vars[e][k-1])
                                    plt.xticks(fontsize=6)
                                else:
                                    plt.semilogy(NoisyResult[:,0],NoisyResult[:,(k-1)+numberOfGenes+1], label = vars[e][k-1])
                                    plt.xticks(tRange, fontsize=6)  
                                plt.xlim(0,tMax)
                                plt.legend(loc='upper right', bbox_to_anchor=(1.13, 1.035), fontsize=6)
                    manager = plt.get_current_fig_manager()
                    manager.window.showMaximized()
                    plt.savefig(filePath + modelName + ' Noisy Simulation Plot.png', dpi=400)  
                    plt.close('all')
            break
        except:
             tries = tries + 1
             print "Failed to simulate your model, generating a new GRN."
             if tries ==10:
                 print('Sorry you network was not able to be simulated, please try again, or read the documentation')
    Output(DataExport, modelName, r, result, NLevel, NoisyResult, filePath)