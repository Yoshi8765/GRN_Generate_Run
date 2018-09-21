# -*- coding: utf-8 -*-
#     Author: Yoshi Goto

import os
import tellurium as te
import numpy as np
import math
from matplotlib import pyplot as plt
import imp

GRAPH_LABEL_FONTSIZE = 8
GRAPH_TITLE_FONTSIZE = 10

# TODO: perturbations?
# TODO: Masspec or rnaseq with all species
# RNASeq
# TODO: Plots do not plot anything above around 50 correctly. Solve this or else fixing plots are useless.
# TDOO: Make a table of appropriate ranges for parameters.

def run_model(antStr,noiseLevel,exportData=[ [0],'P',True,True,True],inputData=[1,40,1],bioTap='',
              savePath='\\model_output\\',showTimePlots=False,seed=0,drawModel=[False,'fdp'],runAttempts=5):
    """Checks if Antimony models will reach steady-state, generates visualizations, and exports data.

    Arguments:
        antStr: (str) The Antimony model.
        noiseLevel: (float) level of noise as a decimal. Ex: 0.05 = 5%
        exportData:
            [genesToExport,speciesDataType,csv,biotapestryCSV,sbml,antimony,seed]
            genesToExport = (int List) Default: [0] . A list object of which genes you want to export, 1-based indexing. Pass in [0] to export all proteins.
            speciesDataType = (str) Default: 'P' . 'P'rotein or 'M'rNA
            csv = (bool) Default:True . This is a flag for the export of the timecourse data of the genesToExport as csv files.
            sbml = (bool) Default:True . This is a flag for the export of the model in SBML format (version based on Tellurium).
            antimony = (bool) Default:True . This is a flag for the export of the Antimony model as a txt.

        inputData:
            [valInput,maxTime,resolution]
        bioTap: (str) If this is not empty, a csv file to use for BioTapestry will be exported.
        savePath:
        showTimePlots:
        seed:
        drawModel: Requires pygraphviz and graphviz installed
        runAttempts:



    Outputs:
        model:
        result:
        resultNoisy:
    """

    # Attempt to execute run_model up to runAttempts times with different generated models.
    for retry in range(runAttempts):
        # Load the Antimony string as a model
        try:
            model = te.loadAntimonyModel(antStr)
            print('Successfully loaded in Antimony string as a model.')
        except Exception as error:
            ErrorPrinting(error)
            raise RuntimeError("Failed to load antimony string due to parsing error. Check that your antStr is correct.")
        # Check that the model reaches steady-state. If the model has any issues running, it will raise an Exception. The user should run get_model again.
        plt.close('all')
        model.resetToOrigin()
        tStep = int(math.ceil(inputData[1]/inputData[2]))
        try:
            model.steadyState()
        except Exception:
            try:
                model.conservedMoietyAnalysis = True
                model.steadyState()
            except Exception as error:
                print "Failed to reach steady-state or there is an Integrator error. Check your model or run get_model again."
                ErrorPrinting(error)
                raise error

        model.resetToOrigin()

        #Specify input
        model.INPUT = inputData[0];
        # Run a simulation for time-course data
        result = model.simulate(0,inputData[1],tStep+1)


#            numGenes = model.getNumFloatingSpecies()
        model_name = model.getInfo().split("'modelName' : ")[1].split("\n")[0]

        # Specify (and make if necessary) a folder to save outputs to
        folderPath = os.getcwd() + '\\Random_GRNs\\'
        filesPath = folderPath + model_name + '\\'
        if os.path.exists(filesPath) == False:
            if os.path.exists(folderPath) == False:
                os.mkdir(folderPath)
            os.mkdir(filesPath)
            print('\nFolder created: ' + filesPath)

        # Create results with artificial noise
        if noiseLevel != '0':
            resultNoisy = np.zeros([len(result[:,0]), len(result[0,:])])
            for k in range(len(result[:,0])):
                for i in range(len(result[0])):
                    if i == 0:
                        resultNoisy[k,i] = result[k,i]
                    else:
                        CurrVal = -1
                        while CurrVal < 0:
                            CurrVal = result[k,i] + np.random.normal(0,noiseLevel*result[k,i])
                        resultNoisy[k,i] = CurrVal


        # Create graphs
        if showTimePlots==True:
            makePlots(exportData,inputData,filesPath,result,noiseLevel,resultNoisy)

        # Export datasets
        Output(exportData,model,seed,result,noiseLevel,resultNoisy, filesPath, antStr,bioTap)

        # Draw the model (requires pygraphviz module)
        if drawModel[0]==True:
            try:
                imp.find_module('pygraphviz')
                model.draw(layout=str(drawModel[1]))
            except ImportError:
                print('pygraphviz is not installed!')

        #returns the model, and result and/or resultNoisy arrays
        return(model,result,resultNoisy)

    raise RuntimeError('Could not create a working model for ' + str(runAttempts) + ' tries.')


###### Other Functions that GetModel uses ######

def Output(exportData,model,seed,result,noiseLevel,resultNoisy,filesPath, antStr,bioTap):
    # export csv of results
    if exportData[2]==True:
        outputSpecies = range(int(np.size(result,1)))
        if exportData[0]!=0:
            if exportData[1]=='P':
                outputSpecies = [gene*2 for gene in exportData[0]]
            if exportData[1]=='M':
                outputSpecies  = [gene*2 for gene in exportData[0]]
        else:
            if exportData[1]=='P':
                outputSpecies = outputSpecies[2::2]
            if exportData[1]=='M':
                outputSpecies = outputSpecies[1::2]
        outputResults = np.array(result[:,outputSpecies])
        outputTime = np.reshape(result[:,0],[len(result[:,0]),1])
        writeResult = np.hstack((outputTime,outputResults))
        writecsvFile(outputSpecies,filesPath + 'Results_Clean.csv',model,np.array(writeResult))
        if noiseLevel != 0:
            outputResultsNoisy = np.array(resultNoisy[:,outputSpecies])
            writeNoisy = np.hstack((outputTime,outputResultsNoisy))
            writecsvFile(outputSpecies,filesPath + 'Noisy_Result.csv',model,np.array(writeNoisy))
    # export csv file for importing into Biotapestry
    if bioTap!='':
        f2 = open(filesPath + "biotapestry.csv", 'w')
        f2.write(bioTap)
        f2.close()
    # export SBML model text
    if exportData[3]==True:
        sbmlStr = model.getSBML()
        te.saveToFile (filesPath + 'OrigModel.xml', sbmlStr)
    # export Antimony model text
    if exportData[4]==True:
        fh = open(filesPath + 'OrigAntimony.txt', 'wb')
        fh.write(str(antStr))

    print('\nData Saved!\n')

# function for exporting to csv
def writecsvFile (outputSpecies,filesPath, model, data):
#    names = model.getFloatingSpeciesIds()
#    if exportData[0]!=0:
#        if exportData[1]=='P':
#            outputSpecies = [gene for gene in exportData[0]]
#        if exportData[1]=='M':
#            exportData = [x.__add__(1) for x in exportData[0]]
#            outputSpecies = [gene for gene in exportData[0]]
    outputSpecies = [i-1 for i in outputSpecies]
    names = []
    for i in outputSpecies:
        names.append(model.getFloatingSpeciesIds()[i] + ',')
    timeStr = list(['time,'])
    #timeStr.extend(names) #fix
    header_string = timeStr + names
    header_string = "".join("".join(header_string))
#    for name in names[1-len(names):]:
#        header_string = header_string + ',' + name

    fh = open(filesPath, "wb")
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

def makePlots(exportData,inputData,filesPath,result,noiseLevel,resultNoisy):
    global GRAPH_LABEL_FONTSIZE
    global GRAPH_TITLE_FONTSIZE

    outputSpecies = range(int(np.size(result,1)))
    if exportData[0]!=0:
        outputSpecies = [gene for gene in exportData[0]]
    else:
        outputSpecies = outputSpecies[2::2]

    # create a plot for clean data
    jet = plt.get_cmap('jet')
    pickColorSpace = [int(np.floor(i)) for i in np.linspace(0,255,len(outputSpecies))]
    colors = jet(pickColorSpace)
    colors = jet(np.linspace(0,1,len(outputSpecies)))

    vars = {'P':[],'M':[]}
    for e in vars.keys():
        plt.figure(1)
        for k in np.arange(0,len(outputSpecies)):
            vars[e].append(e + str(k))
            tStep = int(math.ceil(inputData[1]/inputData[2]))
            tRange = np.arange(0,(int(math.ceil(inputData[1])))+ tStep,tStep)
            if exportData[1]=='P':
                plt.subplot(2,1,1)
                plt.title('Protein Count Vs. Time', fontsize=GRAPH_TITLE_FONTSIZE)
                plt.grid(color='k', linestyle='-', linewidth=.4)
                plt.ylabel('count',fontsize=GRAPH_LABEL_FONTSIZE)
                plt.yticks(fontsize=GRAPH_LABEL_FONTSIZE)
                plt.plot (result[:,0],result[:,outputSpecies[k]], label = vars[e][k],color=colors[k])
                #plt.yscale('log')
                if inputData[1] > 5000:
                    plt.xscale('log')
                    plt.xticks(fontsize=GRAPH_LABEL_FONTSIZE)
                else:
                    plt.xticks(tRange, fontsize=GRAPH_LABEL_FONTSIZE)
                plt.xlim(0,inputData[1])
                plt.legend(loc='upper right',ncol=3, bbox_to_anchor=(1.13, 1.035), fontsize=GRAPH_LABEL_FONTSIZE)
            if exportData[1]=='M':
               plt.subplot(2,1,2)
               plt.title('mRNA Count Vs. Time', fontsize=GRAPH_TITLE_FONTSIZE)
               plt.grid(color='k', linestyle='-', linewidth=.4)
               plt.xlabel('time(s)',fontsize=GRAPH_LABEL_FONTSIZE)
               plt.ylabel('count',fontsize=GRAPH_LABEL_FONTSIZE)
               plt.yticks(fontsize=GRAPH_LABEL_FONTSIZE)
               plt.plot (result[:,0],result[:,outputSpecies[k]+1])#, label = vars[e][k],color=colors[k])
              # plt.yscale('log')
               if inputData[1] > 5000:
                    plt.xscale('log')
                    plt.xticks(fontsize=GRAPH_LABEL_FONTSIZE)
               else:
                   plt.xticks(tRange, fontsize=GRAPH_LABEL_FONTSIZE)
               plt.xlim(0,inputData[1])
               plt.legend(vars[e],loc='upper right',ncol=3, bbox_to_anchor=(1.13, 1.035), fontsize=GRAPH_LABEL_FONTSIZE)
    manager = plt.get_current_fig_manager()
    manager.window.showMaximized()
    manager.window.showMinimized()
    plt.savefig(filesPath + 'Simulation_Plot.png', dpi=400)

    # create a plot for noisy data
    if noiseLevel != 0:
        vars = {'P':[],'M':[]}
        for e in vars.keys():
            plt.figure(2)
            for k in np.arange(0,len(outputSpecies)):
                vars[e].append(e + str(k))
                if exportData[1]=='P':
                    plt.subplot(2,1,1)
                    plt.title('Protein Count Vs. Time (Noisy)', fontsize=GRAPH_TITLE_FONTSIZE)
                    plt.grid(color='k', linestyle='-', linewidth=.4)
                    plt.ylabel('count',fontsize=GRAPH_LABEL_FONTSIZE)
                    plt.xlim(0,inputData[1])
                    plt.yticks(fontsize=GRAPH_LABEL_FONTSIZE)
                    plt.plot (resultNoisy[:,0],resultNoisy[:,outputSpecies[k]], label = vars[e][k],color=colors[k])
                    plt.legend(vars[e])
                    #plt.yscale('log')
                    if inputData[1] > 5000:
                        plt.xscale('log')
                        plt.xticks(fontsize=6)
                    else:
                       plt.xticks(tRange, fontsize=GRAPH_LABEL_FONTSIZE)
                    plt.legend(loc='upper right',ncol=3, bbox_to_anchor=(1.13, 1.035), fontsize=GRAPH_LABEL_FONTSIZE)
                if exportData[1]=='M':
                    #plt.subplot(2,1,2)
                    plt.title('mRNA Count Vs. Time (Noisy)', fontsize=GRAPH_TITLE_FONTSIZE)
                    plt.grid(color='k', linestyle='-', linewidth=.4)
                    plt.xlim(0, inputData[1])
                    plt.yticks(fontsize=GRAPH_LABEL_FONTSIZE)
                    plt.xlabel('time(s)',fontsize=GRAPH_LABEL_FONTSIZE)
                    plt.ylabel('count',fontsize=GRAPH_LABEL_FONTSIZE)
                    plt.plot (resultNoisy[:,0],resultNoisy[:,outputSpecies[k]+1], label = vars[e][k],color=colors[k])
                    plt.legend(vars[e])
                    #plt.yscale('log')
                    if inputData[1] > 5000:
                        plt.xscale('log')
                        plt.xticks(fontsize=GRAPH_LABEL_FONTSIZE)
                    else:
                        plt.xticks(tRange, fontsize=GRAPH_LABEL_FONTSIZE)
                    plt.legend(loc='upper right',ncol=3, bbox_to_anchor=(1.13, 1.035), fontsize=GRAPH_LABEL_FONTSIZE)
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
        manager.window.showMinimized()
        plt.savefig(filesPath + 'Noisy Simulation_Plot.png', dpi=400)
        print('\nPlots created!\n')

    # Attepts to print meaningful Errors
def ErrorPrinting(Error,customHeader=''):
    errorStr = str(Error)
    print('!'*len(errorStr))
    print(customHeader + errorStr)
    print('!'*len(errorStr))
    return()