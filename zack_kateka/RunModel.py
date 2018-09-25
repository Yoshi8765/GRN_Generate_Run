# -*- coding: utf-8 -*-
#     Author: Yoshi Goto

import os
import tellurium as te
import numpy as np
import math
from matplotlib import pyplot as plt
import imp
import pandas as pd

# Global variables used for plot fontsizes
GRAPH_LABEL_FONTSIZE = 8
GRAPH_TITLE_FONTSIZE = 10

# TODO: Plots do not plot anything above around 50 correctly. Solve this or else fixing plots are useless.
# TODO: Make a table of appropriate ranges for parameters.

def run_model(antStr,noiseLevel,inputData=None,exportData=None,bioTap='',
              savePath='\\model_output\\',showTimePlots=False,seed=0,drawModel=None):
    """Checks if Antimony models will reach steady-state, generates visualizations, and exports data.

    Arguments:
        antStr: (str) The Antimony model.
        noiseLevel: (float) level of noise as a decimal. Ex: 0.05 = 5%

        inputData:
            (list) [INPUTval,maxTime,resolution,perturbations,perturbation_params[up/dowm,mean,stdev]]
            - INPUTval: The initial concentration of the input species (name in model: INPUT)
            - maxTime: Duration of minutes to simulate model
            - resolution: The resolution of datapoints. Eg: If resolution = 5, the data will include data every 5 minutes
            - perturbations: What species to perturb, if perturbation is necessary.
            - perturbation_params: (Input: list of numbers) The parameters needed for perturbation of species.
                - up/down: either `UP` or `DOWN` or `KO` used as a flag for upregulation, repression, or knockout, respectively.
                - mean: The mean value used in np.random.normal to generate a random perturbation. The default is recommended for a perturbation of 20-60%.
                - stdev: The stdev value used in np.random.normal to generate a random pertrubation. The default is recommended for a perturbation of 20-60%.

        exportData:
            (list) [genesToExport,speciesDataType,csv,biotapestryCSV,sbml,antimony,seed]
            genesToExport = (int List) Default: [0] . A list object of which genes you want to export, 1-based indexing. Pass in [0] to export all proteins.
            speciesDataType = (str) Default: 'P' . 'P'rotein or 'M'rNA
            csv = (bool) Default:True . This is a flag for the export of the timecourse data of the genesToExport as csv files.
            sbml = (bool) Default:True . This is a flag for the export of the model in SBML format (version based on Tellurium).
            antimony = (bool) Default:True . This is a flag for the export of the Antimony model as a txt.

        bioTap: (str)  Default: '' . If this is not empty, a csv file to use for BioTapestry will be exported.

        savePath: (str) Default: '\\model_output\\' . All output files will be saved to `current_working_directory\savePath\modelName\`

        showTimePlots: (bool) Default: False . Flag for if you want plots to be generated and saved.

        seed: (int) Default: 0 . Seed for reproducibility purposes when generating random data.

        drawModel: (bool) Default: False . Generates a graph with PyGraphViz. Requires pygraphviz and graphviz installed

    Outputs:
        model: (model object) Roadrunner instance of the model.
        result: (np.array) Numpy matrix of the result data
        resultNoisy: (np.array) Numpy matrix of the result data (Noisy)
    """

    #Creating default lists
    if inputData is None:
        inputData = [1,40,1, [0],['UP',35,4]]
    if len(inputData) == 1:
        inputData = [inputData[0],35,4]

    if exportData is None:
        exportData = [ [0],'P',True,True,True]

    if drawModel is None:
        drawModel = [False,'fdp']

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
            print ("Failed to reach steady-state or there is an Integrator error. Check your model or run get_model again.")
            ErrorPrinting(error)
            raise error

    model.resetToOrigin()

    #Specify input
    model.INPUT = inputData[0];

    #specify perturbations
    if len(inputData)> 3:
        if inputData[3] != [0]:
            perturb = inputData[3]
            pertParam = inputData[4]
            for gene in perturb:
                currVm = eval('model.Vm' + str(gene))
                if pertParam[0] == 'UP':
                     newVal = currVm + np.random.normal(pertParam[1],pertParam[2])/100
                     exec('model.Vm' + str(gene) +  ' = ' + str(newVal)) in globals()
                if pertParam[0] == 'DOWN':
                    newVal = currVm - np.random.normal(pertParam[1],pertParam[2])/100
                    exec('model.Vm' + str(gene)  +  ' = ' + str(newVal)) in globals()
                if pertParam[0] == 'KO':
                    exec('model.Vm' + str(gene)  + ' = 0') in globals()
                    exec('model.d_mRNA' + str(gene)  + ' = 0') in globals()
                    exec('model.d_protein' + str(gene)  + ' = 0') in globals()
                    exec('model.mRNA' + str(gene)  + ' = 1E-9') in globals()
                    exec('model.P' + str(gene)  + ' = 1E-9') in globals()
                    #change initVals to 1E-9 instead of 0 to prevent possible solver hanging bug

    # Run a simulation for time-course data
    if exportData[2]==True:
        genesToExport = exportData[0]
        species_type = exportData[1]

        if genesToExport!=0:
            if species_type == 'P':
                selections = ["P" + str(i) for i in genesToExport]
            elif species_type == 'M':
                selections = ["mRNA" + str(i) for i in genesToExport]
        else:
            allGenes = (model.getNumFloatingSpecies() - 1)/2
            if species_type == 'P':
                selections = ["P" + str(i+1) for i in range(allGenes)]
            elif species_type == 'M':
                selections = ["mRNA" + str(i+1) for i in range(allGenes)]

        if species_type not in ['P','M']:
            raise ValueError("Output data type not recognized (must be P or M)")

    selections = ['time'] + selections

    result = model.simulate(0,inputData[1],tStep+1,selections=selections)


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

    # Export datasets
    Output(exportData,model,seed,selections,result,noiseLevel,resultNoisy, filesPath, antStr,bioTap)

    # Create graphs
    if showTimePlots==True:
        data = {next_name:resultNoisy[:,i] for i,next_name in enumerate(selections)}
        df = pd.DataFrame.from_dict(data)
        df.set_index('time', inplace=True)
        df.plot()
        plt.title("Noisy Data")
        plt.xlabel("time")
        plt.savefig(filesPath + 'Simulation_Noisy_Plot2.png', dpi=400)
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
        # do not want to give students the clean data
        #data = {next_name:result[:,i] for i,next_name in enumerate(selections)}
        #df = pd.DataFrame.from_dict(data)
        #df.set_index('time', inplace=True)
        #df.plot()
        #plt.title("Normal Data")
        #plt.xlabel("time")
        #plt.savefig(filesPath + 'Simulation_Clean_Plot2.png', dpi=400)

        plt.show()

    # Draw the model (requires pygraphviz module)
    if drawModel[0]==True:
        try:
            imp.find_module('pygraphviz')
            model.draw(layout=str(drawModel[1]))
        except ImportError:
            print('pygraphviz is not installed!')

    #returns the model, and result and/or resultNoisy arrays
    return(model,result,resultNoisy)

###### Other Functions that GetModel uses ######

def Output(exportData,model,seed,selections,result,noiseLevel,resultNoisy,filesPath, antStr,bioTap):
#     export csv of results
    data = {next_name:result[:,i] for i, next_name in enumerate(selections)}
    df = pd.DataFrame.from_dict(data)
    df.set_index('time', inplace=True)
    df.to_csv(filesPath + "Result_Clean.csv")

    if noiseLevel != 0:
        data = {next_name:resultNoisy[:,i] for i, next_name in enumerate(selections)}
        df = pd.DataFrame.from_dict(data)
        df.set_index('time', inplace=True)
        df.to_csv(filesPath + "Noisy_Result.csv")

#     export csv file for importing into Biotapestry
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
        if np.DataSource().exists(filesPath +  'OrigAntimony.txt'):
            print('Warning: ' + filesPath + 'OrigAntimony.txt already exists! Preventing overwrite.' )
        else:
            fh = open(filesPath + 'OrigAntimony.txt', 'w')
            fh.write(str(antStr))
            fh.close()
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

    fh = open(filesPath, "w")
    fh.write(header_string)
    fh.write('\n')
    for x in range(data.shape[0]):
        for y in range(data.shape[1]):
            val = data[x][y]
            s1 = "{:6.5f}".format(val)
            fh.write (s1)
            if y != data.shape[1] - 1:
                fh.write (',')
        fh.write('\n')
    fh.close()
    return(header_string)

#def makePlots(exportData,inputData,filesPath,header_string,noiseLevel,resultNoisy):
# Create graphs

    #    global GRAPH_LABEL_FONTSIZE
#    global GRAPH_TITLE_FONTSIZE
#
#    outputSpecies = range(int(np.size(result,1)))
#    if exportData[0]!=0:
#        outputSpecies = [gene for gene in exportData[0]]
#    else:
#        outputSpecies = outputSpecies[2::2]
#
#    # create a plot for clean data
#    jet = plt.get_cmap('jet')
#    pickColorSpace = [int(np.floor(i)) for i in np.linspace(0,255,len(outputSpecies))]
#    colors = jet(pickColorSpace)
#    colors = jet(np.linspace(0,1,len(outputSpecies)))
#
#    vars = {'P':[],'M':[]}
#    for e in vars.keys():
#        plt.figure(1)
#        for k in np.arange(0,len(outputSpecies)):
#            vars[e].append(e + str(k))
#            tStep = int(math.ceil(inputData[1]/inputData[2]))
#            tRange = np.arange(0,(int(math.ceil(inputData[1])))+ tStep,tStep)
#            if exportData[1]=='P':
#               # plt.subplot(2,1,1)
#                plt.title('Protein Count Vs. Time', fontsize=GRAPH_TITLE_FONTSIZE)
#                plt.grid(color='k', linestyle='-', linewidth=.4)
#                plt.ylabel('count',fontsize=GRAPH_LABEL_FONTSIZE)
#                plt.yticks(fontsize=GRAPH_LABEL_FONTSIZE)
#                plt.plot (result[:,0],result[:,outputSpecies[k]], label = vars[e][k],color=colors[k])
#                #plt.yscale('log')
#                if inputData[1] > 5000:
#                    plt.xscale('log')
#                    plt.xticks(fontsize=GRAPH_LABEL_FONTSIZE)
#                else:
#                    plt.xticks(tRange, fontsize=GRAPH_LABEL_FONTSIZE)
#                plt.xlim(0,inputData[1])
#                plt.legend(loc='upper right',ncol=3, bbox_to_anchor=(1.13, 1.035), fontsize=GRAPH_LABEL_FONTSIZE)
#            if exportData[1]=='M':
#              # plt.subplot(2,1,2)
#               plt.title('mRNA Count Vs. Time', fontsize=GRAPH_TITLE_FONTSIZE)
#               plt.grid(color='k', linestyle='-', linewidth=.4)
#               plt.xlabel('time(s)',fontsize=GRAPH_LABEL_FONTSIZE)
#               plt.ylabel('count',fontsize=GRAPH_LABEL_FONTSIZE)
#               plt.yticks(fontsize=GRAPH_LABEL_FONTSIZE)
#               plt.plot (result[:,0],result[:,outputSpecies[k]+1])#, label = vars[e][k],color=colors[k])
#              # plt.yscale('log')
#               if inputData[1] > 5000:
#                    plt.xscale('log')
#                    plt.xticks(fontsize=GRAPH_LABEL_FONTSIZE)
#               else:
#                   plt.xticks(tRange, fontsize=GRAPH_LABEL_FONTSIZE)
#               plt.xlim(0,inputData[1])
#               plt.legend(vars[e],loc='upper right',ncol=3, bbox_to_anchor=(1.13, 1.035), fontsize=GRAPH_LABEL_FONTSIZE)
#    manager = plt.get_current_fig_manager()
##    manager.window.showMaximized()
##    manager.window.showMinimized()
#    plt.savefig(filesPath + 'Simulation_Plot.png', dpi=400)
#
#    # create a plot for noisy data
#    if noiseLevel != 0:
#        vars = {'P':[],'M':[]}
#        for e in vars.keys():
#            plt.figure(2)
#            for k in np.arange(0,len(outputSpecies)):
#                vars[e].append(e + str(k))
#                if exportData[1]=='P':
#                    plt.subplot(2,1,1)
#                    plt.title('Protein Count Vs. Time (Noisy)', fontsize=GRAPH_TITLE_FONTSIZE)
#                    plt.grid(color='k', linestyle='-', linewidth=.4)
#                    plt.ylabel('count',fontsize=GRAPH_LABEL_FONTSIZE)
#                    plt.xlim(0,inputData[1])
#                    plt.yticks(fontsize=GRAPH_LABEL_FONTSIZE)
#                    plt.plot (resultNoisy[:,0],resultNoisy[:,outputSpecies[k]], label = vars[e][k],color=colors[k])
#                    plt.legend(vars[e])
#                    #plt.yscale('log')
#                    if inputData[1] > 5000:
#                        plt.xscale('log')
#                        plt.xticks(fontsize=6)
#                    else:
#                       plt.xticks(tRange, fontsize=GRAPH_LABEL_FONTSIZE)
#                    plt.legend(loc='upper right',ncol=3, bbox_to_anchor=(1.13, 1.035), fontsize=GRAPH_LABEL_FONTSIZE)
#                if exportData[1]=='M':
#                    #plt.subplot(2,1,2)
#                    plt.title('mRNA Count Vs. Time (Noisy)', fontsize=GRAPH_TITLE_FONTSIZE)
#                    plt.grid(color='k', linestyle='-', linewidth=.4)
#                    plt.xlim(0, inputData[1])
#                    plt.yticks(fontsize=GRAPH_LABEL_FONTSIZE)
#                    plt.xlabel('time(s)',fontsize=GRAPH_LABEL_FONTSIZE)
#                    plt.ylabel('count',fontsize=GRAPH_LABEL_FONTSIZE)
#                    plt.plot (resultNoisy[:,0],resultNoisy[:,outputSpecies[k]+1], label = vars[e][k],color=colors[k])
#                    plt.legend(vars[e])
#                    #plt.yscale('log')
#                    if inputData[1] > 5000:
#                        plt.xscale('log')
#                        plt.xticks(fontsize=GRAPH_LABEL_FONTSIZE)
#                    else:
#                        plt.xticks(tRange, fontsize=GRAPH_LABEL_FONTSIZE)
#                    plt.legend(loc='upper right',ncol=3, bbox_to_anchor=(1.13, 1.035), fontsize=GRAPH_LABEL_FONTSIZE)
#        manager = plt.get_current_fig_manager()
##        manager.window.showMaximized()
##        manager.window.showMinimized()
#        plt.savefig(filesPath + 'Noisy Simulation_Plot.png', dpi=400)
#        print('\nPlots created!\n')

    # Attepts to print meaningful Errors
def ErrorPrinting(Error,customHeader=''):
    errorStr = str(Error)
    print('!'*len(errorStr))
    print(customHeader + errorStr)
    print('!'*len(errorStr))
    return()