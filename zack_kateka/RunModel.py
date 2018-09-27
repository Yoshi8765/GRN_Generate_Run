# -*- coding: utf-8 -*-
#     Author: Yoshi Goto

import os
import tellurium as te
import numpy as np
import math
from matplotlib import pyplot as plt
import imp
import pandas as pd

# TODO: Make a table of appropriate ranges for parameters.
'''
def run_model2(antStr, noiseLevel, species_type, species_nums, timepoints,
        exportData=False,input_conc=1, perturbs = [], perbParam=[.20,.50], bioTap='',
              save_path = os.getcwd(), filename = "results", showTimePlots=False,seed=0,runAttempts=5):
'''
'''
run_model2(ant_str, noiseLevel=0.05, species_type=species_type, species_nums=selections,
                           timepoints=[200, resolution], exportData=True, perturbs=[(pert, pert_gene)],
                           save_path=savePath, filename=saveName)
'''


def run_model(antStr,noiseLevel,inputData=None,exportData=None,bioTap='',
              savePath='\\results\\',showTimePlots=False,seed=0,drawModel=None):
    """
    Checks if Antimony models will reach steady-state, generates visualizations, and exports data.

    Arguments:
        antStr (str) :
            - The Antimony model.

        noiseLevel (float):
            - level of noise as a decimal. Ex: `0.05` = 5%

        inputData (list: [inputVal,maxTime,resolution,perturbations,per_params]) :
            - **inputVal** (float) : The initial concentration of the input species (name in model: INPUT). Default: 1
            - **maxTime** (int) : Duration of minutes to simulate model. Default: 100
            - **resolution** (int) : The resolution of datapoints. Eg: If resolution = 5, the data will include data every 5 minutes. Default: 1
            - **perturbations** (int list) : What species to perturb, if perturbation is necessary. Default:  0
            - per_params (list: [type, mean, stdev]) :
                - The parameters needed for perturbation of species.
                - **type** (str) : either `UP` or `DOWN` or `KO` used as a flag for upregulation, repression, or knockout, respectively. Default: ['UP']
                - **mean** (float) : The mean value used in `np.random.normal` to generate a random perturbation. The default is recommended for a perturbation of 20-60%. Default: 35
                - **stdev** (float) : The stdev value used in `np.random.normal` to generate a random pertrubation. The default is recommended for a perturbation of 20-60%. Default: 4

        exportData (list: [genesToExport, speciesDataType, sbml, antimony]) :
            - **genesToExport** (int List) : Which genes you want to export. Pass in [0] to export all proteins. Default: [0]
            - **speciesDataType** (str) : `P`rotein or `M`rNA. Default: 'P'
            - **sbml** (bool) : This is a flag for the export of the model in SBML format (version based on Tellurium). Default:True
            - **antimony** (bool) : This is a flag for the export of the Antimony model as a txt. Default:True

        bioTap (str):
            - If this is not empty, a csv file to use for BioTapestry will be exported. Default: ''

        savePath (str):
            - All output files will be saved to `current_working_directory\savePath\modelName\` . Default: '\\results\\'

        showTimePlots (bool):
            - Flag for if you want plots to be generated and saved. Default: False

        seed (int):
            - Seed for reproducibility purposes when generating random data. Default: 0

        drawModel (bool):
            - Generates a graph with PyGraphViz. Requires pygraphviz and graphviz installed Default: False

    Outputs:
        - model (model object): Roadrunner instance of the model.
        - result (np.array): Numpy matrix of the result data
        - resultNoisy (np.array): Numpy matrix of the result data (Noisy)
    """

    # Creating default lists
    if len(inputData) > 3:
        if len(inputData[4])==1:
            inputData[4] = [inputData[4][0],35,4]
    if inputData is None:
        inputData = [1,100,1, 0,['UP']]

    if exportData is None:
        exportData = [0,'P',True,True]

    if drawModel is None:
        drawModel = [False,'fdp']

    # Load the Antimony string as a model
    try:
        model = te.loadAntimonyModel(antStr)
        print('Successfully loaded in Antimony string as a model.')
    except Exception as error:
        ErrorPrinting(error)
        raise RuntimeError("Failed to load antimony string due to parsing error. Check that your antStr is correct.")

    # Check that the model reaches steady-state. If the model has any issues running, it will raise an Exception.
    # If an exception is reached, the user should run get_model again.
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

    # reset model
    model.resetToOrigin()

    # specify input
    model.INPUT = inputData[0];

    # specify perturbations
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

    # RConstruct the selection arguments for simulation
    genesToExport = exportData[0]
    species_type = exportData[1]

    if genesToExport==0:
        allGenes = int((model.getNumFloatingSpecies() - 1)/2)
        if species_type == 'P':
            selections = ["P" + str(i+1) for i in range(allGenes)]
        elif species_type == 'M':
            selections = ["mRNA" + str(i+1) for i in range(allGenes)]
    else:
        if species_type == 'P':
            selections = ["P" + str(i) for i in genesToExport]
        elif species_type == 'M':
            selections = ["mRNA" + str(i) for i in genesToExport]
    if species_type not in ['P','M']:
        raise ValueError("Output data type not recognized (must be P or M)")

    selections = ['time'] + selections

    # Run a simulation for time-course data
    result = model.simulate(0,inputData[1],tStep+1,selections=selections)

    # Retrieve model name
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
        #### do not want to give students the clean data ####
#        data = {next_name:result[:,i] for i,next_name in enumerate(selections)}
#        df = pd.DataFrame.from_dict(data)
#        df.set_index('time', inplace=True)
#        df.plot()
#        plt.title("Clean Data")
#        plt.xlabel("time")
#        plt.savefig(filesPath + 'Simulation_Plot.png', dpi=400)
#        manager = plt.get_current_fig_manager()
#        manager.window.showMaximized()
#        plt.show()
    if showTimePlots==True:
        data = {next_name:resultNoisy[:,i] for i,next_name in enumerate(selections)}
        df = pd.DataFrame.from_dict(data)
        df.set_index('time', inplace=True)
        df.plot()
        plt.title("Noisy Data")
        plt.xlabel("time")
        plt.savefig(filesPath + 'Simulation_Noisy_Plot.png', dpi=400)
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
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

#### END MAIN FUNCTION ####

def Output(exportData,model,seed,selections,result,noiseLevel,resultNoisy,filesPath, antStr,bioTap):
#    export csv of results
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
    if exportData[2]==True:
        sbmlStr = model.getSBML()
        te.saveToFile (filesPath + 'OrigModel.xml', sbmlStr)
    # export Antimony model text
    if exportData[3]==True:
        if np.DataSource().exists(filesPath +  'OrigAntimony.txt'):
            print('Warning: ' + filesPath + 'OrigAntimony.txt already exists! Preventing overwrite.' )
        else:
            fh = open(filesPath + 'OrigAntimony.txt', 'w')
            fh.write(str(antStr))
            fh.close()
    print('\nData Saved!\n')

    # Attepts to print meaningful Errors
def ErrorPrinting(Error,customHeader=''):
    errorStr = str(Error)
    print('!'*len(errorStr))
    print(customHeader + errorStr)
    print('!'*len(errorStr))
    return()
