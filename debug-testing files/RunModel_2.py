# -*- coding: utf-8 -*-
#     Author: Yoshi Goto

import os
import tellurium as te
import numpy as np
import math
from matplotlib import pyplot as plt
import imp
import pandas as pd

def run_model2(antStr, noiseLevel, species_type, species_nums, timepoints,
        exportData=False,input_conc=1, perturbs = [], perbParam=[.20,.50], bioTap='',
              save_path = os.getcwd(), filename = "results", showTimePlots=False,seed=0,runAttempts=5):
    """Checks if Antimony models will reach steady-state, generates visualizations, and exports data.

    Arguments
        antStr = antimony string for the given model
        noiseLevel = (float) level of noise as a decimal. Ex: 0.05 = 5%
        species_type = 'P' for protein, or 'M' for mRNA
        species_nums = The gene numbers for associated to the species you are interested in.
            i.e. species_type = 'P' and species_nums =[1,2,3,4] will return the data for P1, P2, P3, and P4
        timepoints = specifies length and resolution of simulation/experimental data -> [max time value, resolution (i.e. stepsize)]:
        exportData = whether or not to include additional files in the output. These files include the antimony string and the SBML Model
            NOTE: should probably stay at FALSE for student usage
        input_conc = what the value INPUT should be set to before running model
        perturbs = Perturbations to apply -> [(perturbation type, [target1, target2, ...]), (<next perturbation)>)]
            i.e. perturbs =[ ("UP", [1, 2, ...]),  ("DOWN", [4,5,6])]
        perbParam = For Up/down regulation, this specifies the min/max amount you want to change regulation -> perbParm = (min percent change, max percent change)
            i.e. perbParam = (0.20, 0.50) specifies that you want the up/down regulation to change transcription between 20 and 50 percent.
            NOTE: percents are chosen uniformally between the min and max
        bioTap: (str) If this is not empty, a csv file to use for BioTapestry will be exported.
        save_path = directory to save output folder to
        filename = file to save the experimental data to
        showTimePlots: whether or not to output timeplots of the data (as png) in the output file
        seed: random number generator seed
        runAttempts: number of attempts to load/run model before exiting and suggesting to rebuild the model


    Outputs:
        model: roadrunner model for the full working GRN
        result: experimental data without noise
        resultNoisy: experimental data with noise
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
        tStep = int(math.ceil(timepoints[0]/timepoints[1]))
        try:
            model.steadyState()
        except Exception:
            try:
                model.conservedMoietyAnalysis = True
                model.steadyState()
            except Exception as error:
                print( "Failed to reach steady-state or there is an Integrator error. Check your model or run get_model again.")
                ErrorPrinting(error)
                raise error

        model.resetToOrigin()

        #Specify input
        model.INPUT = input_conc

        # set seed
        if seed != 0:
            np.random.seed(seed)

        if len(perturbs) > 0:
        #specify perturbations
            for next_type, targets in perturbs:
                for gene in targets:
                    currVm = eval("model.Vm" + str(gene))
                    if next_type  == 'UP':
                        newVm  = currVm + np.random.uniform(perbParam[0]*currVm, perbParam[1]*currVm)
                        exec("model.Vm" + str(gene) + " = " + str(newVm))
                    elif next_type  == 'DOWN':
                        newVm  = currVm - np.random.uniform(perbParam[0]*currVm, perbParam[1]*currVm)
                        exec("model.Vm" + str(gene) + " = " + str(newVm))
                    elif next_type  == 'KO':
                        exec("model.L" + str(gene) + " = 0")
                        exec('model.Vm' + str(gene)  + ' = 0')
                        exec('model.d_mRNA' + str(gene)  + ' = 0')
                        exec('model.d_protein' + str(gene)  + ' = 0')
                        exec('model.mRNA' + str(gene)  + ' = 1E-9')
                        exec('model.P' + str(gene)  + ' = 1E-9')
                        #change initVals to 1E-9 instead of 0 to prevent possible solver hanging bug
                    else:
                        raise ValueError("Perturbation type is not recognized; expected UP, DOWN, or KO")

        # Run a simulation for time-course data
        if species_type == 'P':
            selections = ["P" + str(i) for i in species_nums]
        elif species_type == 'M':
            selections = ["mRNA" + str(i) for i in species_nums]
        else:
            raise ValueError("Output data type not recognized (must be P or M)")

        selections = ['time'] + selections
        result = model.simulate(0,timepoints[0],tStep+1, selections=selections)

        model_name = model.getInfo().split("'modelName' : ")[1].split("\n")[0]

        # Specify (and make if necessary) a folder to save outputs to
        filesPath = save_path + "/experimental_data_" + model_name + '/'
        if os.path.exists(filesPath) == False:
            os.mkdir(filesPath)
            print('\nFolder created: ' + filesPath)

        # Create results with artificial noise
        noise = np.random.normal(1, noiseLevel, size=result.shape)
        noise[:,0] = 1 # we dont want there to be noise in the time
        resultNoisy = np.multiply(result, noise)

        # Create graphs
        if showTimePlots==True:
            data = {next_name:resultNoisy[:,i] for i,next_name in enumerate(selections)}
            df = pd.DataFrame.from_dict(data)
            df.set_index('time', inplace=True)
            df.plot()
            plt.title("Noisy Data")
            plt.xlabel("time")
            plt.savefig(filesPath + 'Simulation_Noisy_Plot2.png', dpi=400)

            # do not want to give students the clean data
            #data = {next_name:result[:,i] for i,next_name in enumerate(selections)}
            #df = pd.DataFrame.from_dict(data)
            #df.set_index('time', inplace=True)
            #df.plot()
            #plt.title("Normal Data")
            #plt.xlabel("time")
            #plt.savefig(filesPath + 'Simulation_Clean_Plot2.png', dpi=400)

            plt.show()

        # Export datasets
        Output(exportData,model,seed,result,noiseLevel,resultNoisy, filesPath, antStr,bioTap, selections, filename)

        #returns the model, and result and/or resultNoisy arrays
        return(model,result,resultNoisy)

    raise RuntimeError('Could not create a working model for ' + str(runAttempts) + ' tries.')


###### Other Functions that GetModel uses ######

def Output(exportData,model,seed,result,noiseLevel,resultNoisy,filesPath, antStr,bioTap, selections, filename):
    # export csv of results

    # do not want to give students the clean data
    #data = {next_name:result[:,i] for i,next_name in enumerate(selections)}
    #df = pd.DataFrame.from_dict(data)
    #df.set_index('time', inplace=True)
    #df.to_csv(filesPath + "Results_Clean2.csv")

    data = {next_name:resultNoisy[:,i] for i, next_name in enumerate(selections)}
    df = pd.DataFrame.from_dict(data)
    df.set_index('time', inplace=True)
    df.to_csv(filesPath + filename + ".csv")

    # export csv file for importing into Biotapestry
    if bioTap!='':
        f2 = open(filesPath + "biotapestry.csv", 'w')
        f2.write(bioTap)
        f2.close()
    #  export SBML model text
    if exportData==True:
        sbmlStr = model.getSBML()
        te.saveToFile (filesPath + 'OrigModel.xml', sbmlStr)
    # export Antimony model text
    if exportData==True:
        if np.DataSource().exists(filesPath +  'OrigAntimony.txt'):
            print('Warning: ' + filesPath + 'OrigAntimony.txt already exists! Preventing overwrite.' )
        else:
            fh = open(filesPath + 'OrigAntimony.txt', 'w')
            fh.write(str(antStr))

    print('\nData Saved!\n')


    # Attepts to print meaningful Errors
def ErrorPrinting(Error,customHeader=''):
    errorStr = str(Error)
    print('!'*len(errorStr))
    print(customHeader + errorStr)
    print('!'*len(errorStr))
    return()
