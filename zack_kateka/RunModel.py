# -*- coding: utf-8 -*-
#     Author: Yoshi Goto

import os
import tellurium as te
import numpy as np
import math
from matplotlib import pyplot as plt
import pkg_resources

# TODO: check the filesave paths

def run_model(antStr,model_name,noiseLevel,numGenes,runAttempts=5,tmax=60,savePath='\\model_output\\',showTimePlots='y',exportData='y',drawModel=['y','fdp']):
    """Asserts Antimony models will reach steady-state, generates visualizations, and exports data.

    get_model is required to be in the same working directory as run_model.\n
    get_model will be executed if run_model finds out that the current model does not reach steady-state.

    Arguments:
        antStr:
        model_name:
        noiseLevel:
        numGenes:
        runAttempts:
        tmax:
        savePath:
        showTimePlots:
        exportData:
        drawModel: Requires pygraphviz

    Outputs:
        model:
        result:
        resultNoisy:
    """
    # Attempt to run RunModel up to runAttempts times with different generated models.

    # Check that get_model can be used
    try:
        from GetModel import get_model
    except:
        print('Cannot find get_model to import.')

    for retry in range(runAttempts):
        # Load the Antimony string as a model
        try:
            model = te.loadAntimonyModel(antStr)
            print('Success!') # TODO: Change to a better print
        except Exception as error:
            ErrorPrinting(error)
            raise RuntimeError("Failed to load antimony text due to parsing error. Check that your antStr is correct.")
        # Check that the model reaches steady-state
        try:
            plt.close('all')
            model.reset()
            tStep = int(math.ceil(tmax)*5)
            result = model.simulate(0,tmax, tStep)
        # If the model has any issues running, generate a new model and try again.
        except Exception as error:
            print "Failed to solve model."
            ErrorPrinting(error)
            model,antStr,GRNInt = get_model(num_genes, reg_probs, model_name, init_params,seed)
            # TODO: Fix the above
            break
        else:
            # Specify (and make if necessary) a folder to save outputs to
            folderPath = os.getcwd() + '\\Random_GRNs\\'
            modelPath = folderPath + model_name + '\\'
            if os.path.exists(modelPath) == False:
                if os.path.exists(folderPath) == False:
                    os.mkdir(folderPath)
                os.mkdir(modelPath)
                print('\nFolder created: ' + modelPath)

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
            if showTimePlots=='y':
                makePlots(tmax,numGenes,result,resultNoisy)
                # TODO: make plots generalized to how many outputs that the user wants

            # Export data
            if exportData=='y':
                Output(data, Name,model,seed,result,resultNoisy, fileName, antStr)
                # TODO: Fix the above

            # Draw the model (requires pygraphviz)
            if drawModel[0]=='y':
                found = 0
                installed_packages_list = sorted(["%s==%s" % (i.key, i.version) for i in pkg_resources.working_set])
                for item in installed_packages_list:
                    if item.startswith('pygraphviz'):
                        model.draw(layout=str(drawModel[1]))
                        found = 1
                if found == 0:
                    print('pygraphviz is not installed!')

            return(model,result,resultNoisy)

    raise RuntimeError('Could not create a working model for ' + str(runAttempts) + ' tries.')










###### Other Functions that GetModel uses ######

def Output(data, Name,model,seed,result,resultNoisy, fileName, antStr):
    for i in range(len(data)):
        if i == 0 and data[i] == 'y':
            writecsvFile(fileName + Name + '_Results.csv',model,result)
        elif data[i] == 'y':
            writecsvFile(fileName + Name + '_Noisy_Result.csv',model,resultNoisy)
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
    print('\nData Saved!\n')

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

def makePlots(tmax,numGenes,result,resultNoisy):
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
                plt.plot (resultNoisy[:,0],resultNoisy[:,(k-1)+1], label = vars[e][k-1],color=colors[k-1])
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
                plt.plot (resultNoisy[:,0],resultNoisy[:,(k-1)+numGenes+1], label = vars[e][k-1],color=colors[k-1])
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

    #Attepts to show meaningful Errors
def ErrorPrinting(Error,customHeader=''):
    errorStr = str(Error)
    print('!'*len(errorStr))
    print(customHeader + errorStr)
    print('!'*len(errorStr))
    return()