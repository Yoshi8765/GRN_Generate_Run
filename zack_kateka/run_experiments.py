# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 12:13:40 2018

@author: Kateka Seth
"""
from RunModel import run_model


def export_experiments(csv_file="BIOEN 498_ Experiment Request Form.csv", ant_file="pathway_antimony.txt"):
    ant_str = open(ant_file, 'r').read()
        
    f = open(csv_file)
    i = 0
    print("here")
    for line in f:
        if i != 0:
            print("here")
            line = line.replace("\"", "")
            words = line.split(",")
            print(line)
            print(words)
            team = words[2]
            # process pertubations
            if "Up" in words[3]:
                pert = "UP"
            elif "Down" in words[3]:
                pert = "DOWN"
            elif "Deletion" in words[3]:
                pert = "KO"
            else: 
                pert = "Wild"
            if pert != "Wild":
                pert_gene = list_to_ints(words[4].split(";"))
            else:
                pert_gene = [0]
                
            # process experiment
            if "Mass Spectrometry" in words[5]:
                name = "MassSpec"
                selections = list(range(1,9))
                flag = "M"
            elif "RNA" in words[5]:
                name = "RNASeq"
                selections = list(range(1,9))
                flag = "P"
            else: #words[5] == "Fluorescence Tagging (up to 3 proteins)"
                name = "Fl"
                selections = list_to_ints(words[6].split(";"))
                flag = "P"
            
            # select time course points
            if name == "MassSpec" or name == "RNASeq":
                if "Low" in words[5]:
                    resolution = 20
                else:
                    resolution = 10
            else: # flourescence
                resolution = 10
                    
            
            # create file names
            savePath = team + "/"
            savePath += team + "_" + pert + "_" + convert_list(pert_gene) + name
            if name == "Fl":
                savePath += "_" + convert_list(selections)
            savePath = savePath.replace(" ", "_")
            print(savePath)
            inputData = [1, 200, resolution, pert_gene, [pert, 35, 4]]
            exportData = [selections, flag, True, False, True]
            
#           def run_model(antStr,noiseLevel,inputData=None,exportData=None,bioTap='',
#                         savePath='\\model_output\\',showTimePlots=False,seed=0,drawModel=None,runAttempts=5):
            run_model(ant_str, noiseLevel=5, inputData=inputData, exportData=exportData, savePath=savePath)
        i = 1
        
############## Helper functions ##############
def convert_list(genes):
    result = ""
    for i in genes:
        result += str(i) + "_"
    return result

def list_to_ints(genes):
    for i in range(0, len(genes)):
        genes[i] = int(genes[i])
    return genes
# testing code
export_experiments()