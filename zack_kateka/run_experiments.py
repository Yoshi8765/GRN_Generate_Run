# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 12:13:40 2018

@author: Kateka Seth
"""
import os
from RunModel_2 import run_model2

# TODO: add email functionality??
def export_experiments(csv_file="BIOEN 498_ Experiment Request Form.csv", ant_file="pathway_antimony.txt",
                       team_file="team_scores.csv"):
    ant_str = open(ant_file, 'r').read()
    f = open(csv_file)
    i = 0
    for line in f:
        if i != 0:
            line = line.replace("\"", "")
            words = line.split(",")
            print(line)
            team = words[2]
            # process pertubations
            if "Up" in words[3]:
                pert = "UP"
                money = 350
            elif "Down" in words[3]:
                pert = "DOWN"
                money = 350
            elif "Deletion" in words[3]:
                pert = "KO"
                money = 800
            else: 
                pert = "Wild"
                money = 0
            if pert != "Wild":
                pert_gene = list_to_ints(words[4].split(";"))
                money *= len(pert_gene)
            else:
                pert_gene = [0]
                
            # process experiment
            if "Mass Spectrometry" in words[5]:
                name = "MassSpec"
                selections = list(range(1,9))
                species_type = "M"
                money += 1700
            elif "RNA" in words[5]:
                name = "RNASeq"
                selections = list(range(1,9))
                species_type = "P"
                money += 1500
            else: #words[5] == "Fluorescence Tagging (up to 3 proteins)"
                name = "Fl"
                selections = list_to_ints(words[6].split(";"))
                species_type = "P"
                money += 300 * len(selections)
                if len(selections) == 3:
                    money += 50         
            
            # select time course points
            if name == "MassSpec" or name == "RNASeq":
                if "Low" in words[5]:
                    resolution = 20
                else:
                    resolution = 10
                    money += 1500 
            else: # flourescence
                resolution = 10
                
            canBuy = update_money(team_file, team, money)    
            print(selections)
            if canBuy:
                # make team dir
                if os.path.exists(team + "/") == False:
                    os.mkdir(team)
                # create file names
                savePath = team + "/"
                savePath += team + "_" + pert + "_" + convert_list(pert_gene) + "_" + name
                if name == "Fl":
                    savePath += "_" + convert_list(selections)
                savePath = savePath.replace(" ", "_")
                print(savePath)
#                inputData = [1, 200, resolution, pert_gene, [pert, 35, 4]]
#                exportData = [selections, species_type, True, False, True]
                
    #           def run_model(antStr,noiseLevel,inputData=None,exportData=None,bioTap='',
    #                         savePath='\\model_output\\',showTimePlots=False,seed=0,drawModel=None,runAttempts=5):
#                run_model(ant_str, noiseLevel=0.05, inputData=inputData, exportData=exportData, savePath=savePath)
                run_model2(ant_str, noiseLevel=0.05, species_type=species_type, species_nums=selections,
                           timepoints=[200, resolution], exportData=True, perturbs=[(pert, pert_gene)],
                           save_path=savePath, num_genes=8)
        i = 1
 
    
def update_money(team_file, team, money):
    f = open(team_file)
    team = int(team.replace("team ", ""))
    i = 0
    header=""
    words=[]
    canBuy = True;
    for line in f:
        if i == 0:
            header = line
        if i == 1:
            words = line.split(",")
            team_money = int(words[team - 1])
            if team_money - money < 0:
                canBuy=False
                print(team + "only has " + str(team_money) + ". Cannot buy experiment.")
            else:
                words[team-1] = team_money - money
        i += 1
    f.close()
    
    # write new file
    f = open(team_file, 'w')
    f.write(header)
    f.write(str(words[0]))
    for i in range(1, len(words)):
        f.write("," + str(words[i]))
    f.close()
    return canBuy
     
############## Helper functions ##############
def convert_list(genes):
    print(genes)
    result = ""
    result += str(genes[0])
    for i in range(1, len(genes)):
        result += "_" + str(genes[i])
    return result

def list_to_ints(genes):
    for i in range(0, len(genes)):
        genes[i] = int(genes[i])
    return genes
##############################################


# testing code
export_experiments()