# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 12:13:40 2018

@author: Kateka Seth
"""
import os
from RunModel_2 import run_model2
#from RunModel import run_model
import smtplib 
from email.mime.multipart import MIMEMultipart 
from email.mime.text import MIMEText 
from email.mime.base import MIMEBase 
from email import encoders 


"""
Given the csv from google forms, will parse through and run the correct experiments
for each entry. Will update the team's money and send email with the csv of the experiment
results to the student who filled the form.

This function will generate a file called "run_experiments_data.txt" which keeps 
tracks of the last timestamp in the csv so that it does not run experiments twice.

csv_file: File location of the google forms csv. Default is the name downloaded 
          off of google in current directory.
ant_file: File location of a txt that contains the antimony string of the pathway.
          Default is the file given by GetModel with model name "pathway"
team_file: File location of the csv file containing the team scores. 
sendEmail: Set to true to turn on email sending functionalities.
updateMoney: Set to true to update team money by overwritting team_file.
"""
def export_experiments(num_genes, csv_file="BIOEN 498_ Experiment Request Form.csv", ant_file="pathway_antimony.txt",
                       team_file="team_scores.csv", sendEmail=False, updateMoney=False):
    if os.path.isfile(csv_file) == False:
        raise ValueError(csv_file + " does not exist")
    if os.path.isfile(team_file) == False:
        raise ValueError(team_file + " does not exist")
    if os.path.isfile(ant_file) == False:
        raise ValueError(ant_file + " does not exist")
    if type(num_genes) != int:
        raise ValueError("num_genes parameter should be an integer")
    if type(sendEmail) != bool:
        raise ValueError("sendEmail should be a boolean")
    if type(updateMoney) != bool:
        raise ValueError("updateMoney should be a boolean")
        
    ant_str = open(ant_file, 'r').read()
    f = open(csv_file)
    i = 0
    prevTime = getPrevTime()
    timestamp = 0
    for line in f:
        line = line.replace("\"", "")
        words = line.split(",")
        if i != 0:
            team = words[2]
            email = words[1]
            timestamp = words[0]
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
                selections = list(range(1,num_genes+1))
                species_type = "M"
                money += 1700
            elif "RNA" in words[5]:
                name = "RNASeq"
                selections = list(range(1,num_genes+1))
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
             
            canBuy = update_money(team_file, team, money, updateMoney)    
            money_left = canBuy[1]
            if canBuy[0]: # comment out line to run experiment regardless of money (i.e. for testing)
                savePath = team
                savePath = savePath.replace(" ", "_")
                # make team dir if it doesn't exist
                if os.path.exists(savePath + "/") == False:
                    os.mkdir(savePath)
                saveName = team + "_" + pert + "_" + convert_list(pert_gene) + "_" + name
                if name == "Fl":
                    saveName += "_" + convert_list(selections)
                saveName = saveName.replace(" ", "_")
#                inputData = [1, 200, resolution, pert_gene, [pert, 35, 4]]
#                exportData = [selections, species_type, True, False, True]
                
#                def run_model(antStr,noiseLevel,inputData=None,exportData=None,bioTap='',
#                             savePath='\\model_output\\',showTimePlots=False,seed=0,drawModel=None,runAttempts=5):
#                run_model(ant_str, noiseLevel=0.05, inputData=inputData, exportData=exportData, savePath=savePath)
    
                run_model2(ant_str, noiseLevel=0.05, species_type=species_type, species_nums=selections,
                           timepoints=[200, resolution], exportData=True, perturbs=[(pert, pert_gene)],
                           save_path=savePath, filename=saveName)
                
                path = savePath + "/experimental_data_pathway/" + saveName + ".csv"
                saveName = saveName + ".csv"
                if sendEmail:
                    body = ("Here is your experiment results. You have spent " + str(money) + 
                           ". Your team has " + str(money_left) + " credits left.")
                    send_email(email, body, saveName, path, True)
                    print("Success! Emailed " + email)
            else:
                if sendEmail:
                    body = ("Lacking funds -- experiment has not been run. You have " 
                           + str(money_left) + " credits leftover.")
                    send_email(email, body)
                    print("Emailed " + email)
        if prevTime == "":
            i = 1
        else:
            time = words[0]
            if prevTime == time:
                i = 1
    f.close()
    f = open("run_experiments_data.txt", "w")
    f.write(timestamp)
    f.close()
        
    
"""
Will update the given csv of team credits by subtracting the experiment cost. 
Overwrites the given file with the new updated money.

team_file: File location of the csv file containing the team scores.
team: Name of team purchasing experiment.
money: Cost of the experiment.
updateMoney: Set to true to update team money by overwritting team_file.
"""
def update_money(team_file, team, money, updateMoney):
    f = open(team_file)
    #team = int(team.replace("team ", ""))
    i = 0
    teamNum = 0
    header=""
    words=[]
    canBuy = True
    money_left = 0
    for line in f:
        words = line.split(",")
        if i == 0:
            header = line
            for j in range(0, len(words)):
                if words[j].strip() == team.strip():
                    teamNum = j
        if i == 1:
            team_money = int(words[teamNum])
            money_left = team_money
            if team_money - money < 0:
                canBuy=False
                print("Team " + str(team) + " only has " + str(team_money) + 
                      ". Cannot buy experiment that costs " + str(money) + ".")
            else:
                words[teamNum] = team_money - money
                money_left = words[teamNum]
        i += 1
    f.close()
    if updateMoney:
        # write new file
        f = open(team_file, 'w')
        f.write(header)
        f.write(str(words[0]))
        for i in range(1, len(words)):
            f.write("," + str(words[i]))
        f.close()
    return [canBuy, money_left]


"""
Sends an email to students with data for their requested experiment using 
bioen498@gmail.com as the sender. Will tell students how much money they spent
and how much they have left.

toaddr: student's email
body: a string of the message to be sent in the email
filename: name of attachment file (including extension)
path: path to filename (including filename)
attachment: set to true to send email with attachments.
"""
def send_email(toaddr, body, filename=None, path=None, attachment=False):
    fromaddr = "bioen498@gmail.com"
    # instance of MIMEMultipart 
    msg = MIMEMultipart() 
    
    msg['From'] = fromaddr 
    msg['To'] = toaddr 
    msg['Subject'] = "BIOEN 498 Experiment Data"
    msg.attach(MIMEText(body, 'plain')) 
    
    if attachment:  
        # open the file to be sent  
        filename = filename
        attachment = open(path, "r") 
        # instance of MIMEBase and named as p 
        p = MIMEBase('application', 'octet-stream') 
        # To change the payload into encoded form 
        p.set_payload((attachment).read())   
        # encode into base64 
        encoders.encode_base64(p) 
        p.add_header('Content-Disposition', "attachment; filename= %s" % filename) 
        # attach the instance 'p' to instance 'msg' 
        msg.attach(p) 
    
    # creates SMTP session 
    s = smtplib.SMTP('smtp.gmail.com', 587) 
    # start TLS for security 
    s.starttls() 
    # Authentication 
    s.login(fromaddr, "uwbioenrules!") 
    # Converts the Multipart msg into a string 
    text = msg.as_string() 
    # sending the mail 
    s.sendmail(fromaddr, toaddr, text) 
    # terminating the session 
    s.quit() 

 
"""
Returns the previous time stamp or "" if the timestamp file doesn't exist.
"""
def getPrevTime():
    if os.path.isfile("run_experiments_data.txt"):
        f = open("run_experiments_data.txt", 'r')
        time = f.read()
        f.close()
        return time
    else:
        return "" 
    
    
############## Helper functions ##############
def convert_list(genes):
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