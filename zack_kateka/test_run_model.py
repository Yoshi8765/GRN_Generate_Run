# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 10:46:56 2018

@author: Yoshi
"""

from RunModel import run_model
from GetModel import get_model
from GetModel import convert_to_biotapestry
#from RunModel_2 import run_model2

model_name='exp5'
numGenes = 8
seednum = 123468

antStr,biotap_str = get_model(numGenes,model_name=model_name,seed=seednum,export=False)

#antStr= open('C:\Users\Yoshi\Documents\GitHub\DREAM-work\zack_kateka\Random_GRNs\exp3\OrigAntimony_up7.txt','r').read()
noiseLevel = 0.00 # put in a percentage. 0.05 = 5%
tmax=200 # minutes. The complete data will have tmax*5 datapoints
resolution = 1 # The number of minutes between each timepoint you want in your output.
<<<<<<< HEAD
r,res,resN = run_model(antStr,noiseLevel,exportData=[0,'P',True,True,True,True],inputData=[1,tmax,resolution],showTimePlots=True,bioTap=biotap_str)
=======
#perturb = np.random.normal(35,4)/100
r,res,resN = run_model(antStr,noiseLevel,exportData=[0,'M',True,True,True,True],inputData=[1,tmax,resolution],seed=seednum, showTimePlots=True,bioTap=biotap_str)
>>>>>>> 391963d9507708405e72632ec2867435a6ed0569

#r,res,resN = run_model(antStr,noiseLevel,exportData=[[2,8],'P',True,False,False],inputData=[1,tmax,resolution],showTimePlots=True)

print('done!')

#%%

from change_biotapestry import remove_biotapestry


"""
Given a biotapestry csv file, removes connections specified by remove. Prints
a new csv file to given location.

:param remove: A list containing int tuples of connections to be remove. The tuple
            should be formatted as (source,target). Ex: (1,3) to remove connection
           starting at Gene 1 going to Gene 3. Use "INPUT" for "INPUT" box.
:param csv_filename: Location with name of the original csv file
:param csv_newfile: Location with new file name
"""

#rem = [(8,1),(1,5),(4,7),(5,3),(6,6)]
rem = [(5,3),(6,6),(8,1),(1,5),(4,7)]
#rem = [(8,8),(1,5),(5,3)]

#remove_biotapestry(rem,'C:\Users\Yoshi\Documents\GitHub\DREAM-work\zack_kateka\Random_GRNs\exp3\\biotapestry.csv','stestWorking.csv')
print('done')

#I don't think we were going to have as an assumption that 3 connections are gone.
