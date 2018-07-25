# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 14:52:12 2018

@author: Alan
"""
import tellurium as te
import numpy as np
import math
import matplotlib.pyplot as plt
import os

#def prints(what): Used when making prog
#    print(what.keys())
#    for j in what.keys():
#        print(what[j])
        
def RandomGRN(genes,Probabilities, InitVals, Perturbation):
    
    
    Interactions = ["SingleAct", "SingleRep", "and", "or", "Nand", "Nor", "Counter"]
    GRN = {'Rate': [], 'mRNA':[], 'Prot': [],'PDeg':[],'mDeg': [],'Leak': [],'TRa':[], 'TRb':[], 'DualDissoc':[], 'Dissoc':[], 'HillCoeff': []}
    GRNInt = {'TF':[],'TFs':[]}
    StringList=[]
    def ModelInitNames():
        for k in range(genes):
            num = str(k+1)
            GRN['Leak'].append("L" + num)
            GRN['Rate'].append('R' + num +' := ')
            GRN['mDeg'].append('Dm' + num)
            GRN['PDeg'].append('Dp' + num)
            GRN['mRNA'].append('M' + num)
            GRN['Prot'].append('P' + num)
            GRN['TRa'].append('TRa' + num)
            GRN['TRb'].append('TRb' + num)
            GRN['HillCoeff'].append('H' + num)
            GRN['Dissoc'].append('K' + num)
            GRN['DualDissoc'].append('Kd' + num)
        return
    def AssignProteins():
        nums =[int(math.ceil(float(genes)/2)), int(math.floor(float(genes)/2))]
        RanAct = np.random.choice(GRN['Prot'], nums[np.random.randint(2)])
        GRNInt.update({'activators': RanAct})
        index1 = np.arange(len(GRNInt['activators']))
        index2 = np.sort(-index1)
        for m in range(len(index1)):
            ActivatorCur =GRNInt['activators'][m]
            Activators1 = GRNInt['activators'][index2[m]:]
            Activators2 =  GRNInt['activators'][:index1[m]]
            if m > len(index1)-2:
                Activators1 = []
            if ActivatorCur in Activators1 or ActivatorCur in Activators2:
                AssignProteins()
        return
    def DetermineInteractors():
        Proteins = []
        ProbabilityMatrix = []
        ProbSingle = Probabilities[0]
        ProbDouble = Probabilities[1]
        ProbCounter = Probabilities[2]
        for h in range(len(GRN['Prot'])):
            Proteins.append(GRN['Prot'][h])
        for s in GRNInt['activators']:#Issue generates too many interactions when activators has multiple regulatory interactions per gene
            if s in Proteins:
                Proteins.remove(s)
            GRNInt.update({'repressors': Proteins})
        Psing = ProbSingle/2
        Pdub = ProbDouble/4
        ProbabilityMatrix = [Psing, Psing, Pdub, Pdub, Pdub, Pdub, ProbCounter]
        GRNInt.update({'interactions': np.random.choice(Interactions,genes,True,ProbabilityMatrix)})
        for Inter in range(len(GRNInt['interactions'])):  
            reg = GRNInt['interactions'][Inter]    
            RegType = 'repressors'
            if 'Single' in reg :
                if 'Act' in reg:
                    RegType = 'activators'
                GRNInt['TF'].append([reg, [GRNInt[RegType][np.random.randint(len(GRNInt[RegType]))]]])
            else:
                Tps={'activators':[], 'repressors':[]}
                for m in Tps.keys():
                    n=1
                    Tps[m].append(np.random.choice(GRNInt['activators'],2))
                    while Tps[m][0][n-1] == Tps[m][0][n]:
                        Tps[m][0][n] =np.random.choice(GRNInt['activators'])
                TypeA = 'activators'
                TypeR = 'repressors'
                if reg == "and" or reg == "or":
                    GRNInt['TF'].append([reg, [Tps[TypeA][0][0],Tps[TypeA][0][1]]])
                elif reg == "Nand" or reg == "Nor":
                    GRNInt['TF'].append([reg, [Tps[TypeR][0][0],Tps[TypeR][0][1]]])
                else:
                    GRNInt['TF'].append([reg, [Tps[TypeA][0][0],Tps[TypeR][0][1]]] )
        return 
    def ValueGeneration():
        Species=''
        Values = {'M':[InitVals[0]], 'P':[InitVals[1]],'Dp':[InitVals[2]],'Dm':[InitVals[3]],'L':[InitVals[4]],'TRa':[InitVals[5]], 'TRb':[InitVals[6]],'Kd':[InitVals[7]], 'K':[InitVals[8]], 'H': [InitVals[9]]}
        Keys = Values.keys()
                    
            
        for i in Keys:
            for n in np.arange(1, genes+1):
                val = 0
                if i == 'M' or i=='4P':
                    Values[i].append(val)
                    Species += i + str(n) + ' = ' + str(val) + '; \n'
                else:
                    while val <= 0:
                        val = round(np.random.normal(Values[i][0],.25),3)
#                        val = np.random.normal(Values[i][0],.25)
#                        if i == 'M' or i=='P':
#                            val = round(val,0)
#                        else:
#                            val = round(val,3)
                Values[i].append(val)
                Species += i + str(n) + ' = ' + str(val) + '; \n'
        StringList.append(Species)
        return
    #%% Identify proteins regulating genes
    def ChooseProtein():
        for i in range(genes):
            TF=[]
            string = GRNInt['TF'][i][1]
            for n in string:
                TF.append(int(n.replace("P","")))
            GRNInt['TFs'].append(TF)
    #print(GRNInt['TFs'])
        return
    
    
    #%%Generate Rxn Rates
    def GetReactionRates():
        ReactionRates=[]
        WorkingString=''
        for i in range(genes):
            Reg = GRNInt['TF'][i][0]
            
            WorkingString += 'Rm_' + str(i+1) + ' := '
            TF1 = GRNInt['TFs'][i][0] - 1
            P1 = GRN['Prot'][TF1]
            K1 =  GRN['Dissoc'][TF1]
            HCoeff = GRN['HillCoeff'][i]
            num = ''
            denom = ''
            frac = ''
            if 'Single' in Reg:
                frac = '(' + P1 + '*' + K1 + ')'
                if 'Act' in Reg:
                    Hill = '(' + frac + '/(1 + ' + frac + '))'
                else:
                    Hill = '(1 /(1 + ' + frac + '))'
                WorkingString = GRN['Leak'][i] + '+' + GRN['TRa'][i] + '*' + Hill + '-' + GRN['mRNA'][i] + '*' + GRN['mDeg'][i] + '\n'
                
            elif len(GRNInt['TFs'][i]) >1 :
                TF2 = GRNInt['TFs'][i][1] - 1
                
                P2 = GRN['Prot'][TF2]
                K2 =  GRN['Dissoc'][TF2]
                K3 = GRN['DualDissoc'][TF2]
                if 'Counter' in Reg:
                    num = K1 + '*' + P1
                    denom = '(1+' + K1 + '*' + P1 + '+' + K2 + '*' + P2 + '+' + K3 + '*' + P1 + '*' + P2 + ')'#K3 is a place holder.. need to generate dissoc associated with dual reg 
                    eq = '(' + num  + '/' + denom + ')'
                            
                elif  not('N' in Reg) and("and" in Reg or 'or' in Reg):
                    if 'and' in Reg:
                        num = '(' + K1 + '*' + K2 + '*' + P1 + '*' + P2  + ')'
                        denom = '(1+' + K1 + '*' + P1 + '+' + K2 + '*' + P2 + '+' + K1 + '*' + K2 + '*' + P1 + '*' + P2 + ')'
                        eq = '(' + num  + '/' + denom + ')'
                    else:
                        num = '(' + K1 + '*' + P1 + '+' + K2 + '*' + P2 + ')'
                        denom = '(1+' + K1 + '*' + P1 + '+' + K2 + '*' + P2  + ')'
                        eq = '(' + num  + '/' + denom + ')'
                else:
                    denom = '(1+' + K1 + '*' + P1 + '+' + K2 + '*' + P2 + '+' + K3+ '*' + P1 + '*' + P2 + ')'
                    if 'and' in Reg:
                        num = '(1+' + K1 + '*' + P1 + '+' + K2 + '*' + P2 + ')'
                        eq = '(' + num  + '/' + denom + ')'
                    else:
                        eq =  '(1' + '/' +  denom + ')'
                WorkingString = GRN['Leak'][i] + '+' + GRN['TRa'][i] + '*' + eq + '-' + GRN['mRNA'][i] + '*' + GRN['mDeg'][i]+ '\n' 
        
            ReactionRates.append(WorkingString)
        return ReactionRates
    #%%Events-- attempt to perturb a user defined genes trancription at the associated time
    def GetEvents():
        EventString = ''
#        for m in range(genes):
#            EventString += 'Em' + str(m+1) + ': at(' + GRN['mRNA'][m] + '<0): ' + GRN['mRNA'][m] + '=0;\n'
#            EventString += 'Ep' + str(m+1) + ': at(' + GRN['Prot'][m] + '<0): ' + GRN['Prot'][m] + '=0;\n'    Old Catch. Previous issue.
        counter = 0
        for k in Perturbation.keys():
            counter +=1
            EventString += 'Ed' + str(counter) + ': at(time >'  + str(Perturbation[k][0]) +'): ' + k + '= ' + str(39) +';\n' 
        StringList.append(EventString)
        return
    #%% Reaction Model Generation
    def GetReactions():
        ReactionM =''
        ReactionP =''
        ReactionSumString= 'model Random_GRN()\n' 
        
        for i in np.arange(1, genes+1):
            M = GRN['mRNA'][i-1]
            P = GRN['Prot'][i-1]
            ReactionM += 'Rm' + str(i) + ':' + '$N' + '=> ' + M + ';' + Rates[i-1] + '\n'
            ReactionP += 'Rp' + str(i) + ':' + '$A' + '=> ' + P + ';' + GRN['TRb'][i-1] + '*' + M + '-' + P + '*' + GRN['PDeg'][i-1] + '\n'
        ReactionSumString += ReactionM +ReactionP
        StringList.append(ReactionSumString)
          
        
    
    def GetModelString():
        antStr = ''
        for i in StringList:
            antStr += i
        antStr += 'end'
        return antStr 
    #%%
    ModelInitNames()
    AssignProteins()
    DetermineInteractors()
    ChooseProtein()
    Rates = GetReactionRates()
    GetReactions()
    ValueGeneration()
#    GetEvents()
    antStr = GetModelString()
    return antStr
#%%
def GetModel(genes,Probabilities, InitVals, Name):
    tries = 0
    while tries < 10:
       try:
          antStr = RandomGRN(genes,Probabilities, InitVals,{'Tra1':[10]})
          r = te.loadAntimonyModel(antStr)
          print r.steadyState()
          break
       except:
           tries = tries + 1
           print "Failed to find steady state"
    return ((r, antStr))

#%%
# fileName == entire path
def writecsvFile (fileName, r, data):
    names = r.getFloatingSpeciesIds()
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
def Output(data, Name,r,rSeed,result,NoisyResult, fileName, antStr):
    for i in range(len(data)):
        if i == 0 and data[i] == 'y':
            writecsvFile(fileName + Name + '_Results.csv',r,result)
        elif data[i] == 'y':
            writecsvFile(fileName + Name + '_Noisy Result.csv',r,NoisyResult)
    if rSeed != 0:
        fh = open(fileName + Name + '_Seed.txt', 'wb')
        fh.write('Random Seed = ' + str(rSeed))
    sbmlStr = r.getSBML()
    te.saveToFile (fileName + Name + '_Model.xml', sbmlStr)
    fh = open(fileName + Name + '_Antimony.txt', 'wb')
    fh.write(str(antStr))
    plt.close('all')
# RunModel takes the number of genes, the duration of the simulation, a list of three probabilities for regulation types, the noise level, 
# the name of the simulations, the user defined responses that reflect which data will be output, a seed for reproducibile simulations,
# the number of steps for the simulation, and a string associated with the directory the main output directory for all GRN generations.
#
# Given the correct inputs, this sub routine will call and sub routine GetModel, which will result in the creation of a random GRN(if no seed input)
# The model is simulated using user defined inputs(if the model fails a new GRN is generated up to 10 times). A noisy dataset is created by 
# adding noise to the simulated dataset. Figures containing Protein and mRNA Vs time plots are created and saved for regular and noisy simulations.
# Finally all the Outputs sub-routine is called which w
def RunModel(genes,tmax,InitProb,Percent,modelName,DataOut,rSeed, Steps, filePath):
    
    if rSeed != 0:
        np.random.seed (rSeed)
    Percent /= 100
    InitParams = [0,0,0.5,0.9,0.8,30,30,0.2,0.5,1]
    r,antStr = GetModel(genes,InitProb, InitParams, modelName) 
    tries = 0
    while tries < 10:
        try:
            plt.close('all')
            r.reset()
            result = r.simulate(0,tmax, Steps)
            NoisyResult = np.zeros([len(result[:,0]), len(result[0,:])])
            for k in range(len(result[:,0])):
                for i in range(len(result[0])):
                    if i == 0:
                        NoisyResult[k,i] = result[k,i]
                    else:    
                        CurrVal = -1
                        while CurrVal < 0:
                            CurrVal = result[k,i] + np.random.normal(0,Percent*result[k,i])
                        NoisyResult[k,i] = CurrVal
            vars = {'P':[],'M':[]}
            for e in vars.keys():
                plt.figure(1)    
                for k in np.arange(1,genes+1):
                    vars[e].append(e + str(k))
                    tStep = int(math.ceil(tmax)/20)
                    tRange = np.arange(0,(int(math.ceil(tmax)))+ tStep,tStep)
                    
#                    for i in range(len((result[0,:]-1)/2)):
                    if e=='P':
                        plt.subplot(2,1,1)
                        plt.title('Protein Count Vs. Time', fontsize=8)
                        plt.grid(color='k', linestyle='-', linewidth=.4)
#                        plt.xlabel('time(s)')
                        plt.ylabel('count',fontsize=6)
                        plt.yticks(fontsize=6)
                        plt.plot (result[:,0],result[:,(k-1)+1], label = vars[e][k-1])
                        plt.yscale('log')
                        if tmax > 5000: 
                            plt.xscale('log')
                            plt.xticks(fontsize=6)
                        else:
                            plt.xticks(tRange, fontsize=6) 
                        plt.xlim(0,tmax)
                        plt.legend(loc='upper right', bbox_to_anchor=(1.13, 1.035), fontsize=6)
                    else:
                       plt.subplot(2,1,2) 
                       plt.title('mRNA Count Vs. Time', fontsize=8)
                       plt.grid(color='k', linestyle='-', linewidth=.4)
                       plt.xlabel('time(s)',fontsize=6)
                       plt.ylabel('count',fontsize=6)
                       plt.yticks(fontsize=6)
                       plt.plot (result[:,0],result[:,(k-1)+genes+1], label = vars[e][k-1])
                       plt.yscale('log')
                       if tmax > 5000: 
                            plt.xscale('log')
                            plt.xticks(fontsize=6)
                       else:
                           plt.xticks(tRange, fontsize=6) 
                       plt.xlim(0,tmax)
                       plt.legend(vars[e],loc='upper right', bbox_to_anchor=(1.13, 1.035), fontsize=6)
            manager = plt.get_current_fig_manager()
            manager.window.showMaximized()
            plt.savefig(filePath + modelName + '_Simulation Plot.png', dpi=400)             
            vars = {'P':[],'M':[]}             
            for e in vars.keys():
                plt.figure(2)    
                for k in np.arange(1,genes+1):
                    vars[e].append(e + str(k))
                    
                    if e=='P':
                        plt.subplot(2,1,1)
                        plt.title('Protein Count Vs. Time (Noisy)', fontsize=8)
                        plt.grid(color='k', linestyle='-', linewidth=.4)
                        plt.ylabel('count',fontsize=6)
                        plt.xlim(0,tmax)
                        plt.yticks(fontsize=6)
                        plt.plot (NoisyResult[:,0],NoisyResult[:,(k-1)+1], label = vars[e][k-1])
                        plt.legend(vars[e])
                        plt.yscale('log')
                        if tmax > 5000: 
                            plt.xscale('log')
                            plt.xticks(fontsize=6)
                        else:
                           plt.xticks(tRange, fontsize=6) 
                        plt.legend(loc='upper right', bbox_to_anchor=(1.13, 1.035), fontsize=6)
                    else:
                        plt.subplot(2,1,2)
                        plt.title('mRNA Count Vs. Time (Noisy)', fontsize=8)
                        plt.grid(color='k', linestyle='-', linewidth=.4)
                        plt.xlim(0,tmax)
                        plt.yticks(fontsize=6)
                        plt.xlabel('time(s)',fontsize=6)
                        plt.ylabel('count',fontsize=6)
                        plt.plot (NoisyResult[:,0],NoisyResult[:,(k-1)+genes+1], label = vars[e][k-1])
                        plt.legend(vars[e])
                        plt.yscale('log')
                        if tmax > 5000: 
                            plt.xscale('log')
                            plt.xticks(fontsize=6)
                        else:
                            plt.xticks(tRange, fontsize=6)  
                        plt.legend(loc='upper right', bbox_to_anchor=(1.13, 1.035), fontsize=6)
            manager = plt.get_current_fig_manager()
            manager.window.showMaximized()
            plt.savefig(filePath + modelName + '_Noisy Simulation Plot.png', dpi=400)  
            break
        except:
             tries = tries + 1
             print "Failed to solve Model"
             if tries ==10:
                 print('Sorry you network was not able to be simulated, please try again, or read the documentation')
    
    Output(DataOut, modelName, r, rSeed, result, NoisyResult, filePath, antStr)


try:
    initPath = os.getcwd()
    print('Your current directory is ' + initPath)
    folderName = str(raw_input('Please choose a folder for saving your generated data(press enter for default): '))
    #    folderPath = 'c:\\' + str(raw_input('Please choose a folder for saving your generated data(press enter for default): ')) + '\\' 
    forbiddenChar = ['<','>',':','"','/','\'','|','?','*']
    if any(i in folderName for i in forbiddenChar):
          folderName = str(raw_input('You seem to have included a "forbidden" character in your folder name, please try another, or simply press enter for the default: '))      
          if any(i in folderName for i in forbiddenChar):
              int('i')
    folderPath = initPath + folderName + '\\' 
    if folderName == '':
        folderPath = initPath + '\\' + 'GRN generator output\\'
        if os.path.exists(folderPath) == False:
            os.mkdir(folderPath)
    else:
        if os.path.exists(folderPath) == False:
            os.mkdir(folderPath)
    modelName = str(raw_input("Please enter the name of your network(str): ")) 
    if any(i in modelName for i in forbiddenChar):
          modelName = str(raw_input('You seem to have included a "forbidden" character in your folder name, please try another, or simply press enter for the default: '))      
          if any(i in modelName for i in forbiddenChar):
              int('i')
    filePath = folderPath + modelName + '\\'
    if os.path.exists(filePath) == True:
        response = str(raw_input('You have a model ouput with this name. Would you like to replace the files?(y/n): '))
        if response != 'n' and response!='y':
            response = str(raw_input('Your response must be a y or n. Please try again.(y/n): '))
            if response != 'n' and response!='y':
                int('i')
        if response == 'n':
            modelName = str(raw_input("Please enter a different name for your network(str): "))
            filePath = folderPath + modelName + '\\'
            os.mkdir(filePath)
#        elif response == 'y':
#            
    elif os.path.exists(filePath) == False: 
        os.mkdir(filePath)  
    genes = int(raw_input("Please enter the number of genes in your network(2<int<16): "))
    if genes <= 2:
        print('You seem to have entered two, or less genes.\n')
        genes = int(raw_input("Please enter the number of genes in your network(int): "))
        if genes <= 2:
            print('The number of genes must larger than two.\n')
            int('i')
    if genes > 15:
        print('You seem to have entered a larger number of genes than the program can handle.\n')
        genes = int(raw_input("Please enter the number of genes in your network(2<int<15): "))
        if genes > 15:
            int('i')
#    elif genes:
    InitProb=[]
    InitProb.append(float(raw_input("Please enter the probability of Single regulation.(float): ")))
    InitProb.append(float(raw_input("Please enter the probability of Double regulation.(float): ")))
    CountProb = round(1-np.sum(InitProb),3)
    print('The probabilty of Counter regulation will be ' + str(CountProb) + '\n' )
    InitProb.append(CountProb)
    tries = 0
    while any(i<0 for i in InitProb) == True or round(np.sum(InitProb),3) !=1.0:
        if tries==2:
            int('i')
        else:   
            print("Sorry, the sum of the probabilities must equal 1 or at least one is negative.")
            InitProb[0] = float(raw_input("Please enter the probability (<1) of having single regulated gene (float): "))
            InitProb[1] = float(raw_input("Please enter the probability (<1) of double regulation.(float): "))
            CountProb = round(1-np.sum(InitProb[0:2]),3)
            InitProb[2] = CountProb
            print('\n Your new probabilty of counter regulation is ' + str(CountProb) + '\n' )
            tries = tries + 1
    tMax = float(raw_input("Please enter how long you'd like to simulate the network for(float): "))
    if tMax <=0:
        tMax = float(raw_input('Your simulation time must be positive, please enter a positive value.'))
        if tMax <=0:
            int('i')
    Steps = int(raw_input("Please enter number of steps to generate for your simulation output (int): "))
    if Steps <= 0 or Steps >=50:
        Steps = int(raw_input("You entered a number of steps that is outside of the allowable range, please enter an integer between 0 and 50: "))
        if Steps <= 0 or Steps >=50:
            int('i')
    NLevel = float(raw_input("Please enter the percentage of the noise level (float): "))
    if NLevel < 0 or NLevel >100:
        NLevel = float(raw_input("You entered a number of steps that is outside of the allowable range, please enter a number between 0 and 100: "))
        if NLevel < 0 or NLevel >100:
            int('i')
    DataOut = []
    DataOut.append(str(raw_input("Would you like simulated data exported(y/n): ")))
    if DataOut[0] != 'n' and DataOut[0]!='y':
            DataOut[0] = str(raw_input('Your response must be a y or n. Please try again.(y/n)'))
            if DataOut[0] != 'n' and DataOut[0]!='y':
                int('i')
    if NLevel < 0.01:
        DataOut.append ('n')
    else:
        DataOut.append(str(raw_input("Would you like noisy simulated data exported(y/n): ")))
        if DataOut[0] != 'n' and DataOut[0]!='y':
            DataOut[0] = str(raw_input('Your response must be a y or n. Please try again.(y/n)'))
            if DataOut[0] != 'n' and DataOut[0]!='y':
                int('i')    
    rSeed = int (raw_input ('Please enter an optional random seed (enter 0 for no seed): '))
    if rSeed <0 or rSeed > 2**32 - 1:
        rSeed = int (raw_input ('Your response was outside of the allowable range, please enter an integer between 0 and 2**32 - 1: '))
        if rSeed <0 or rSeed > 2**32 - 1:
            int('i')
    RunModel(genes,tMax,InitProb,NLevel,modelName,DataOut,rSeed, Steps, filePath)

except WindowsError:
    print('folder error') 
        
except ValueError:
    print("There was an issue with your last entry, check your formatting.")
except:
    print ("Somethng bad happenned")
    
#RunModel(15,5000.0,[.4,.4,.2], 10.0,'foo',['y','y'],15,400,'c:\\GRN generator output\\foo\\'  )