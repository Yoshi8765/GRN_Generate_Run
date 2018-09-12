# -*- coding: utf-8 -*-
"""
Created on Wed Aug 01 16:38:48 2018

@author: Yoshi
"""

def BiFanGenerator():
    import tellurium as te
    import numpy as np
    import os
    #import matplotlib.pyplot as plt
    
    np.set_printoptions(linewidth=160)
    
    #%%
    # generalized_hill := Vm_A*( (K_A*p^H) / (K_A + p^H) ) #Activation# + 
    #                     Vm_R*( K_R / (K_R + p^H) ) #Repression# ;
    #                     To turn either/or off, set Vm_A or Vm_R to 0.
    
    r=te.loada('''
    model *BiFanMotif()
    
      // Compartments and Species:
      const N, AA; # Nucleus, Amino Acids //must specify how many genes for randomizer code to work!
      species m1, m2, m3, m4;
      species p_i, p2, p3, p_o;
      
      // Assignment Rules:
      //hill1: Regulation of p_i to p3
      //hill2: Regulation of p2 to p3
      //hill3: Regulation of p_i to p_o
      //hill4: Regulation of p2 to p_o
      
      // Reactions:
      ts1: => m1        ; L1 + a_m1 - d_m1*m1
      ts2: => m2        ; L2 + a_m2 - d_m2*m2
      ts3: => m3        ; L3 + Vm_A1*( (K_A1*p_i^H1) / (K_A1 + p_i^H1) ) + Vm_R1*( K_R1 / (K_R1 + p_i^H1) ) * Vm_A2*( (K_A2*p2^H2) / (K_A2 + p2^H2) ) + Vm_R2*( K_R2 / (K_R2 + p2^H2) ) - d_m3*m3
      ts4: => m4        ; L4 + Vm_A2*( (K_A2*p2^H2) / (K_A2 + p2^H2) ) + Vm_R2*( K_R2 / (K_R2 + p2^H2) ) * Vm_A2*( (K_A2*p2^H2) / (K_A2 + p2^H2) ) + Vm_R2*( K_R2 / (K_R2 + p2^H2) ) - d_m4*m4
      tl1: => p_i  ; a_p1*m1 - d_p1*p_i
      tl2: => p2       ; a_p2*m2 - d_p2*p2
      tl3: => p3       ; a_p3*m3 - d_p3*p3
      tl4: => p_o ; a_p4*m4 - d_p4*p_o
      // Species initializations:
      N = 1;
      AA = 1;
      m1 = 0;
      m2 = 0;
      m3 = 0;
      m4 = 0;
      p_i = 0;
      p2 = 0;
      p3 = 0;
      p_o = 0;
      
      // Parameter initializations:
         L1 = .01;    L2 = .01;  L3 = .01;  L4 = .01;
       K_A1 = .65;  K_A2 = .65; 
       K_R1 = .65;  K_R2 = .65; 
      Vm_A1 =  15; Vm_A2 =  15; 
      Vm_R1 =  15; Vm_R2 =  15; 
        d_m1 = .5;   d_m2 = .5; d_m3 = .5; d_m4 = .5;
        d_p1 = .5;   d_p2 = .5; d_p3 = .5; d_p4 = .5;
        a_m1 = 15;   a_m2 = 15;
        a_p1 = .5;   a_p2 = .5; a_p3 = .5; a_p4 = .5;
          H1 =  1;     H2 =  1;  
    end
    ''')
    Params = r.getGlobalParameterIds()[:-r.getNumFloatingSpecies()/2] 
    # we don't need to randomize var objects (hill expressions) so take that out for Params
    
    # Relative abundance of FFL types in yeast + e. coli. 
    # Numbers totaled together between the two organisms.
    # From Mangan & Alon (2003,PNAS), Table 1+2.
    
    #X->Y , Y->Z, X->Z : Total abundance = 98
    FFL_types = {
            #coherent types
            "+++": 54/98.0,
            "--+": 1/98.0,
            "-+-": 7/98.0,
            "+--": 4/98.0,
            #incoherent types
            "+-+": 26/98.0,
            "-++": 1/98.0,
            "++-": 2/98.0,
            "---": 3/98.0,
            }
    types = FFL_types.keys()
    freq = np.array(FFL_types.values())
    
    picks = np.random.choice(types,1,p=freq)
    
    ##Check if picks is working
    #a = []
    #import collections
    #for n in np.arange(10000):
    #    a.extend(np.random.choice(types,1,p=freq))
    #typesFreq = collections.Counter(a)
    
    counter = 0
    for n in range(len(Params)):
        param = Params[n]
        randVal = 0
        while randVal <= 0:
            val = r.getValue(param)
            randVal = np.random.normal(val,val*.25)
            if picks[counter]=='+':
                if param[0:3]=='Vm_A':
                    randVal = 1
                if param[0:3]=='Vm_R':
                    randVal = 0
                counter += 1
            if picks[counter]=='-':
                if param[0:3]=='Vm_A':
                    randVal = 0
                if param[0:3]=='Vm_R':
                    randVal = 1
                counter += 1            
            else:
                randVal = round(randVal,4)
                    
            setattr(r, param,randVal)
    
#    tmax=200
    
#    result = r.simulate(0, tmax, tmax*2,)
    
#    plt.figure()
#    plt.grid(color='k', linestyle='-', linewidth=.4)
#    plt.ylim(0,np.max(result[:,4:7])*1.1)
#    plt.xlim(0,tmax)
#    plt.yticks(np.arange(0,np.max(result[:,4:7])*1.1,np.max(result[:,4:7])/12))
#    #M1 , = plt.plot (result[:,0],result[:,1], label = 'M1')
#    #M2 , = plt.plot (result[:,0],result[:,3], label = 'M2')
#    #M2 , = plt.plot (result[:,0],result[:,6], label = 'M3')
#    p_i , = plt.plot (result[:,0],result[:,4], label = 'p_i')
#    P2 , = plt.plot (result[:,0],result[:,5], label = 'P2')
#    p_o , = plt.plot (result[:,0],result[:,6], label = 'p_o')
#    plt.legend([p_i, P2, p_o], ['p_i', 'P2', 'p_o'])
        
    r.reset()
    #plt.close("all")
    #res = r.simulate(0,50,1000)
    #r.plot()
    #r.draw()
    #r.reset()
    print('Saving model...\n')
    r.exportToAntimony('BiFan.txt') #export as antimony
    return str(os.getcwd())+('\\BiFan.txt')