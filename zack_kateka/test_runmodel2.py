from GetModel import get_model
from RunModel_2 import run_model2 
from RunModel import run_model

#seed = 19443232

seed = 123468
ant_str, biotap = get_model(8, seed=seed,export=True)

#f = open("debug_biotap.csv", 'w')
#f.write(biotap)
#f.close()
#
#run_model2(ant_str, 0.05,seed=seed, species_type='P', species_nums = [1,2,3,4,5,6,7,8], timepoints = [200,5],
#        showTimePlots=True)
#
#run_model2(ant_str, 0.05,seed=seed, species_type='M', species_nums = [1,2,3,4,5,6,7,8], timepoints = [200,5],
#        showTimePlots=True)
#run_model2(ant_str, 8, 0.05,seed=seed, species_type='M', species_nums = [1,2,3,4,5,6,7,8], timepoints = [200,5],
#        showTimePlots=True, perturbs= [("DOWN", [3]), ("KO", [5,6,7])])

#RunModel.py
#run_model(ant_str, 0.05,seed=seed, exportData=[ [1,2,3,4,5,6,7,8], 'M', True, True, True], showTimePlots=True, inputData=[1,200,5, [0],[0]])
