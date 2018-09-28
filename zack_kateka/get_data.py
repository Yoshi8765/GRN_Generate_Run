from run_experiments import export_experiments
import os
import numpy as np


"""
This section will create a team score file if it doesn't exist.
***Be sure to change team names in the `team_names` array.***
"""
if os.path.isfile("team_scores.csv") == False:
    # change team names here. Make sure it matches the google form team names
    team_names = np.array(["team 1", "team 2", "team 3", "team 4", "team 5", "team 6"])
    money = np.array(["14000"]*len(team_names))
    result = np.vstack((team_names, money))
    np.savetxt("team_scores.csv", result, fmt="%s",delimiter=",")
    print("Team score file made!")

"""
This file will send out the next batch of experimental data.
Simply ensure that num_genes is the number of genes in the network, and run the script.
"""
num_genes = 8

export_experiments(num_genes, sendEmail=True, updateMoney=True)
