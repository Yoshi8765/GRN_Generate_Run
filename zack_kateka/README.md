# Usage Instructions
* NOTE: Do not rearrange the provided files; they are required to have a specific relative path to each other, so
  rearranging them will break the code. Feel free to place the overall folder wherever you wish however.

## Generating A Model
* In order to generate a model, you will run GetModel.py. I have created a helper file called make_model.py that makes this process easy.
  Open this file, and there will be some comments describing some options you have to set certain features of the model. When you are satisfied
  with these, run the script. You might have to run it multiple times to generate a model that is sufficiently interesting. This will create a
  biotapestry file and antimony file in your current working directory. Hold on to these.

## Breaking the model
* In order to break the model, first open Biotapestry. From there, click File > Import > Import Full Model Hierarchy from CSV and select the
  biotapestry CSV generated in the previous step. This will give you a visualization of the model, and from here select which connections to remove.
  Deleting the connections off the Biotapestry file will not work, as Biotapestry does not support this. As such, we have written our own code for
  removing/adding connections, found in change_biotapestry.py. I have written a helper file called break_model.py that includes instructions for 
  breaking the model. This will output a new CSV of the form model_name_broken.csv. I recommend this to be the file you give to the students.

## Ordering Experimental Data
* To collect data orders from students, we have created a google form [BIOEN 498: Experiment Request Form](https://docs.google.com/forms/d/1OFsoRf8hEJw4d3bpdQHlR1wrq_fUVGD6PmKRf3d1TdY/edit) which students should be given access to.
  







# Function Descriptions
*(see documentation in files for more details)

## Main Folder

### GetModel.py
	Allows you to randomly generate a gene regulatory network that meets certain criteria.
	Has functions for converting this network to an antimony string or to a 
	CSV format that Biotapestry can read. Returns a tuple (antimony_string, biotapestry_string)
	and can also has an option to export this information to files in the working directory
	Relavent functions:
		get_model()
		convert_to_antimony()
		convert_to_biotapestry()

### Biotapestry.py
	Contains the method 'convert_biotapestry_to_antimony' which is a convenience method
	for converting between the two formats. Since biotapestry does not store parameter values,
	you must provide these manually. 

### RunModel.py
	Allows you to run a model (generated using GetModel.py) and generate "experimental"
	data from it. Can add noise to this data as well.

### change_biotapestry.py
	Has methods for automatically adding and removing connections in a gene network
	from a Biotapestry CSV format. Convenient for breaking gene network, or for
	trying out new possible connections in attempt to fix the broken network.
	Relevant Functions:
		add_biotapestry() = adds the given connections to model
		remove_biotapestry() = removes the given connections from model

### compare_biotapestry.py
	Compares two biotapestry CSV formats and outputs their differences.
	Might be useful for assessing how well students captured the true network
	at end of quarter

### run_experiments.py
	Given the csv from google forms, will parse through and run the correct experiments for each entry. Will update the team's money and send email with the csv of the experiment results to the student who filled the form.

### Other
	experiments.md = list of available experiments and pertubations, as well as their costs


## Playing the Game folder

### data_analysis.py
	Used to make sense of experimental data. This includes parameter estimation, and probing
	for possible missing connections.

### the_game.py
	Plots current working model vs experimental data to help spot shortcomings in 
	current model. Helps us manually decide which connections to consider, and
	which tests to order

### Other
	assumptions.md = list of assumptions made about the true network when playing the game
			We tried to capture all the information we felt the students might need
			to fairly play the game
	what we think.md = describes the strategy we used while playing the game
	model_files = all the experimental data we bought
