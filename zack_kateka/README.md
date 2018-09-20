# Function Descriptions
*(see documentation in files for more details + use instructions)

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
	contains the method 'convert_biotapestry_to_antimony' which is a convenience method
	for converting between the two formats


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
