* You will be given rough, biologically reasonable values for each parameter
* There are no outside players; all genes involved are already present
* There can be at most two regulators for each gene
* A single protein cannot act on both regulatory sites of another gene
* The number of regulatory sites a protein can act on is not limited
* INPUT does not act on any genes besides the one it is currently attached to
* A protein might act on nothing
* self-feedback loops are allowed
* The true network is fully connected
* Nothing can act on INPUT
* Your experimental data will have error
* All of the currently provided connections are correct
* INPUT is constant at 1 unit



ADVICE
* Biotapestry has a convenient CSV format for building models (see the tutorials on their website
	for more information). If you write code to convert this CSV format to an antimony string,
	you will save yourself a lot of time trying to simulate potential connections. We also wrote
	some helper methods that can quickly add/remove connections from the CSV format.

