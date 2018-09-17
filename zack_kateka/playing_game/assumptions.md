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


COSTS
* RNA seq (low/high res)			900/1300
* Mass Spec (low/high res)			1000/1500
* GFP Protein (up to 3)				200/500/850

* Full gene knockout (Gene Deletion)		1000 
* Partial gene knockout (CRISPR-based)		350/550
	* arbitrary (20%-50%)
	* choice (20%,30%,40%,50%)
* Partial gene upregulation (CRISPR-based)	350/550
        * arbitrary (20%-50%)
        * choice (20%,30%,40%,50%)
