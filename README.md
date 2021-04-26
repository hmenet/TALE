A python rewrite of ALE, easy to use, part of a three level method in progress

requirement : python 3

How to : 

	python3 main.py upper_dir lower_dir
	
Example, command line from the repo :

	python3 main.py "ex_pylori/pylori_species_tree/" "ex_pylori/genes_sub2/" -nre 1 -mp

Usage: main.py [-h] [-mdir MATCHING DIR] [-o OUTPUT] [-ns N_REC_SAMPLE]
               [-nre N_RATES_ESTIMATION_STEPS]
               [-nres N_RATES_ESTIMATION_REC_SAMPLE] [-b] [-v]
               upper_dir lower_dir

positional arguments:
  upper_dir             directory with newick unrooted files for each lower
                        level tree
  lower_dir             directory with newick unrooted files for each lower
                        level tree, possibility of list of trees to use
                        amalgamation

optional arguments:
  -h, --help            show this help message and exit
  -mdir MATCHING DIR, --matching dir MATCHING DIR
                        provide a directory with matching files for each lower
                        level tree to match their leaves to upper tree ones.
                        Without it, lower tree and upper tree leaves must have
                        the same name (default: None)
  -o OUTPUT, --output OUTPUT
                        output recphyloxml file (default: None)
  -ns N_REC_SAMPLE, --n_rec_sample N_REC_SAMPLE
                        number of reconciliation scenarios sampled to compute
                        events frequency (default: 100)
  -nre N_RATES_ESTIMATION_STEPS, --n_rates_estimation_steps N_RATES_ESTIMATION_STEPS
                        number of steps in the rates estimation process, for
                        each step we compute likelihood of reconciliation and
                        set rates to the observed frequency for each events
                        (default: 5)
  -nres N_RATES_ESTIMATION_REC_SAMPLE, --n_rates_estimation_rec_sample N_RATES_ESTIMATION_REC_SAMPLE
                        in the rates estimation process, number of scenarios
                        sampled to estimate event frequencies (default: 100)
  -b, --best_rec        return the best (maximum likelihood) reconciliation
                        scenario. Require output file specified with -o option
                        (default: False)
  -v, --verbose         increase output verbosity (default: False)


todo :

	multiprocessing : parallel for gene families
	
	save amalgamation
	
	more verbose
	
	possibility for one lower leaf to be match to multiple upper ones with uniform prior
	
	add 3 level


Some included examples :
	ex_pylori : an example with 1 gene tree, and 1 species tree
		
	ex_ALE : example given in ALE git, with amalgamation
	
	ex_erwinia : example with a matching file for leaves, from Manzano Marin et al.


