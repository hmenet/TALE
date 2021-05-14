A python rewrite of ALE part of a three level method in progress

# requirement : 
python 3

# Usage : 

input : lower tree and upper tree

output : reconciliation scenarios in recphyloxml, see https://github.com/simonpenel/rectree2svg to get a visual representation in svg, and event frequencies.


	python3 src/main.py upper_dir lower_dir
	
Example, command line from the repo :

	python3 src/main.py ex_pylori/pylori_species_tree/ ex_pylori/genes_sub2/ -mpf

Usage: 


	python3 src/main.py [-h] [-mdir MATCHING_DIR] [-mfile MATCHING_FILE] [-o OUTPUT][-of OUTPUT_FREQ] [-ns N_REC_SAMPLE][-nre N_RATES_ESTIMATION_STEPS][-nres N_RATES_ESTIMATION_REC_SAMPLE] [-b][-nrxml N_RECPHYLOXML] [-mp] [-mpf] [-tl THIRD_UPPER_LEVEL][-imd INTER_MATCH_DIR] [-imf INTER_MATCH_FILE] [-v] 

	positional arguments:
  	upper_dir             directory with newick rooted files for each upper level tree
  	lower_dir             directory with newick unrooted files for each lower level tree, possibility of list of trees to use amalgamation

	optional arguments:
	  -h, --help            show this help message and exit
	  -mdir MATCHING_DIR, --matching_dir MATCHING_DIR
		                provide a directory with matching files for each lower
		                level tree to match their leaves to upper tree ones.
		                Without it, lower tree and upper tree leaves must have
		                the same name. A lower leaf can be matched to multiple
		                upper leaves, which will be interpreted as uncertainty
		                on the match (not failure to diverge) (default: None)
	  -mfile MATCHING_FILE, --matching_file MATCHING_FILE
		                provide a file with all matching infos for each lower
		                level tree to match their leaves to upper tree ones.
		                Alternative to -mdir (default: None)
	  -o OUTPUT, --output OUTPUT
		                output recphyloxml file. If multiple samples, multiple
		                recphyloxml are generated, if best rec option, then
		                first rec is the best one (default: output/rec)
	  -of OUTPUT_FREQ, --output_freq OUTPUT_FREQ
		                output event frequency file. (default:
		                output/event_frequency)
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
		                scenario. (default: False)
	  -nrxml N_RECPHYLOXML, --n_recphyloxml N_RECPHYLOXML
		                number of sampled scenarios stored as recphyloxml
		                file. (default: 1)
	  -mp, --multiprocess   enable multiprocessing for the sampling parts
		                (default: False)
	  -mpf, --multiprocess_fam
		                enable multiprocessing with one process for each lower
		                tree family. If chosen, stop from using multiprocess
		                for sampling -mp (default: False)
	  -tl THIRD_UPPER_LEVEL, --third_upper_level THIRD_UPPER_LEVEL
		                add a directory with an upper level on top of the two
		                previous ones, the upper become intermediate (default:
		                None)
	  -imd INTER_MATCH_DIR, --inter_match_dir INTER_MATCH_DIR
		                add a directory for the upper intermediate matching,
		                same format as the default levels matchings (default:
		                None)
	  -imf INTER_MATCH_FILE, --inter_match_file INTER_MATCH_FILE
		                add a file for the upper intermediate matching, same
		                format as for the default levels matchings (default:
		                None)
	  -v, --verbose         increase output verbosity (default: False)




# todo :

+ save amalgamation
+ more verbose
+ commenting and renaming pass (for upper/intermediate/lower vocab) 
+ test 3 level
+ 3 level rates inference
+ 3 level: P transfer in E computation
+ possibility to match from lower to upper and erwinia example


# Some included examples :

+ ex_pylori : an example with 1 gene tree, and 1 strains tree, and with gene tree and population tree with matching file and multiple match for one leaf, from Alexia Nguyen Trung
		
		python3 src/main.py example_data/ex_pylori/pylori_species_tree/ example_data/ex_pylori/genes_sub2/
		
		python3 src/main.py example_data/ex_pylori/pop_pylori/ example_data/ex_pylori/genes_sub2/ -mfile matching_pop_genes

+ ex_pylori 3 level : 

		python3 src/main.py example_data/ex_pylori/pop_pylori/ example_data/ex_pylori/genes_sub2/ -nre 0 -mfile example_data/ex_pylori/matching_pop_genes -tl example_data/ex_pylori/upper_level/ -imf example_data/ex_pylori/matching_pop_geo -mp
		
+ ex_ALE : example given in ALE git (https://github.com/ssolo/ALE), with amalgamation

		python3 src/main.py example_data/ex_ALE/species/ example_data/ex_ALE/gene
	
+ ex_erwinia : example with a matching file for leaves, from Manzano Marin et al. https://doi.org/10.1038/s41396-019-0533-6
		
		python3 src/main.py example_data/ex_erwinia/symbiont/ example_data/ex_erwinia/genes/ -mdir example_data/ex_erwinia/lower_matching/
