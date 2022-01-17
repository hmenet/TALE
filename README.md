Three level Amalgamated Likelihood Estimation (TALE) is a three-level Duplication Transfer Loss (DTL) phylogenetic reconciliation method integrating gene, symbiont and host in a single framework. It is based on a reimplementation of ALE undated (https://github.com/ssolo/ALE) in python 3, requiring basic python packages. 

Duplication Transfer Loss (DTL) phylogenetic reconciliation is a set of method to compare the evolution of two levels, such as host and symbiont or gene and species, through three specific evolutionary events, duplication, horizontal transfer and loss. 



# Input
TALE takes as input three sets of trees in newick format, and without taking into account branch lengths:
+ host : binary and rooted
+ symbiont : binary, rooted or unrooted, or a distribution of binary trees to account for uncertainty, to be amalgamated
+ gene : binary and unrooted, or a distribution of binary trees to account for uncertainty, to be amalgamated


# Requirement : 
Python 3 with numpy, random, os, time, argparse, multiprocessing, sys

# Usage : 

## 2-level :

Input : a directory with a (or a set of) binary rooted upper tree (host, species) in newick format, and a directory with a lower tree (symbiont, or gene) in newick format

Output : reconciliation scenarios in recphyloxml, see https://github.com/simonpenel/rectree2svg to get a visual representation in svg, and event frequencies for each gene families.


	python3 src/main.py upper_dir lower_dir

## 3-level : 

	python3 src/main.py symbiont_dir gene_dir -tl host_dir



		usage: main.py [-h] [-mdir MATCHING_DIR] [-mfile MATCHING_FILE] [-o OUTPUT]
		       [-of OUTPUT_FREQ] [-ns N_REC_SAMPLE]
		       [-nre N_RATES_ESTIMATION_STEPS]
		       [-nres N_RATES_ESTIMATION_REC_SAMPLE] [-b]
		       [-nrxml N_RECPHYLOXML] [-mp] [-mpf] [-tl THIRD_UPPER_LEVEL]
		       [-imd INTER_MATCH_DIR] [-imf INTER_MATCH_FILE]
		       [-tlh THREE_LEVEL_HEURISTIC] [-tlMCs THREE_LEVEL_MC_SAMPLE]
		       [-inre INTER_N_RATES_ESTIMATION_STEPS]
		       [-inres INTER_N_RATES_ESTIMATION_REC_SAMPLE] [-v]
		       [-dr DUPLICATION_RATE] [-tr TRANSFER_RATE] [-lr LOSS_RATE]
		       upper_dir lower_dir

	positional arguments:
	  upper_dir             directory with newick unrooted files for each lower
		                level tree
	  lower_dir             directory with newick unrooted files for each lower
		                level tree, possibility of list of trees to use
		                amalgamation

	optional arguments:
	  -h, --help            show this help message and exit
	  -mdir MATCHING_DIR, --matching_dir MATCHING_DIR
		                provide a directory with matching files for each
		                lower level tree to match their leaves to upper tree
		                ones. Without it, lower tree and upper tree leaves
		                must have the same name. A lower leaf can be matched
		                to multiple upper leaves, which will be interpreted
		                as uncertainty on the match (not failure to diverge)
		                (default: None)
	  -mfile MATCHING_FILE, --matching_file MATCHING_FILE
		                provide a file with all matching infos for each lower
		                level tree to match their leaves to upper tree ones.
		                Alternative to -mdir (default: None)
	  -o OUTPUT, --output OUTPUT
		                output recphyloxml file. If multiple samples,
		                multiple recphyloxml are generated, if best rec
		                option, then first rec is the best one (default:
		                output/rec)
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
		                enable multiprocessing with one process for each
		                lower tree family. If chosen, stop from using
		                multiprocess for sampling -mp (default: False)
	  -tl THIRD_UPPER_LEVEL, --third_upper_level THIRD_UPPER_LEVEL
		                add a directory with an upper level on top of the two
		                previous ones, the upper become intermediate
		                (default: None)
	  -imd INTER_MATCH_DIR, --inter_match_dir INTER_MATCH_DIR
		                add a directory for the upper intermediate matching,
		                same format as the default levels matchings (default:
		                None)
	  -imf INTER_MATCH_FILE, --inter_match_file INTER_MATCH_FILE
		                add a file for the upper intermediate matching, same
		                format as for the default levels matchings (default:
		                None)
	  -tlh THREE_LEVEL_HEURISTIC, --three_level_heuristic THREE_LEVEL_HEURISTIC
		                heuristic for 3 level rec, can be either dec for
		                decoupled or MC for montecarlo (default: MC)
	  -tlMCs THREE_LEVEL_MC_SAMPLE, --three_level_MC_sample THREE_LEVEL_MC_SAMPLE
		                number of samples of upper intermediate
		                reconciliation for monte carlo heuristic of 3 level
		                reconciliation (default: 10)
	  -inre INTER_N_RATES_ESTIMATION_STEPS, --inter_n_rates_estimation_steps INTER_N_RATES_ESTIMATION_STEPS
		                number of steps in the rates estimation process for
		                the inter upper rec, for each step we compute
		                likelihood of reconciliation and set rates to the
		                observed frequency for each events (default: 5)
	  -inres INTER_N_RATES_ESTIMATION_REC_SAMPLE, --inter_n_rates_estimation_rec_sample INTER_N_RATES_ESTIMATION_REC_SAMPLE
		                in the rates estimation process of inter and upper,
		                number of scenarios sampled to estimate event
		                frequencies (default: 100)
	  -v, --verbose         increase output verbosity (default: False)
	  -dr DUPLICATION_RATE, --duplication_rate DUPLICATION_RATE
		                initial duplication rate. (default: 0.01)
	  -tr TRANSFER_RATE, --transfer_rate TRANSFER_RATE
		                initial transfer rate. (default: 0.01)
	  -lr LOSS_RATE, --loss_rate LOSS_RATE
		                initial loss rate. (default: 0.01)



# Included examples

## 2-level examples

+ ex_pylori : an example with 1 gene tree, and 1 strains tree, and with gene tree and population tree with matching file and multiple match for one leaf, from Alexia Nguyen Trung
		
		python3 src/main.py example_data/ex_pylori/pylori_species_tree/ example_data/ex_pylori/genes_sub2/
		
		python3 src/main.py example_data/ex_pylori/pop_pylori/ example_data/ex_pylori/genes_sub2/ -mfile matching_pop_genes
+ ex_ALE : example given in ALE git (https://github.com/ssolo/ALE), with amalgamation

		python3 src/main.py example_data/ex_ALE/species/ example_data/ex_ALE/gene

+ ex_erwinia : example with a matching directory for leaves, from Manzano Marin et al. https://doi.org/10.1038/s41396-019-0533-6
		
		python3 src/main.py example_data/ex_erwinia/symbiont/ example_data/ex_erwinia/genes/ -mdir example_data/ex_erwinia/lower_matching/
		
## 3-level examples

+ ex_pylori 3 level : 

		python3 src/main.py example_data/ex_pylori/pop_pylori/ example_data/ex_pylori/genes_sub2/ -nre 0 -mfile example_data/ex_pylori/matching_pop_genes -tl example_data/ex_pylori/upper_level/ -imf example_data/ex_pylori/matching_pop_geo -mpf

		
+ ex_erwinia 3-Level : example with a matching directory for leaves, and matching file for upper level and free living symbionts, from Manzano Marin et al. https://doi.org/10.1038/s41396-019-0533-6
		
		python3 src/main.py example_data/ex_erwinia/symbiont/ example_data/ex_erwinia/genes/ -mdir example_data/ex_erwinia/lower_matching/ -tl example_data/ex_erwinia/host/ -imf example_data/ex_erwinia/matching_symbiont_host -tlh dec -mpf -ncpu 4
