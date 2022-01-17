Three level Amalgamated Likelihood Estimation (TALE) is a three-level Duplication Transfer Loss (DTL) phylogenetic reconciliation method integrating gene, symbiont and host in a single framework. It is based on a reimplementation of ALE undated (https://github.com/ssolo/ALE) in python 3, requiring basic python packages, with a command line interface and recphyloxml output.

Duplication Transfer Loss (DTL) phylogenetic reconciliation is a set of method to compare the evolution of two levels, such as host and symbiont or gene and species, through three specific evolutionary events, duplication, horizontal transfer and loss. 


# Input
TALE takes as input three sets of trees in newick format, and without taking into account branch lengths:
+ host : binary and rooted
+ symbiont : binary, rooted or unrooted, or a distribution of binary trees to account for uncertainty, to be amalgamated
+ gene : binary and unrooted, or a distribution of binary trees to account for uncertainty, to be amalgamated

The matching between the leaves of the gene and symbiont tree and between the host and symbiont one can be given in different fashion with directory containing files for each gene family or with a file valid for every genes. By default the program will try to match a leaf to the leaf with the same name in the other tree. See example data.

# Output 

Multiple files are given as output

+ Recphyloxml files are reconciliation scenarios, in recphyloxml format, see Thirdkind (https://github.com/simonpenel/thirdkind) to get a visual representation in svg. In 3-level, upper correspond to the host/symbiont reconciliation and lower to the symbiont/gene one.
+ frequency files are the observed frequencies of evolutionary events over multiple sampled scenarios for each gene family
+ newick files are trees with renamed internal nodes (corresponding to the frequency files) or the trees obtained through amalgamation
+ additional info give information on the reconciliation and reconciliation scenario such as likelihood.

# Requirement : 
Python 3 with numpy, random, os, time, argparse, multiprocessing, sys

# Usage : 

## 2-level :


	python3 src/main.py upper_dir lower_dir

## 3-level : 

	python3 src/main.py symbiont_dir gene_dir -tl host_dir

## Options :

	usage: main.py [-h] [-mdir MATCHING_DIR] [-mf MATCHING_FILE] [-o OUTPUT] [-ns N_REC_SAMPLE]
		       [-nre N_RATES_ESTIMATION_STEPS] [-nres N_RATES_ESTIMATION_REC_SAMPLE] [-b]
		       [-nrxml N_RECPHYLOXML] [-inrxml INTER_N_RECPHYLOXML] [-mpf]
		       [-ncpu N_CPU_MULTIPROCESS] [-tl THIRD_UPPER_LEVEL] [-imd INTER_MATCH_DIR]
		       [-imf INTER_MATCH_FILE] [-tlh THREE_LEVEL_HEURISTIC]
		       [-tlMCs THREE_LEVEL_MC_SAMPLE] [-inre INTER_N_RATES_ESTIMATION_STEPS]
		       [-inres INTER_N_RATES_ESTIMATION_REC_SAMPLE] [-ins INTER_N_REC_SAMPLE]
		       [-insr INTER_N_REC_SAMPLE_RATES] [-v] [-dr DUPLICATION_RATE]
		       [-tr TRANSFER_RATE] [-lr LOSS_RATE] [-idr INTER_DUPLICATION_RATE]
		       [-itr INTER_TRANSFER_RATE] [-ilr INTER_LOSS_RATE] [-ilo] [-ia] [-dd]
		       [-ir INCOMPLETE_SORTING_RATE] [-iir INTER_INCOMPLETE_SORTING_RATE] [-geo]
		       [-slm SECOND_LEVEL_MODEL]
		       upper_dir lower_dir

	positional arguments:
	  upper_dir             directory with newick unrooted files for each lower level tree
	  lower_dir             directory with newick unrooted files for each lower level tree,
		                possibility of list of trees to use amalgamation

	optional arguments:
	  -h, --help            show this help message and exit
	  -mdir MATCHING_DIR, --matching_dir MATCHING_DIR
		                provide a directory with matching files for each lower level tree to
		                match their leaves to upper tree ones. Without it, lower tree and
		                upper tree leaves must have the same name. A lower leaf can be
		                matched to multiple upper leaves, which will be interpreted as
		                uncertainty on the match (not failure to diverge) (default: None)
	  -mf MATCHING_FILE, --matching_file MATCHING_FILE
		                provide a file with all matching infos for each lower level tree to
		                match their leaves to upper tree ones. Alternative to -mdir (default:
		                None)
	  -o OUTPUT, --output OUTPUT
		                output recphyloxml file. If multiple samples, multiple recphyloxml
		                are generated, if best rec option, then first rec is the best one
		                (default: output/rec)
	  -ns N_REC_SAMPLE, --n_rec_sample N_REC_SAMPLE
		                number of reconciliation scenarios sampled to compute events
		                frequency (default: 100)
	  -nre N_RATES_ESTIMATION_STEPS, --n_rates_estimation_steps N_RATES_ESTIMATION_STEPS
		                number of steps in the rates estimation process, for each step we
		                compute likelihood of reconciliation and set rates to the observed
		                frequency for each events (default: 5)
	  -nres N_RATES_ESTIMATION_REC_SAMPLE, --n_rates_estimation_rec_sample N_RATES_ESTIMATION_REC_SAMPLE
		                in the rates estimation process, number of scenarios sampled to
		                estimate event frequencies (default: 100)
	  -b, --best_rec        return the best (maximum likelihood) reconciliation scenario.
		                (default: False)
	  -nrxml N_RECPHYLOXML, --n_recphyloxml N_RECPHYLOXML
		                number of sampled scenarios stored as recphyloxml file (if multiple
		                upper scenario sampled, number of lower scenarios for each upper
		                one). (default: 1)
	  -inrxml INTER_N_RECPHYLOXML, --inter_n_recphyloxml INTER_N_RECPHYLOXML
		                number of sampled inter upper scenarios stored as recphyloxml file.
		                (default: 1)
	  -mpf, --multiprocess_fam
		                enable multiprocessing with one process for each lower tree family.
		                If chosen, stop from using multiprocess for sampling -mp (default:
		                False)
	  -ncpu N_CPU_MULTIPROCESS, --n_cpu_multiprocess N_CPU_MULTIPROCESS
		                number of cpu used for multiprocessing (default: 4)
	  -tl THIRD_UPPER_LEVEL, --third_upper_level THIRD_UPPER_LEVEL
		                add a directory with an upper level on top of the two previous ones,
		                the upper become intermediate (default: None)
	  -imd INTER_MATCH_DIR, --inter_match_dir INTER_MATCH_DIR
		                add a directory for the upper intermediate matching, same format as
		                the default levels matchings (default: None)
	  -imf INTER_MATCH_FILE, --inter_match_file INTER_MATCH_FILE
		                add a file for the upper intermediate matching, same format as for
		                the default levels matchings (default: None)
	  -tlh THREE_LEVEL_HEURISTIC, --three_level_heuristic THREE_LEVEL_HEURISTIC
		                heuristic for 3 level rec, can be either seq (or dec) for sequential
		                or MC for montecarlo (default: MC)
	  -tlMCs THREE_LEVEL_MC_SAMPLE, --three_level_MC_sample THREE_LEVEL_MC_SAMPLE
		                number of samples of upper intermediate reconciliation for monte
		                carlo heuristic of 3 level reconciliation (default: 10)
	  -inre INTER_N_RATES_ESTIMATION_STEPS, --inter_n_rates_estimation_steps INTER_N_RATES_ESTIMATION_STEPS
		                number of steps in the rates estimation process for the inter upper
		                rec, for each step we compute likelihood of reconciliation and set
		                rates to the observed frequency for each events (default: 5)
	  -inres INTER_N_RATES_ESTIMATION_REC_SAMPLE, --inter_n_rates_estimation_rec_sample INTER_N_RATES_ESTIMATION_REC_SAMPLE
		                in the rates estimation process of inter and upper, number of
		                scenarios sampled to estimate event frequencies (default: 100)
	  -ins INTER_N_REC_SAMPLE, --inter_n_rec_sample INTER_N_REC_SAMPLE
		                number of reconciliation scenarios sampled to compute events
		                frequency for the inter upper reconciliation (default: 100)
	  -insr INTER_N_REC_SAMPLE_RATES, --inter_n_rec_sample_rates INTER_N_REC_SAMPLE_RATES
		                number of reconciliation scenarios sampled to compute events
		                frequency for the inter upper reconciliation for rate inference
		                (default: 100)
	  -v, --verbose         increase output verbosity (default: False)
	  -dr DUPLICATION_RATE, --duplication_rate DUPLICATION_RATE
		                initial duplication rate. (default: 0.01)
	  -tr TRANSFER_RATE, --transfer_rate TRANSFER_RATE
		                initial transfer rate. (default: 0.01)
	  -lr LOSS_RATE, --loss_rate LOSS_RATE
		                initial loss rate. (default: 0.01)
	  -idr INTER_DUPLICATION_RATE, --inter_duplication_rate INTER_DUPLICATION_RATE
		                initial inter duplication rate. (default: 0.01)
	  -itr INTER_TRANSFER_RATE, --inter_transfer_rate INTER_TRANSFER_RATE
		                initial inter transfer rate. (default: 0.01)
	  -ilr INTER_LOSS_RATE, --inter_loss_rate INTER_LOSS_RATE
		                initial inter loss rate. (default: 0.01)
	  -ilo, --inter_less_output
		                for three level, output only event frequency by gene family and no
		                summation files (default: False)
	  -ia, --inter_amalgamation
		                for three level, the intermediate level can be input as a list of
		                trees as a density for amalgamation (default: False)
	  -dd, --distance_dependent
		                add dependence on distance in the tree for transfers (default: False)
	  -ir INCOMPLETE_SORTING_RATE, --incomplete_sorting_rate INCOMPLETE_SORTING_RATE
		                add the possibility of I event, some kind of incomplete sorting
		                useful in some setting, speciation but one of the child do not
		                descend. (default: 0.0)
	  -iir INTER_INCOMPLETE_SORTING_RATE, --inter_incomplete_sorting_rate INTER_INCOMPLETE_SORTING_RATE
		                add the possibility of I event for the inter reconciliation, some
		                kind of incomplete sorting useful in some setting, speciation but one
		                of the child do not descend. (default: 0.0)
	  -geo, --geo_rates     geographic null events (S, D, I) get all the same rates at each
		                inference step (default: False)
	  -slm SECOND_LEVEL_MODEL, --second_level_model SECOND_LEVEL_MODEL
		                model for the two level reconciliation, upper one in three level. Can
		                be, l for likelihood, compute likelihood (and margin ml if best),
		                joint_ml for joint maximum likelihood to get the maximum likelihood
		                scenario, and tree_ml for maximum likelihood amalgamated tree
		                (default: l)



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
