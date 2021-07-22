#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 16:02:46 2020

@author: hmenet
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 15/02/2021

@author: hmenet
"""
from os import walk, path, makedirs

import time

import argparse


from reconciliation import reconciliation
from rates_inference import gene_rates_ml
from arbre import save_tree
from out_recphyloxml import save_recphyloxml_from_l_event
from read_input import read_input
from event_frequency_output import output_frequency_for_all_family
from read_clade_frequencies import clade_to_name_by_fam

### arbre espèce raciné, arbre genes non raciné, possibilité de liste d'arbre pour amalgamation, possibilité plusieurs arbres de gènes (plusieurs fichiers pour plusieurs famille de gènes)

### pour l'amalgamation, je considère tout les arbres de la liste (je fais sans burn in), je ne sais pas exactement ce que fait ALE (s'il discard les premiers arbres ou non)

"""
### pour exemple ALE : matching des feuilles implicite avec feuilles des arbres d'espèces et de gènes de même nom, liste d'arbre avec amalgamation
data_dir="/home/hmenet/Documents/rewrite_ALE/ex_ALE/"
symbiont_dir="species/"
gene_dir="gene/"
file_links=[symbiont_dir, gene_dir]
for i in range(len(file_links)):
    file_links[i]=data_dir+file_links[i]
symbiont_dir, gene_dir=file_links
symbiont_list, clades_data_list, c_match_list=read_input(symbiont_dir, gene_dir)
"""


parser=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("upper_dir", type=str, help="directory with newick unrooted files for each lower level tree")
parser.add_argument("lower_dir", type=str, help="directory with newick unrooted files for each lower level tree, possibility of list of trees to use amalgamation")
parser.add_argument("-mdir", "--matching_dir", type=str, help="provide a directory with matching files for each lower level tree to match their leaves to upper tree ones. Without it, lower tree and upper tree leaves must have the same name. A lower leaf can be matched to multiple upper leaves, which will be interpreted as uncertainty on the match (not failure to diverge)")
parser.add_argument("-mfile", "--matching_file", type=str, help="provide a file with all matching infos for each lower level tree to match their leaves to upper tree ones. Alternative to -mdir")
parser.add_argument("-o", "--output", default="output/rec",type=str, help="output recphyloxml file. If multiple samples, multiple recphyloxml are generated, if best rec option, then first rec is the best one")
parser.add_argument("-of", "--output_freq", default="output/event_frequency",type=str, help="output event frequency file.")
parser.add_argument("-ns", "--n_rec_sample", type=int, default=100, help="number of reconciliation scenarios sampled to compute events frequency")
parser.add_argument("-nre", "--n_rates_estimation_steps", type=int, default=5, help="number of steps in the rates estimation process, for each step we compute likelihood of reconciliation and set rates to the observed frequency for each events")
parser.add_argument("-nres", "--n_rates_estimation_rec_sample", type=int, default=100, help="in the rates estimation process, number of scenarios sampled to estimate event frequencies")
parser.add_argument("-b", "--best_rec", action="store_true", help="return the best (maximum likelihood) reconciliation scenario.")
parser.add_argument("-nrxml", "--n_recphyloxml", type=int, default=1, help="number of sampled scenarios stored as recphyloxml file.")
parser.add_argument("-mp", "--multiprocess", default=False, help="enable multiprocessing for the sampling parts", action="store_true")
parser.add_argument("-mpf", "--multiprocess_fam", default=False, help="enable multiprocessing with one process for each lower tree family. If chosen, stop from using multiprocess for sampling -mp", action="store_true")
parser.add_argument("-tl", "--third_upper_level", default=None, help="add a directory with an upper level on top of the two previous ones, the upper become intermediate")
parser.add_argument("-imd", "--inter_match_dir", default=None, help="add a directory for the upper intermediate matching, same format as the default levels matchings")
parser.add_argument("-imf", "--inter_match_file", default=None, help="add a file for the upper intermediate matching, same format as for the default levels matchings")
parser.add_argument("-tlh", "--three_level_heuristic", default="MC", help="heuristic for 3 level rec, can be either dec for decoupled or MC for montecarlo")
parser.add_argument("-tlMCs", "--three_level_MC_sample", default=10, type=int, help="number of samples of upper intermediate reconciliation for monte carlo heuristic of 3 level reconciliation")
parser.add_argument("-inre", "--inter_n_rates_estimation_steps", type=int, default=5, help="number of steps in the rates estimation process for the inter upper rec, for each step we compute likelihood of reconciliation and set rates to the observed frequency for each events")
parser.add_argument("-inres", "--inter_n_rates_estimation_rec_sample", type=int, default=100, help="in the rates estimation process of inter and upper, number of scenarios sampled to estimate event frequencies")
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
parser.add_argument("-dr", "--duplication_rate", default=0.01, type=float, help="initial duplication rate.")
parser.add_argument("-tr", "--transfer_rate", default=0.01, type=float, help="initial transfer rate.")
parser.add_argument("-lr", "--loss_rate", default=0.01, type=float, help="initial loss rate.")

args=parser.parse_args()

if args.verbose:
    print("verbosity turned on")

if not args.third_upper_level:
    symbiont_list, clades_data_list, c_match_list, gene_file_list=read_input(args.upper_dir, args.lower_dir, leaf_matching_directory=args.matching_dir, leaf_matching_file=args.matching_file)
    upper_input=None
else:
    symbiont_list, clades_data_list, c_match_list, gene_file_list, host_list, inter_clades_data_list, inter_c_match_list, inter_file_list, inter_clade_to_tree=read_input(args.upper_dir, args.lower_dir, leaf_matching_directory=args.matching_dir, leaf_matching_file=args.matching_file, host_directory=args.third_upper_level, host_matching_file=args.inter_match_file, host_matching_dir=args.inter_match_dir)
    rates_inter=dict()
    rates_inter["T"]=0.05
    rates_inter["D"]=0.05
    rates_inter["L"]=0.05


    #test for free living

    name_to_tree=dict()
    for symbiont in symbiont_list:
        for u in symbiont.leaves():
            name_to_tree[u.name]=u
    for i_clade_inter in range(len(inter_c_match_list)):
        for c in inter_clades_data_list[i_clade_inter][2].keys():
            if len(inter_clades_data_list[i_clade_inter][2][c])==1 and not c in inter_c_match_list[i_clade_inter]:
                own_match=name_to_tree[inter_clades_data_list[i_clade_inter][2][c][0]]
                symbiont=own_match.root
                if not symbiont in host_list:
                    host_list.append(symbiont)
                inter_c_match_list[i_clade_inter][c]=[own_match]

    upper_post_order=[]
    for upper in host_list:
        upper_post_order+=upper.post_order_traversal()



    upper_input=upper_post_order, inter_clades_data_list, inter_c_match_list, rates_inter, args.three_level_heuristic, inter_clade_to_tree, args.three_level_MC_sample, host_list, args.inter_n_rates_estimation_steps, args.inter_n_rates_estimation_rec_sample

init_rates=[args.duplication_rate, args.loss_rate, args.transfer_rate]

def rec_and_output(symbiont_list, clades_data_list, c_match_list, gene_file_list, out_file="output/rec", n_sample=1, n_steps=5, n_rec_sample_rates=100, best_rec=False, n_recphyloxml=0, multiprocess=False, multiprocess_fam=False, python_output=False, out_freq_file="output/event_frequency", upper_input=None, init_rates=[0.01,0.01,0.01]):
    parasite_post_order=[]
    for symbiont in symbiont_list:
        parasite_post_order+=symbiont.post_order_traversal()

    if upper_input:
        upper_post_order, inter_clades_data_list, inter_c_match_list, rates_inter, args.three_level_heuristic, inter_clade_to_tree, args.three_level_MC_sample, host_list, inter_n_steps, inter_n_rec_sample=upper_input
        init_rates_inter=[rates_inter["D"], rates_inter["L"], rates_inter["T"]]
        t1=time.perf_counter()
        rates_inter=gene_rates_ml(upper_post_order,inter_clades_data_list,inter_c_match_list,inter_n_steps,init_rates_g=init_rates_inter, n_rec_sample=inter_n_rec_sample,multi_process=multiprocess, multi_process_family=multiprocess_fam, upper_input=None)
        cmpt_time=time.perf_counter()-t1
        print("Upper inter rates estimated in ", cmpt_time, " s")
        print("Rates", rates_inter)
        upper_input=upper_post_order, inter_clades_data_list, inter_c_match_list, rates_inter, args.three_level_heuristic, inter_clade_to_tree, args.three_level_MC_sample, host_list


    t1=time.perf_counter()
    rates=gene_rates_ml(parasite_post_order,clades_data_list,c_match_list, n_steps, init_rates_g=init_rates, n_rec_sample=n_rec_sample_rates, multi_process=multiprocess, multi_process_family=multiprocess_fam, upper_input=upper_input)
    cmpt_time=time.perf_counter()-t1
    print("Rates estimated in ", cmpt_time, " s")
    print("Rates", rates)



    t1=time.perf_counter()
    out_rec=reconciliation(parasite_post_order, clades_data_list, c_match_list, rates, sample=n_sample>0, n_sample=n_sample, best=best_rec, n_recphyloxml=n_recphyloxml, multi_process=multiprocess, multi_process_family=multiprocess_fam, upper_input=upper_input)
    if n_sample>0:
        if upper_input != None:
            likelihood, l_event_gene, l_scenarios, l_scenarios_upper, log_likelihood_list=out_rec
        else:
            likelihood,l_event_gene, l_scenarios=out_rec
    else:
        likelihood=out_rec
    cmpt_time=time.perf_counter()-t1
    print("Reconciliation ended in ", cmpt_time, " s")
    print("Log Likelihood: ", likelihood)
    print("Rates: ", rates)

    out_dir1=out_file[:out_file.rfind("/")+1]
    if not path.isdir(out_dir1):
        makedirs(out_dir1)
    out_dir2=out_freq_file[:out_freq_file.rfind("/")+1]
    if not path.isdir(out_dir2):
        makedirs(out_dir2)

    if upper_input != None:
        upper_post_order, inter_clades_data_list,inter_c_match_list, rates_inter, heuristic, inter_clade_to_tree_list, n_sample_MC, host_list=upper_input
        if heuristic=="dec":
            n_sample_MC=1
    else:
        n_sample_MC=1
    for upper_rec in range(n_sample_MC):
        if upper_input != None :
            out_file_name=out_file+str(upper_rec)+"upper"+".recphyloxml"

            inter_clade_to_name_list=clade_to_name_by_fam(inter_clades_data_list, symbiont_list)

            save_recphyloxml_from_l_event(host_list, l_scenarios_upper[upper_rec], out_file_name, c_match_list=inter_c_match_list, clade_data_list=inter_clades_data_list, clade=True, clade_to_name_list=inter_clade_to_name_list,inter_symbiont_list=symbiont_list)

        for i_recphyloxml in range(min(n_sample, n_recphyloxml)):
            scenario_by_family=l_scenarios[i_recphyloxml]
            if best_rec and i_recphyloxml==0:
                out_file_name=out_file+str(i_recphyloxml)+"_best"+".recphyloxml"
            else:
                out_file_name=out_file+str(i_recphyloxml)+".recphyloxml"
            save_recphyloxml_from_l_event(symbiont_list, scenario_by_family, out_file_name, c_match_list=c_match_list, clade_data_list=clades_data_list, clade=True)
    if n_sample>0:
        output_frequency_for_all_family(l_event_gene, gene_file_list, output_file=out_freq_file)


rec_and_output(symbiont_list, clades_data_list, c_match_list, gene_file_list, out_file=args.output, n_sample=args.n_rec_sample, n_steps=args.n_rates_estimation_steps, n_rec_sample_rates=args.n_rates_estimation_rec_sample, best_rec=args.best_rec, n_recphyloxml=args.n_recphyloxml, multiprocess=args.multiprocess, multiprocess_fam=args.multiprocess_fam, out_freq_file=args.output_freq, upper_input=upper_input, init_rates=init_rates)

