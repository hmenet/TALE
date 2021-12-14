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
from out_recphyloxml import save_recphyloxml_from_rec
from read_input import read_input
from event_frequency_output import output_frequency_for_all_family,save_likelihood, save_end_likelihood

from rec_classes import*

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
parser.add_argument("-mf", "--matching_file", type=str, help="provide a file with all matching infos for each lower level tree to match their leaves to upper tree ones. Alternative to -mdir")
parser.add_argument("-o", "--output", default="output/rec",type=str, help="output recphyloxml file. If multiple samples, multiple recphyloxml are generated, if best rec option, then first rec is the best one")
parser.add_argument("-ns", "--n_rec_sample", type=int, default=100, help="number of reconciliation scenarios sampled to compute events frequency")
parser.add_argument("-nre", "--n_rates_estimation_steps", type=int, default=5, help="number of steps in the rates estimation process, for each step we compute likelihood of reconciliation and set rates to the observed frequency for each events")
parser.add_argument("-nres", "--n_rates_estimation_rec_sample", type=int, default=100, help="in the rates estimation process, number of scenarios sampled to estimate event frequencies")
parser.add_argument("-b", "--best_rec", action="store_true", help="return the best (maximum likelihood) reconciliation scenario.")
parser.add_argument("-nrxml", "--n_recphyloxml", type=int, default=1, help="number of sampled scenarios stored as recphyloxml file (if multiple upper scenario sampled, number of lower scenarios for each upper one).")
parser.add_argument("-inrxml", "--inter_n_recphyloxml", type=int, default=1, help="number of sampled inter upper scenarios stored as recphyloxml file.")
parser.add_argument("-mpf", "--multiprocess_fam", default=False, help="enable multiprocessing with one process for each lower tree family. If chosen, stop from using multiprocess for sampling -mp", action="store_true")
parser.add_argument("-ncpu", "--n_cpu_multiprocess", default=4, type=int,help="number of cpu used for multiprocessing")
parser.add_argument("-tl", "--third_upper_level", default=None, help="add a directory with an upper level on top of the two previous ones, the upper become intermediate")
parser.add_argument("-imd", "--inter_match_dir", default=None, help="add a directory for the upper intermediate matching, same format as the default levels matchings")
parser.add_argument("-imf", "--inter_match_file", default=None, help="add a file for the upper intermediate matching, same format as for the default levels matchings")
parser.add_argument("-tlh", "--three_level_heuristic", default="MC", help="heuristic for 3 level rec, can be either dec for decoupled or MC for montecarlo")
parser.add_argument("-tlMCs", "--three_level_MC_sample", default=10, type=int, help="number of samples of upper intermediate reconciliation for monte carlo heuristic of 3 level reconciliation")
parser.add_argument("-inre", "--inter_n_rates_estimation_steps", type=int, default=5, help="number of steps in the rates estimation process for the inter upper rec, for each step we compute likelihood of reconciliation and set rates to the observed frequency for each events")
parser.add_argument("-inres", "--inter_n_rates_estimation_rec_sample", type=int, default=100, help="in the rates estimation process of inter and upper, number of scenarios sampled to estimate event frequencies")
parser.add_argument("-ins", "--inter_n_rec_sample", type=int, default=100, help="number of reconciliation scenarios sampled to compute events frequency for the inter upper reconciliation")
parser.add_argument("-insr", "--inter_n_rec_sample_rates", type=int, default=100, help="number of reconciliation scenarios sampled to compute events frequency for the inter upper reconciliation for rate inference")
parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
parser.add_argument("-dr", "--duplication_rate", default=0.01, type=float, help="initial duplication rate.")
parser.add_argument("-tr", "--transfer_rate", default=0.01, type=float, help="initial transfer rate.")
parser.add_argument("-lr", "--loss_rate", default=0.01, type=float, help="initial loss rate.")
parser.add_argument("-idr", "--inter_duplication_rate", default=0.01, type=float, help="initial inter duplication rate.")
parser.add_argument("-itr", "--inter_transfer_rate", default=0.01, type=float, help="initial inter transfer rate.")
parser.add_argument("-ilr", "--inter_loss_rate", default=0.01, type=float, help="initial inter loss rate.")
parser.add_argument("-ilo", "--inter_less_output", action="store_true", help="for three level, output only event frequency by gene family and no summation files")
parser.add_argument("-ia","--inter_amalgamation",action="store_true",help="for three level, the intermediate level can be input as a list of trees as a density for amalgamation")
parser.add_argument("-dd","--distance_dependent",action="store_true",help="add dependence on distance in the tree for transfers")
parser.add_argument("-ir","--incomplete_sorting_rate",default=0.0, type=float, help="add the possibility of I event, some kind of incomplete sorting useful in some setting, speciation but one of the child do not descend.")
parser.add_argument("-iir","--inter_incomplete_sorting_rate",default=0.0, type=float, help="add the possibility of I event for the inter reconciliation, some kind of incomplete sorting useful in some setting, speciation but one of the child do not descend.")
parser.add_argument("-geo","--geo_rates",action="store_true",help="geographic null events (S, D, I) get all the same rates at each inference step")
parser.add_argument("-slm", "--second_level_model", default="l", help="model for the two level reconciliation, upper one in three level. Can be, l for likelihood, compute likelihood (and margin ml if best), joint_ml for joint maximum likelihood to get the maximum likelihood scenario, and tree_ml for maximum likelihood amalgamated tree")



args=parser.parse_args()

if args.verbose:
    print("verbosity turned on")

if args.third_upper_level and args.three_level_heuristic!="unaware":
    symbiont_list, am_tree_list, host_list, inter_am_tree_list=read_input(args.upper_dir, args.lower_dir, leaf_matching_directory=args.matching_dir, leaf_matching_file=args.matching_file, host_directory=args.third_upper_level, host_matching_file=args.inter_match_file, host_matching_dir=args.inter_match_dir, inter_amalgamation=args.inter_amalgamation)


    #test for free living

    d=symbiont_list.name_to_tree()
    for inter_am_tree in inter_am_tree_list:
        for u in inter_am_tree.leaves:
            u.corresponding_tree=d[u.name]

    for inter_am_tree in inter_am_tree_list:
        for c in inter_am_tree.reverse_post_order:
            if c.is_leaf() and c.match is None:
                if args.verbose:
                    print("FREE LIVING detected")
                c.match = [c.corresponding_tree]
                symbiont=c.corresponding_tree.root
                if not symbiont in host_list:
                    host_list.append(symbiont)
                    symbiont.added_for_free_living=True


    rec_upper_problem=Rec_problem(symb_list=host_list,amal_genes=inter_am_tree_list)
    rec_upper_problem.output_path = args.output
    rec_upper_problem.n_sample = args.inter_n_rec_sample
    rec_upper_problem.n_steps = args.inter_n_rates_estimation_steps
    rec_upper_problem.n_rec_sample_rates = args.inter_n_rec_sample_rates
    rec_upper_problem.rates = Event_rates(tr=args.inter_transfer_rate,lr=args.inter_loss_rate,dr=args.inter_duplication_rate,ir=args.inter_incomplete_sorting_rate)
    rec_upper_problem.ncpu = args.n_cpu_multiprocess
    #rec_upper_problem.multiprocess_fam=args.multiprocess_fam
    rec_upper_problem.n_output_scenario=args.inter_n_recphyloxml
    rec_upper_problem.dd=args.distance_dependent
    rec_upper_problem.geo=args.geo_rates
    rec_upper_problem.slm=args.second_level_model

    rec_problem=Rec_problem(symb_list=symbiont_list,amal_genes=am_tree_list)
    rec_problem.third_level=True
    rec_problem.heuristic=args.three_level_heuristic
    if args.distance_dependent:
        rec_problem.heuristic="dd_"+args.three_level_heuristic

    rec_problem.mc_samples=args.three_level_MC_sample
    rec_problem.upper_rec=rec_upper_problem

    """
    d=0
    d2=0
    l=[]
    print(len(inter_am_tree_list))
    for am_tree in inter_am_tree_list:
        for clade in am_tree.reverse_post_order:
            d+=len(clade.child_frequencies)
            for u in clade.child_frequencies:
                l.append(clade.child_frequencies[u])
            d2+=len(clade.child_frequencies)/len(am_tree.reverse_post_order)


    import matplotlib.pyplot as plt
    #plt.figure()
    #plt.hist(l,range=(0,1))
    #plt.show()
    print("diametre amalgamation", d, d2,len(inter_am_tree_list[0].reverse_post_order))


    #inter_am_tree_list[0].pruning()

    d=0
    d2=0
    l=[]
    print(len(inter_am_tree_list))
    for am_tree in inter_am_tree_list:
        for clade in am_tree.reverse_post_order:
            d+=len(clade.child_frequencies)
            for u in clade.child_frequencies:
                l.append(clade.child_frequencies[u])
            d2+=len(clade.child_frequencies)/len(am_tree.reverse_post_order)

    import matplotlib.pyplot as plt
    plt.figure()
    plt.hist(l,range=(0,1))
    plt.show()
    print("diametre amalgamation", d, d2,len(inter_am_tree_list[0].reverse_post_order),len(inter_am_tree_list[0].leaves))
    #for u in inter_am_tree_list[0].leaves:
    #    print(u.name)

    seen=dict()
    """
    #for u in inter_am_tree_list[0].reverse_post_order:
    #    print(u.clade_leaves)
else:
    symbiont_list, am_tree_list=read_input(args.upper_dir, args.lower_dir, leaf_matching_directory=args.matching_dir, leaf_matching_file=args.matching_file)
    rec_problem=Rec_problem(symb_list=symbiont_list,amal_genes=am_tree_list)
    rec_problem.dd=args.distance_dependent
    rec_problem.geo=args.geo_rates
    rec_problem.slm=args.second_level_model





#initializing method parameters

rec_problem.output_path = args.output
rec_problem.n_sample = args.n_rec_sample
rec_problem.n_steps = args.n_rates_estimation_steps
rec_problem.n_rec_sample_rates = args.n_rates_estimation_rec_sample
rec_problem.best_rec = args.best_rec
rec_problem.rates = Event_rates(tr=args.transfer_rate,lr=args.loss_rate,dr=args.duplication_rate,ir=args.incomplete_sorting_rate)
rec_problem.ncpu = args.n_cpu_multiprocess
rec_problem.multiprocess_fam=args.multiprocess_fam
rec_problem.n_output_scenario=args.n_recphyloxml

"""

RESTE À AJOUTER

less_output=args.inter_less_output
verbose=args.verbose

"""



def rec_and_output(rec):

    if rec.third_level:
        t1=time.perf_counter()
        gene_rates_ml(rec.upper_rec)
        cmpt_time=time.perf_counter()-t1
        print("Upper inter rates estimated in ", cmpt_time, " s")
        print("Upper rates", rec.upper_rec.rates.pp())


    t1=time.perf_counter()
    rates=gene_rates_ml(rec)
    cmpt_time=time.perf_counter()-t1
    print("Rates estimated in ", cmpt_time, " s")
    print("Rates", rec.rates.pp())



    t1=time.perf_counter()
    rec_sol=reconciliation(rec)
    cmpt_time=time.perf_counter()-t1
    rec.lower_tree_computation.time=cmpt_time

    if rec.third_level:
        if rec.heuristic=="MC":
            print("Heuristic: ", rec.heuristic, rec.mc_samples)
        else:
            print("Heuristic: ", rec.heuristic)
    else:
        rec.heuristic="2l"
        print("Heuristic: 2 Level")
    print("Reconciliation ended in ", cmpt_time, " s")
    print("Log Likelihood: ", rec_sol.log_likelihood)
    print("Rates: ", rec.rates.pp())



    #recphyloxml output
    out_file=rec.output_path

    if not path.isdir(out_file):
        makedirs(out_file)


    if rec.third_level:
        n_upper_scenario=rec.upper_rec.n_output_scenario
    else:
        n_upper_scenario=1
    likelihood_list=[]
    for i_upper_scenario in range(n_upper_scenario):
        if rec.third_level:
            rec_sol_to_use=rec_sol.upper_divided_sol[i_upper_scenario]
            out_file_name=out_file+str(i_upper_scenario)+"upper"+".recphyloxml"
            save_recphyloxml_from_rec(rec.upper_rec,rec_sol.upper_divided_sol[i_upper_scenario].upper_scenario, out_file_name,rec_sol)
            lower_scenario_list=rec_sol.upper_divided_sol[i_upper_scenario].scenario_list
            likelihood_list.append(rec_sol.upper_divided_sol[i_upper_scenario].log_likelihood)
        else:
            lower_scenario_list=rec_sol.scenario_list
            rec_sol_to_use=rec_sol
        i_recphylo=0
        for lower_scenario in lower_scenario_list:
            if rec.best_rec and i_recphylo==0:
                out_file_name=out_file+str(i_upper_scenario)+"u_"+str(i_recphylo)+"_best"+".recphyloxml"
            else:
                out_file_name=out_file+str(i_upper_scenario)+"u_"+str(i_recphylo)+".recphyloxml"
            save_recphyloxml_from_rec(rec, lower_scenario, out_file_name,rec_sol_to_use)
            i_recphylo+=1
    out_file_name=out_file+"lower_log_likelihood_list"
    save_likelihood(likelihood_list,out_file_name)
    save_end_likelihood(rec_sol.log_likelihood, out_file+"lower_log_likelihood")


    #frequency per gene output
    if rec.n_sample>0:
        output_frequency_for_all_family(rec_sol.event_list_by_fam, rec, rec_sol.log_likelihood_by_gene,output_file=out_file+"freq")

################$

rec_and_output(rec_problem)

