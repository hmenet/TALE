import numpy as np
import multiprocessing as mp

from gene_species_rec import compute_upper_gene_E, compute_upper_gene_P
from gene_sample_rec import sample_gene_upper_rec
from transfer_prob_pre_process import prob_transfer_sequential
from rec_aux_func import log_add, log_minus, log_add_list, is_mult_match

from event_frequency_output import output_3level_transfer_info, output_3level_matching_info

from rec_classes import Tree_list

def aggregate_scenarios(scenario_list):
    d=dict()
    for scenario in scenario_list:
        for event in scenario.event_list:
            ekey=event.key()
            d.setdefault(ekey,0)
            d[ekey]+=1
    for ekey in d:
        d[ekey]/=len(scenario_list)
    return d


#if parallelizing on gene family
def family_job(rec):

    ### P computation
    P,P_TL, log_l,corr_size=compute_upper_gene_P(rec)
    rec.lower_tree_computation.P=P
    rec.lower_tree_computation.P_TL=P_TL
    rec.lower_tree_computation.log_l=log_l
    rec.lower_tree_computation.corr_size=corr_size

    ### Sampling
    if rec.rates_inference:
        n_sample=rec.n_sample_rates
    else:
        n_sample=rec.n_sample
    scenario_list=[]
    if n_sample > 0 :
        for i in range(n_sample):
            if i==0:
                best_here=rec.best
            else:
                best_here=False
            sampled_scenario=sample_gene_upper_rec(rec,best=best_here)
            scenario_list.append(sampled_scenario)
        aggregate_dict=aggregate_scenarios(scenario_list)
    return scenario_list, aggregate_dict,rec



def two_level_rec(rec):
    n_lower_tree=len(rec.lower)
    n_scenario=rec.n_output_scenario
    if n_scenario>0 :
        r_list_by_family=[]
    if rec.n_sample>0:
        l_event_by_family=[dict() for i in range(n_lower_tree)]
        l_matching_by_family=[dict() for i in range(n_lower_tree)]
        l_scenarios=[[[] for j in range(n_lower_tree)] for i in range(n_scenario)]
    log_likelihood=0
    log_likelihood_list=[]
    compute_upper_gene_E(rec)

    #map output is an iterable map object
    if multi_process_family:
        with mp.Pool(rec.ncpu) as pool:
            output_list=pool.map(family_job, rec)
    else:
        output_list=map(family_job, rec)

    i_clade=0


    #aggregation of results:
    for scenario_list,single_rec in output_list:
        l_events_aggregate=l_event_by_family[i_clade]
        log_likelihood+=single_rec.lower_tree_computation.log_l
        log_likelihood_list.append(single_rec.lower_tree_computation.log_l)
        i_clade+=1


    rec_sol=Rec_sol()
    rec_sol.log_likelihood=log_likelihood
    rec_sol.log_likelihood_by_gene=log_likelihood_list
    rec_sol.l_event_by_family=l_event_by_family
    rec_sol.l_scenario=l_scenarios


    return rec_sol

#rec is a rec_problem instance, and rate inference is True if called during rate inference
def reconciliation(rec):

    if rec.third_level:
        if rec.heuristic in ["dec","dec_no_ghost"]:
            rec.upper_rec.best=True
            rec.upper_rec.n_sample=1
            n_sample_MC=1
        if rec.heuristic=="MC":
            best=False
            rec.upper_rec.n_sample=rec.n_mc_samples
        rec.upper.n_output_scenario=rec.n_sample
        upper_rec_sol= two_level_rec(rec.upper_rec)
        compute_upper_gene_E(rec.upper_rec)

        l_scenario=upper_rec_sol.l_scenario

        ######### done till there

        rec_sol_global=Rec_sol()
        rec_sol_global.upper_divided_sol=[]

        for am_tree in rec.lower:
            am_tree.save_match()

        i_upper_sample=-1

        for by_fam_scenario in l_scenario:
            reconstructed_inter_list=[]
            i_upper_sample+=1

            for scenario in by_fam_scenario:
                reconstructed_inter=scenario.reconstructed_lower
                reconstructed_inter_list.append(reconstructed_inter)
            reconstructed_inter_list=Tree_list(reconstructed_inter_list)

            P_transfer=prob_transfer_sequential(rec.upper_rec, reconstructed_inter_list)
            rec.upper_tree_computation.P_transfer=P_transfer
            rec.upper=reconstructed_inter_list
            #we match lower to reconstructed
            for am_tree in rec.lower:
                for u in am_tree.leaves:
                    for scenario in by_fam_scenario:
                        if u.constant_match in scenario.am_tree_to_reconstructed:
                            u.match=scenario.am_tree_to_reconstructed_tree(u.constant_match)
            rec_sol_this_upper=two_level_rec(rec)
            rec_sol_this_upper.upper_scenario=by_fam_scenario
            rec_sol_global.upper_divided_sol.append(rec_sol_this_upper)
        log_likelihood_list=[]
        l_event_aggregate_global=dict()
        for rec_sol in rec_sol_global.upper_divided_sol:
            l_event_aggregate=rec_sol.l_event_aggregate
            log_likelihood.append(rec_sol.log_likelihood)
            for lower_family_l_event:
                for event in l_event_aggregate:
                    l_event_aggregate[event]/=n_sample_MC
        log_likelihood=log_add_list(log_likelihood_list) - np.log(n_sample_MC)
        return
    else:
        return two_level_rec(rec)
