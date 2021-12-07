import numpy as np
import multiprocessing as mp

from gene_species_rec import compute_upper_gene_E, compute_upper_gene_P
from gene_sample_rec import sample_gene_upper_rec
from transfer_prob_pre_process import prob_transfer_sequential, distance_dependent
from rec_aux_func import log_add, log_minus, log_add_list, is_mult_match

from rec_classes import Tree_list, Rec_sol

#aggregate scenarios and create new events common to all genes
def aggregate_scenarios(scenario_list):
    d=dict()
    seen_events=dict()
    for scenario in scenario_list:
        for event in scenario.event_list:

            ekey=event.key()
            if not ekey in seen_events:
                seen_events[ekey]=event
            else:
                event=seen_events[ekey]


            d.setdefault(event,0)
            d[event]+=1

    for event in d:
        d[event]/=len(scenario_list)
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
    if rec.rate_inference:
        n_sample=rec.n_sample_rates
    else:
        n_sample=rec.n_sample
    scenario_list=[]
    if n_sample > 0 :
        for i in range(n_sample):
            if i==0:
                best_here=rec.best_rec
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
    compute_upper_gene_E(rec)

    #map output is an iterable map object
    if rec.multiprocess_fam:
        with mp.Pool(rec.ncpu) as pool:
            output_list=pool.map(family_job, rec)
    else:
        output_list=map(family_job, rec)


    l_event_by_family=[]
    log_likelihood=0
    log_likelihood_list=[]
    l_scenarios=[[] for i in range(rec.n_output_scenario)]



    #aggregation of results:
    for scenario_list, aggregate_dict, single_rec in output_list:
        l_event_by_family.append(aggregate_dict)
        log_likelihood+=single_rec.lower_tree_computation.log_l
        log_likelihood_list.append(single_rec.lower_tree_computation.log_l)

        for i in range(rec.n_output_scenario):
            l_scenarios[i].append(scenario_list[i])


    rec_sol=Rec_sol(log_likelihood,l_event_by_family,l_scenarios)
    rec_sol.log_likelihood_by_gene=log_likelihood_list


    return rec_sol

#rec is a rec_problem instance, and rate inference is True if called during rate inference
def reconciliation(rec):
    if rec.third_level:
        if rec.heuristic in ["dec","dec_no_ghost","dd_dec"]:
            rec.upper_rec.best=True
            rec.upper_rec.n_sample=1
            n_sample_MC=1
        if rec.heuristic=="MC":
            best=False
            rec.upper_rec.n_sample=rec.mc_samples
            rec.upper_rec.n_output_scenario=rec.mc_samples
            n_sample_MC=rec.mc_samples
        #rec.upper_rec.n_output_scenario=rec.n_sample
        upper_rec_sol= two_level_rec(rec.upper_rec)
        compute_upper_gene_E(rec.upper_rec)

        l_scenario=upper_rec_sol.scenario_list

        ######### done till there

        rec_sol_global=Rec_sol(None,None,None)
        rec_sol_global.upper_divided_sol=[]
        rec_sol_global.upper_log_likelihood=upper_rec_sol.log_likelihood

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

            tree_from_name=reconstructed_inter_list.name_to_tree()
            for am_tree in rec.lower:
                for u in am_tree.leaves:
                    u.match=[tree_from_name[v.name] for v in u.constant_match]

            rec_sol_this_upper=two_level_rec(rec)
            rec_sol_this_upper.upper_scenario=by_fam_scenario
            rec_sol_global.upper_divided_sol.append(rec_sol_this_upper)
        log_likelihood_list=[]
        l_event_aggregate_global=[dict() for u in rec.lower]
        rec_sol_global.event_list_by_fam=l_event_aggregate_global


        seen_events=dict()
        log_likelihood_by_gene_list=[0 for i in rec.lower]
        for rec_sol in rec_sol_global.upper_divided_sol:
            l_event_by_fam=rec_sol.event_list_by_fam
            log_likelihood_list.append(rec_sol.log_likelihood)
            for i_gene in range(len(rec_sol.log_likelihood_by_gene)):
                log_likelihood_by_gene_list[i_gene]+=rec_sol.log_likelihood_by_gene[i_gene]

            for (l_event,l_event_global) in zip(l_event_by_fam,l_event_aggregate_global):
                for event in l_event:
                    ekey=event.key()
                    if not ekey in seen_events:
                        seen_events[ekey]=event
                        new_event=event
                    else:
                        new_event=seen_events[ekey]


                    l_event_global.setdefault(new_event,0)
                    l_event_global[new_event]+=l_event[event]




        for l_event in l_event_aggregate_global:
            for event in l_event:
                l_event[event]/=n_sample_MC
        log_likelihood=log_add_list(log_likelihood_list) - np.log(n_sample_MC)
        rec_sol_global.log_likelihood=log_likelihood
        rec_sol_global.log_likelihood_by_upper=log_likelihood_list

        for i_gene in range(len(log_likelihood_by_gene_list)):
            log_likelihood_by_gene_list[i_gene]/=n_sample_MC
        rec_sol_global.log_likelihood_by_gene=log_likelihood_by_gene_list

        return rec_sol_global
    else:
        if rec.dd:
            P_transfer=distance_dependent(rec)
            rec.upper_tree_computation.P_transfer=P_transfer
        return two_level_rec(rec)
