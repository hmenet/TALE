import numpy as np
import multiprocessing as mp

from gene_species_rec import compute_upper_gene_E, compute_upper_gene_P
from gene_sample_rec import sample_gene_upper_rec
from transfer_prob_pre_process import prob_transfer_sequential
from rec_aux_func import log_add, log_minus, log_add_list, is_mult_match

#parallel job, if parallelizing on samples
def job_todo(entry_tuple):
    P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list_i_clade, log_l, corr_size,best,i, return_r, P_transfer=entry_tuple
    if i==0:
        best_here=best
    else:
        best_here=False
    return sample_gene_upper_rec(P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list_i_clade, log_l, corr_size, best=best_here, return_r=return_r, P_transfer=P_transfer)

#if parallelizing on gene family
def family_job(entry_tuple):
    parasite_post_order,E, Eavg_no_log, clades_data, rates_g, c_match_list_i_clade, best, multi_process, multi_process_family, sample, n_sample, return_r, P_transfer = entry_tuple
    P,P_TL, log_l,corr_size=compute_upper_gene_P(parasite_post_order,E, Eavg_no_log, clades_data, rates_g, c_match_list_i_clade,more_output=True, P_transfer=P_transfer)
    if sample :
        l_events_list=[]
        list_r=[]
        if multi_process_family or not multi_process:
            #l_events_list=map(job_todo, job_list)
            for i in range(n_sample):
                if i==0:
                    best_here=best
                else:
                    best_here=False
                l_events_tmp, r_tmp=sample_gene_upper_rec(P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list_i_clade, log_l, corr_size, best=best_here, return_r=return_r, P_transfer=P_transfer)
                l_events_list.append(list(l_events_tmp))
                list_r.append(r_tmp)
        else:
            job_list=[(P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list_i_clade, log_l, corr_size,best,i, return_r, P_transfer) for i in range(n_sample)]
            #parallel version only if not already parrallel by family
            with mp.Pool(mp.cpu_count()) as pool:
                output_tmp=pool.map(job_todo, job_list)
                for le, r in output_tmp:
                    l_events_list.append(le)
                    list_r.append(r)
    if sample:
        if return_r:
            return log_l, l_events_list, list_r
        else:
            return log_l, l_events_list
    else:
        return log_l


def two_level_rec(parasite_post_order, clades_data_list, c_match_list,rates_g, sample=True, n_sample=100, best=False, n_recphyloxml=0, multi_process=False, multi_process_family=False, return_r=False, P_transfer=None):
    if return_r:
        r_list_by_family=[]
    if sample:
        l_event_by_family=[dict() for i in range(len(clades_data_list))]
        l_matching_by_family=[dict() for i in range(len(clades_data_list))]
        l_scenarios=[[[] for j in range(len(clades_data_list))] for i in range(n_recphyloxml)]
    log_likelihood=0

    E,Eavg_no_log=compute_upper_gene_E(parasite_post_order, rates_g)

    job_fam_list=[(parasite_post_order,E, Eavg_no_log, clades_data_list[i_clade], rates_g, c_match_list[i_clade], best, multi_process, multi_process_family, sample, n_sample, return_r, P_transfer) for i_clade in range(len(clades_data_list))]
    #map output is an iterable map object
    if not multi_process_family:
        fam_output=map(family_job, job_fam_list)
    else:
        with mp.Pool(mp.cpu_count()) as pool:
            fam_output=pool.map(family_job, job_fam_list)

    if sample:
        i_clade=0
        #aggregation of results:
        for one_fam_output in fam_output:
            if return_r:
                log_l, l_events_list, list_r=one_fam_output
                r_list_by_family.append(list_r)
            else:
                log_l, l_events_list=one_fam_output
            l_events_aggregate=l_event_by_family[i_clade]
            log_likelihood+=log_l
            i_sample=0
            for l_event in l_events_list:
                if i_sample < n_recphyloxml:
                    l_scenarios[i_sample][i_clade]=list(l_event)
                    #warning : scenario is not renamed after parallelization, so the trees referenced in it may be copies
                for event in l_event:
                    if event in l_events_aggregate:
                        l_events_aggregate[event]+=1
                    else:
                        l_events_aggregate[event]=1
                i_sample+=1
            i_clade+=1
    else:
        for log_l in fam_output:
            log_likelihood+=log_l

    if sample:
        name_to_tree=dict()
        for u in parasite_post_order:
            name_to_tree[u.name]=u
        if multi_process or multi_process_family:
            new_l_event_by_family=[]
            #we go through all events to assign the trees in the input list in the events instead of the copies created for multiprocessing
            for l_events_aggregate in l_event_by_family:
                new_l_events_aggregate=dict()
                for e in l_events_aggregate:
                    new_event=(e[0],name_to_tree[e[1].name],e[2], name_to_tree[e[3].name], e[4], name_to_tree[e[5].name], e[6])
                    new_l_events_aggregate.setdefault(new_event, 0)
                    new_l_events_aggregate[new_event]+=l_events_aggregate[e]
                new_l_event_by_family.append(new_l_events_aggregate)
            l_event_by_family=new_l_event_by_family
        if return_r:
            new_r_list_by_sample=[[dict() for i in range(len(r_list_by_family))] for i in range(n_sample)]
            for i_clade in range(len(r_list_by_family)):
                r_list=r_list_by_family[i_clade]
                for i_sample in range(n_sample):
                    r=r_list[i_sample]
                    new_r=new_r_list_by_sample[i_sample][i_clade]
                    for c in r:
                        new_r[c]=[]
                        for e in r[c]:
                            new_r[c].append(name_to_tree[e.name])

        #effectif -> frequences
        for l_events_aggregate in l_event_by_family:
            for event in l_events_aggregate:
                    l_events_aggregate[event]/=n_sample
    if sample:
        if return_r:
            return log_likelihood, l_event_by_family, l_scenarios, new_r_list_by_sample
        else:
            return log_likelihood, l_event_by_family, l_scenarios
    else:
        return log_likelihood


def reconciliation(parasite_post_order, clades_data_list, c_match_list,rates_g, upper_input=None, sample=True, n_sample=100, best=False, n_recphyloxml=0, multi_process=False, multi_process_family=False):

    if upper_input:
        upper_post_order, inter_clades_data_list,inter_c_match_list, rates_inter, heuristic, inter_clade_to_tree_list, n_sample_MC=upper_input

        if heuristic=="dec":
            best=True
            n_sample_inter=1
            n_sample_MC=1
        if heuristic=="MC":
            best=False
            n_sample_inter=n_sample_MC

        log_likelihood, l_event_by_family, l_scenarios_upper, r_list_by_sample= two_level_rec(upper_post_order, inter_clades_data_list, inter_c_match_list,rates_inter, sample=sample, n_sample=n_sample_inter, best=best, n_recphyloxml=1, multi_process=multi_process, multi_process_family=multi_process_family, return_r=True)
        E,Eavg_no_log=compute_upper_gene_E(upper_post_order, rates_inter)
        E_no_log=dict()
        for h in E:
            E_no_log[h]=np.exp(E[h])

        log_likelihood_list=[]
        l_scenarios=[]
        l_event_by_family_aggregate=[dict() for i_clade in range(len(clades_data_list))]
        for r_by_fam in r_list_by_sample:
            match_hp=dict()
            for i_clade in range(len(r_by_fam)):
                r=r_by_fam[i_clade]
                clade_to_tree=inter_clade_to_tree_list[i_clade]
                clade_to_tree_rev=dict()
                for tree in clade_to_tree:
                    clade_to_tree_rev[clade_to_tree[tree]]=tree
                for c in r:
                    match_hp[clade_to_tree_rev[c]]=r[c]
            host_info=upper_post_order, rates_inter,match_hp, E_no_log

            P_transfer=prob_transfer_sequential(host_info, parasite_post_order)
            log_l, l_event_by_family, l_scenar = two_level_rec(parasite_post_order, clades_data_list, c_match_list,rates_g, sample=sample, n_sample=n_sample, best=best, n_recphyloxml=n_recphyloxml, multi_process=multi_process, multi_process_family=multi_process_family, P_transfer=P_transfer)
            log_likelihood_list.append(log_l)
            l_scenarios+=l_scenar
            for i_clade in range(len(l_event_by_family)):
                l_event=l_event_by_family[i_clade]
                l_event_aggregate=l_event_by_family_aggregate[i_clade]
                for event in l_event:
                    l_event_aggregate.setdefault(event,0)
                    l_event_aggregate[event]+=l_event[event]
        for i_clade in range(len(l_event_by_family_aggregate)):
            l_event_aggregate=l_event_by_family_aggregate[i_clade]
            for event in l_event:
                l_event[event]/=n_sample_MC
        log_likelihood=log_add_list(log_likelihood_list) - np.log(n_sample_MC)
        return log_likelihood, l_event_by_family_aggregate, l_scenarios
    else:
        return two_level_rec(parasite_post_order, clades_data_list, c_match_list,rates_g, sample=sample, n_sample=n_sample, best=best, n_recphyloxml=n_recphyloxml, multi_process=multi_process, multi_process_family=multi_process_family)
