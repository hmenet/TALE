import multiprocessing as mp
from gene_species_rec import compute_upper_gene_E, compute_upper_gene_P
from gene_sample_rec import sample_gene_upper_rec


#parallel job, if parallelizing on samples
def job_todo(entry_tuple):
    P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list_i_clade, log_l, corr_size,best,i=entry_tuple
    if i==0:
        best_here=best
    else:
        best_here=False
    return sample_gene_upper_rec(P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list_i_clade, log_l, corr_size, best=best_here)

#if parallelizing on gene family
def family_job(entry_tuple):
    parasite_post_order,E, Eavg_no_log, clades_data, rates_g, c_match_list_i_clade, best, multi_process, multi_process_family, sample, n_sample= entry_tuple
    P,P_TL, log_l,corr_size=compute_upper_gene_P(parasite_post_order,E, Eavg_no_log, clades_data, rates_g, c_match_list_i_clade,more_output=True)
    if sample :
        if multi_process_family or not multi_process:
            #l_events_list=map(job_todo, job_list)
            l_events_list=[]
            for i in range(n_sample):
                if i==0:
                    best_here=best
                else:
                    best_here=False
                l_events_list.append(list(sample_gene_upper_rec(P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list_i_clade, log_l, corr_size, best=best_here)))
        else:
            job_list=[(P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list_i_clade, log_l, corr_size,best,i) for i in range(n_sample)]
            #parallel version only if not already parrallel by family
            with mp.Pool(mp.cpu_count()) as pool:
                l_events_list=pool.map(job_todo, job_list)
    if sample:
        return log_l, l_events_list
    else:
        return log_l


def reconciliation(parasite_post_order, clades_data_list, c_match_list,rates_g, sample=True, n_sample=100, best=False, n_recphyloxml=0, multi_process=False, multi_process_family=False):

    if sample:
        l_event_by_family=[dict() for i in range(len(clades_data_list))]
        l_matching_by_family=[dict() for i in range(len(clades_data_list))]
        l_scenarios=[[[] for j in range(len(clades_data_list))] for i in range(n_recphyloxml)]
    log_likelihood=0

    E,Eavg_no_log=compute_upper_gene_E(parasite_post_order, rates_g)

    job_fam_list=[(parasite_post_order,E, Eavg_no_log, clades_data_list[i_clade], rates_g, c_match_list[i_clade], best, multi_process, multi_process_family, sample, n_sample) for i_clade in range(len(clades_data_list))]
    #map output is an iterable map object
    if not multi_process_family:
        fam_output=map(family_job, job_fam_list)
    else:
        with mp.Pool(mp.cpu_count()) as pool:
            fam_output=pool.map(family_job, job_fam_list)

    if sample:
        i_clade=0
        #aggregation of results:
        for (log_l, l_events_list) in fam_output:
            l_events_aggregate=l_event_by_family[i_clade]
            log_likelihood+=log_l
            i_sample=0
            for l_event in l_events_list:
                if i_sample < n_recphyloxml:
                    l_scenarios[i_sample][i_clade]=list(l_event)
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

    print(len(l_event_by_family), len(l_event_by_family[0]), len(l_event_by_family[1]))


    if sample:
        if multi_process or multi_process_family:
            new_l_event_by_family=[]
            #we go through all events to assign the trees in the input list in the events instead of the copies created for multiprocessing
            name_to_tree=dict()
            for u in parasite_post_order:
                name_to_tree[u.name]=u
            for l_events_aggregate in l_event_by_family:
                new_l_events_aggregate=dict()
                for e in l_events_aggregate:
                    new_event=(e[0],name_to_tree[e[1].name],e[2], name_to_tree[e[3].name], e[5], name_to_tree[e[5].name], e[6])
                    new_l_events_aggregate.setdefault(new_event, 0)
                    new_l_events_aggregate[new_event]+=l_events_aggregate[e]
                new_l_event_by_family.append(new_l_events_aggregate)
            l_event_by_family=new_l_event_by_family

        #effectif -> frequences
        for l_events_aggregate in l_event_by_family:
            for event in l_events_aggregate:
                    l_events_aggregate[event]/=n_sample
    print(len(l_event_by_family), len(l_event_by_family[0]), len(l_event_by_family[1]))
    if sample:
        return log_likelihood, l_event_by_family, l_scenarios
    else:
        return log_likelihood
