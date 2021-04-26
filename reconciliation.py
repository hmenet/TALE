import multiprocessing as mp



from gene_species_rec import compute_upper_gene_E, compute_upper_gene_P
from gene_sample_rec import sample_gene_upper_rec

def job_todo(entry_tuple):
                    P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list_i_clade, log_l, corr_size,best,i=entry_tuple
                    if i==0:
                        return sample_gene_upper_rec(P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list_i_clade, log_l, corr_size, best=best)
                    else:
                        return sample_gene_upper_rec(P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list_i_clade, log_l, corr_size, best=False)



def reconciliation(parasite_post_order, clades_data_list, c_match_list,rates_g, sample=True, n_sample=100, best=False, n_recphyloxml=0, multi_process=False):

    if sample:
        l_event_by_family=[dict() for i in range(len(clades_data_list))]
        l_matching_by_family=[dict() for i in range(len(clades_data_list))]
        l_scenarios=[[[] for j in range(len(clades_data_list))] for i in range(n_recphyloxml)]
    log_likelihood=0

    E,Eavg_no_log=compute_upper_gene_E(parasite_post_order, rates_g)
    for i_clade in range(len(clades_data_list)):
        clades_data=clades_data_list[i_clade]
        #pour un seul gène
        P,P_TL, log_l,corr_size=compute_upper_gene_P(parasite_post_order,E, Eavg_no_log, clades_data, rates_g, c_match_list[i_clade],more_output=True)
        log_likelihood += log_l
        if sample :
            l_events_aggregate=l_event_by_family[i_clade]

            if not multi_process:
                for i_sample in range(n_sample):
                    #only the first scenarios is sampled as the best one (if best selected), the following ones are always sampled randomly
                    if i_sample==0:
                        l_events=sample_gene_upper_rec(P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list[i_clade], log_l, corr_size, best=best)
                    else:
                        l_events=sample_gene_upper_rec(P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list[i_clade], log_l, corr_size, best=False)
                    if i_sample < n_recphyloxml:
                        l_scenarios[i_sample][i_clade]=list(l_events)
                    for event in l_events:
                        if event in l_events_aggregate:
                            l_events_aggregate[event]+=1
                        else:
                            l_events_aggregate[event]=1
            else:
                ###parallel version with multiprocessing for sampling
                #input_tuple=(P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list[i_clade], log_l, corr_size,best)
                #job_list=[i for i in range(n_sample)]
                job_list=[(P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list[i_clade], log_l, corr_size,best,i) for i in range(n_sample)]
                #def job_todo(i):
                #    if i==0:
                #        return sample_gene_upper_rec(P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list[i_clade], log_l, corr_size, best=best)
                #    else:
                #        return sample_gene_upper_rec(P, P_TL,E, parasite_post_order, clades_data, rates_g,c_match_list[i_clade], log_l, corr_size, best=False)
                if True:
                    pool=mp.Pool(mp.cpu_count())
                    l_events_list=pool.map(job_todo, job_list)
                    pool.close()
                #aggregation of results:
                for i_sample in range(n_sample):
                    l_event=l_events_list[i_sample]
                    if i_sample < n_recphyloxml:
                        l_scenarios[i_sample][i_clade]=list(l_event)
                    for event in l_event:
                        if event in l_events_aggregate:
                            l_events_aggregate[event]+=1
                        else:
                            l_events_aggregate[event]=1



    if sample:
        #effectif -> frequences
        for l_events_aggregate in l_event_by_family:
            for event in l_events_aggregate:
                    l_events_aggregate[event]/=n_sample
    if sample:
        return log_likelihood, l_event_by_family, l_scenarios
    else:
        return log_likelihood
