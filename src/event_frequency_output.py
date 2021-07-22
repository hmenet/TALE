import numpy as np


#for one gene family
def output_event_frequency(l_event_aggregate, output_file):
    s=""
    #sort by type of event
    event_sorted_by_types=dict()
    for event in l_event_aggregate:
        if not event[0] in event_sorted_by_types:
            event_sorted_by_types[event[0]]=[]
        event_sorted_by_types[event[0]].append(event)
    for event_type in event_sorted_by_types:
        event_list=event_sorted_by_types[event_type]
        event_list.sort(key=lambda e: l_event_aggregate[e], reverse=True)
    for event_type in ["T","TL","D","SL","S"]:
        if event_type in event_sorted_by_types:
            if event_type in ["S", "D"] :
                s+=event_type+"\tparent species\tobserved frequency\n"
                for event in event_sorted_by_types[event_type]:
                    s+=event_type+"\t"+event[1].name+"\t"+str(l_event_aggregate[event])+"\n"
            if event_type=="SL":
                s+="SL\tparent species\tchild species keeping the gene\tchild species losing the gene\tobserved frequency\n"
                for event in event_sorted_by_types["SL"]:
                    s+="SL\t"+event[1].name+"\t"+event[3].name+"\t"+event[5].name+"\t"+str(l_event_aggregate[event])+"\n"
            if event_type in ["T","TL"]:
                s+=event_type+"\tgiving species\treceiving species\tobserved frequency\n"
                for event in event_sorted_by_types[event_type]:
                    s+=event[0]+"\t"+event[1].name+"\t"+event[3].name+"\t"+str(l_event_aggregate[event])+"\n"
    f=open(output_file, "w")
    f.write(s)
    f.close()

def output_frequency_for_all_family(l_event_by_family, gene_family_list, output_file="event_frequency"):
    for i in range(len(l_event_by_family)):
        output_event_frequency(l_event_by_family[i], output_file+gene_family_list[i])



def compute_likelihood_upper(upper_scenar_by_fam, rates_inter, E):
    log_l=1
    log_rates=dict()
    for e in rates_inter:
        log_rates[e]=np.log(rates_inter[e])
    if not "S" in rates_inter:
        log_rates["S"]=np.log(1-sum([rates_inter[e] for e in rates_inter]))
    for l_e in upper_scenar_by_fam:
        for u in l_e:
            for event in u[0]:
                if event=="L":
                    log_l+=log_rates["L"]
                    #log_l+=E[u[5]]
                else:
                    log_l+=log_rates[event]
    return log_l

def output_3level_transfer_info(l_event_by_family, output_file, lower_log_l, match_hp, rates_inter, E, upper_scenar_by_fam):

    upper_log_l=compute_likelihood_upper(upper_scenar_by_fam, rates_inter, E)

    s="log likelihood host symbiont: " + str(upper_log_l) + "\n" + "log likelihood symbiont gene:" + str(lower_log_l) + "\n"
    aggregate_transfer=dict()
    aggregate_transfer_by_fam=dict()
    i_fam=-1
    for l_event in l_event_by_family:
        i_fam+=1
        aggregate_transfer_this_fam=dict()
        for u in l_event:
            if "T" in u[0]:
                e_list=match_hp[u[1]]
                h_list=match_hp[u[3]]
                for e in e_list:
                    for h in h_list:
                        if not e in aggregate_transfer_this_fam:
                            aggregate_transfer_this_fam[e]=dict()
                        if not h in aggregate_transfer_this_fam[e]:
                            aggregate_transfer_this_fam[e][h]=0
                        if not e in aggregate_transfer:
                            aggregate_transfer[e]=dict()
                            aggregate_transfer_by_fam[e]=dict()
                        if not h in aggregate_transfer[e]:
                            aggregate_transfer[e][h]=0
                            aggregate_transfer_by_fam[e][h]=[0 for i in range(len(l_event_by_family))]
                        if l_event[u]>0.02:#1/n_sample
                            aggregate_transfer[e][h]+=l_event[u]
                            aggregate_transfer_this_fam[e][h]+=l_event[u]
        for e in aggregate_transfer_by_fam:
            for h in aggregate_transfer_by_fam[e]:
                if e in aggregate_transfer_this_fam and h in aggregate_transfer_this_fam[e]:
                    aggregate_transfer_by_fam[e][h][i_fam]=aggregate_transfer_this_fam[e][h]
    for e in aggregate_transfer:
        for h in aggregate_transfer[e]:
            s+=e.name
            s+="->"
            s+=h.name
            s+=":"
            s+=str(aggregate_transfer[e][h])
            s+="\n"
    f=open(output_file, "w")
    f.write(s)
    f.close()



    s=""
    for e in aggregate_transfer_by_fam:
        for h in aggregate_transfer_by_fam[e]:
            if e!=h:
                if not h in aggregate_transfer_by_fam or not e in aggregate_transfer_by_fam[h]:
                    list_to_use=[0 for i in range(len(aggregate_transfer_by_fam[e][h]))]
                else:
                    list_to_use=aggregate_transfer_by_fam[h][e]
                l=[aggregate_transfer_by_fam[e][h][fam]-list_to_use[fam] for fam in range(len(aggregate_transfer_by_fam[e][h]))]
                l.sort()
                n_fam=len(aggregate_transfer_by_fam[e][h])
                s+=e.name + "(-> - <-)"+h.name+" min "+str(l[0])+" q1 "+str(l[int(n_fam/4)])+ " med "+ str(l[int(n_fam/2)])+" q2 "+str(l[int(3*n_fam/4)]) + " max "+str(l[-1])+"\n"

    f=open(output_file+"transfer_direction", "w")
    f.write(s)
    f.close()


    aggregate_donor=dict()
    aggregate_receiver=dict()
    for l_event in l_event_by_family:
        for u in l_event:
            if "T" in u[0]:
                e=u[1]
                h=u[3]
                aggregate_donor.setdefault(e,0)
                aggregate_receiver.setdefault(h,0)
                if l_event[u]>0.02:
                    aggregate_donor[e]+=l_event[u]
                    aggregate_receiver[h]+=l_event[u]

    s=""
    i=0
    for l_aggregate in [aggregate_donor, aggregate_receiver]:
        l=[]
        for e in l_aggregate:
            for upper_host in match_hp[e]:
                l.append((upper_host,e,l_aggregate[e]))
        l.sort(key=lambda tup: tup[0].name)
        for upper_host,e,freq in l:
            s+=upper_host.name +"\t"+e.name
            if i==0:
                s+="->"
            else:
                s+="<-"
            s+=str(freq)+"\n"
        i+=1

    f=open(output_file+"transfer_freq","w")
    f.write(s)
    f.close()

def output_3level_matching_info(file, match_hp):
    s=""
    for p in match_hp:
        for h in match_hp[p]:
            s+=p.name + "\t" + h.name +"\n"
    f=open(file, 'w')
    f.write(s)
    f.close()


