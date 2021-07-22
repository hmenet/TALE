#read output

from os import walk

import matplotlib.pyplot as plt

def construct_file_list(directory):
    liste_fichiers = []
    for (repertoire, sousRepertoires, fichiers) in walk(directory):
        liste_fichiers.extend(fichiers)
    return liste_fichiers


def read_output_files(directory):
    liste_fichiers=construct_file_list(directory)
    l_events_aggregate=dict()
    l_events_by_fam=[]
    for file_name in liste_fichiers:
        l_events=dict()
        f=open(directory+file_name, "r")
        s=f.read()
        f.close()
        list_events=s.split(sep="\n")
        seen_events=dict()
        for u in list_events:
            v=u.split(sep="\t")
            event=tuple(v[:-1])
            if len(event)>0:
                if not event[0] in seen_events:
                    seen_events[event[0]]=1
                else:
                    event_freq=float(v[-1])
                    l_events_aggregate.setdefault(event,0)
                    l_events.setdefault(event,0)
                    l_events_aggregate[event]+=event_freq
                    l_events[event]+=event_freq
        l_events_by_fam.append(l_events)
    return l_events_aggregate, l_events_by_fam





path_dir="/home/hmenet/Documents/Stage_M2/These/script/h_pylori/test_010721/resultats/events_by_gene/"

l_event_aggregate,l_events_by_fam=read_output_files(path_dir)

path_file_match_hp="/home/hmenet/Documents/Stage_M2/These/script/h_pylori/test_010721/resultats/3level_info0_hp"

def read_output_match_hp(match_hp_file):
    f=open(match_hp_file)
    s=f.read()
    f.close()
    match_hp=dict()
    l1=s.split(sep="\n")
    for u in l1:
        v=u.split(sep="\t")
        if len(v)>1:
            match_hp.setdefault(v[0],[])
            match_hp[v[0]].append(v[1])
    return match_hp

match_hp=read_output_match_hp(path_file_match_hp)

def transfer_frequency_only_ext(l_event, match_hp, output_file):

    aggregate_donor=dict()
    aggregate_receiver=dict()
    for u in l_event:
        if "T" in u[0]:
            e=u[1]
            h=u[2]
            if len(set(match_hp[e]) & set(match_hp[h])) ==0: #True:
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
        l.sort(key=lambda tup: tup[0])
        for upper_host,e,freq in l:
            s+=upper_host +"\t"+e
            if i==0:
                s+="->"
            else:
                s+="<-"
            s+=str(freq)+"\n"
        i+=1

    f=open(output_file+"transfer_freq_only_ext","w")
    f.write(s)
    f.close()

#transfer_frequency_only_ext(l_event_aggregate, match_hp, "/home/hmenet/Documents/Stage_M2/These/script/h_pylori/test_010721/resultats/3level_info0_onlyext")


def transfer_direction(l_event_by_family, match_hp):
    aggregate_transfer=dict()
    aggregate_transfer_by_fam=dict()
    i_fam=-1
    for l_event in l_event_by_family:
        i_fam+=1
        aggregate_transfer_this_fam=dict()
        for u in l_event:
            if "T" in u[0]:
                e_list=match_hp[u[1]]
                h_list=match_hp[u[2]]
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
    #for e in aggregate_transfer:
    #    for h in aggregate_transfer[e]:
    #        s+=e.name
    #        s+="->"
    #        s+=h.name
    #        s+=":"
    #        s+=str(aggregate_transfer[e][h])
    #        s+="\n"
    #f=open(output_file, "w")
    #f.write(s)
    #f.close()
    l_delta_transfer=dict()
    s=""
    for e in aggregate_transfer_by_fam:
        l_delta_transfer[e]=dict()
        for h in aggregate_transfer_by_fam[e]:
            if e!=h:
                if not h in aggregate_transfer_by_fam or not e in aggregate_transfer_by_fam[h]:
                    list_to_use=[0 for i in range(len(aggregate_transfer_by_fam[e][h]))]
                else:
                    list_to_use=aggregate_transfer_by_fam[h][e]
                l=[aggregate_transfer_by_fam[e][h][fam]-list_to_use[fam] for fam in range(len(aggregate_transfer_by_fam[e][h]))]
                l_delta_transfer[e][h]=l
                l.sort()
                #n_fam=len(aggregate_transfer_by_fam[e][h])
                #s+=e.name + "(-> - <-)"+h.name+" min "+str(l[0])+" q1 "+str(l[int(n_fam/4)])+ " med "+ str(l[int(n_fam/2)])+" q2 "+str(l[int(3*n_fam/4)]) + " max "+str(l[-1])+"\n"

    #f=open(output_file+"transfer_direction", "w")
    #f.write(s)
    #f.close()

    return l_delta_transfer

l_delta_transfer=transfer_direction(l_events_by_fam, match_hp)

l=l_delta_transfer["hpAfrica1"]["hpEurope"]
plt.figure()
plt.hist(l)
plt.title("hpAfrica1(-> - <-)hpEurope dist")
plt.show()
