
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
