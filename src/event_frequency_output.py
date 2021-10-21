import numpy as np

def upper_match_print(upper_match):
    s=""
    s+="\t"
    for u in upper_match:
        s+="|"+u
    return s

def intersection(match1,match2):
    return len(set(match1).intersection(set(match2)))

#for one gene family
def output_event_frequency(l_event_aggregate, rec,output_file,likelihood):
    s="log likelihood:"+str(likelihood)+"\n"
    #sort by type of event
    event_sorted_by_types=dict()
    for event in l_event_aggregate:
        if rec.third_level and "T" in event.name:
            if intersection(event.upper_match,event.upper_left_match):
                event_name=event.name+"_intra"
            else:
                event_name=event.name+"_inter"
        else:
            event_name=event.name
        if not event_name in event_sorted_by_types:
            event_sorted_by_types[event_name]=[]
        event_sorted_by_types[event_name].append(event)
    for event_type in event_sorted_by_types:
        event_list=event_sorted_by_types[event_type]
        event_list.sort(key=lambda e: l_event_aggregate[e], reverse=True)
    if rec.third_level:
        event_type_list=["T_inter","T_intra","TL_inter","TL_intra","D","SL","S","I"]
    else:
        event_type_list=["T","TL","D","SL","S","I"]
    for event_type in event_type_list:
        if event_type in event_sorted_by_types:
            if event_type in ["S", "D"] :
                s+=event_type+"\tparent species\tobserved frequency\tlower node"
                if rec.third_level:
                    s+="\tupper match"
                s+="\n"
                for event in event_sorted_by_types[event_type]:
                    s+=event_type+"\t"+event.upper+"\t"+str(l_event_aggregate[event])+"\t"+event.lower
                    if rec.third_level:
                        s+=upper_match_print(event.upper_match)

                    s+="\n"
            if event_type=="SL":
                s+="SL\tparent species\tchild species keeping the gene\tchild species losing the gene\tobserved frequency\t lower node"
                if rec.third_level:
                    s+="\tupper match\tupper left match\tupper right match"
                s+="\n"
                for event in event_sorted_by_types["SL"]:
                    s+="SL\t"+event.upper+"\t"+event.upper_left_or_keeper_or_receiver+"\t"+event.upper_right_or_loser_or_donor+"\t"+str(l_event_aggregate[event])+"\t"+event.lower
                    if rec.third_level:
                        s+=upper_match_print(event.upper_match)
                        s+=upper_match_print(event.upper_left_match)
                        s+=upper_match_print(event.upper_right_match)


                    s+="\n"
            if event_type in ["T","T_inter","T_intra","I"]:
                s+=event_type+"\tgiving species\treceiving species\tobserved frequency\tlower node\tlower staying\tlower leaving"
                if rec.third_level:
                    s+="\tupper giving match\tupper receiving match"
                s+="\n"
                for event in event_sorted_by_types[event_type]:
                    s+=event_type+"\t"+event.upper+"\t"+event.upper_left_or_keeper_or_receiver+"\t"+str(l_event_aggregate[event])+"\t"+event.lower+"\t"+event.lower_right+"\t"+event.lower_left
                    if rec.third_level:
                        s+=upper_match_print(event.upper_match)
                        s+=upper_match_print(event.upper_left_match)
                    s+="\n"
            if event_type in ["TL","TL_inter","TL_intra"]:
                s+=event_type+"\tgiving species\treceiving species\tobserved frequency\tlower node"
                if rec.third_level:
                    s+="\tupper giving match\tupper receiving match"
                s+="\n"
                for event in event_sorted_by_types[event_type]:
                    s+=event_type+"\t"+event.upper+"\t"+event.upper_left_or_keeper_or_receiver+"\t"+str(l_event_aggregate[event])+"\t"+event.lower
                    if rec.third_level:
                        s+=upper_match_print(event.upper_match)
                        s+=upper_match_print(event.upper_left_match)
                    s+="\n"

    f=open(output_file, "w")
    f.write(s)
    f.close()

def output_frequency_for_all_family(l_event_by_family, rec, likelihood_by_gene,output_file="event_frequency"):
    for i in range(len(l_event_by_family)):
        output_event_frequency(l_event_by_family[i],rec, output_file+rec.lower[i].tree_name,likelihood_by_gene[i])

def save_likelihood(likelihood_list,out_file):
    s=""
    for likelihood in likelihood_list:
        s+=str(likelihood)+"\n"
    f=open(out_file, "w")
    f.write(s)
    f.close()




