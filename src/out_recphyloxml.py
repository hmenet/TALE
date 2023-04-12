#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:17:16 2020

@author: hmenet
"""

#######################
#save recphyloxml
#######################

import arbre
import numpy as np


from rec_aux_func import is_mult_match

def tree_to_string(tree):
    s=""
    s+="<clade>\n"
    s+=" <name>"+tree.name+"</name>\n"
    if not tree.isLeaf():
        s+=tree_to_string(tree.left)
        s+=tree_to_string(tree.right)
    s+="</clade>\n"
    return s

def rec_to_string(tree_node,transfer_back=False,transfer_back_specie=None):
    n_loss=0
    s=""
    s+="<clade>\n"
    innername=tree_node.name
    s+=" <name>"+innername+"</name>\n"
    s+="<eventsRec>\n"
    if transfer_back:
        s+="<transferBack destinationSpecies=\"" + transfer_back_specie +"\"/>"
    transfer_back=False
    le=tree_node.event_list
    #print(len(le),[(e.name,e.upper) for e in le])


    for e in le:
        e1=e.upper
        e0=e.name
        if e0=="I":
            e0="T"



        if e0 == "S":
            s+="<speciation speciesLocation=\"" + e1 +"\"/>"
        if e0 == "D":
            s+="<duplication speciesLocation=\"" + e1 +"\"/>"
        if e0 == "T":
            s+="<branchingOut speciesLocation=\"" + e1 +"\"/>"
            transfer_back=True
        if e0 =="E":
            s+="<leaf speciesLocation=\"" + e1 +"\"/>"
        if e0 =="SL":
            s+="<speciation speciesLocation=\"" + e1 +"\"/>"
            s+="\n"
            s+="</eventsRec>\n"
            s+="<clade>\n"
            n_loss+=1
            s+=" <name>"+innername+"LOSS"+"</name>\n"
            s+="<eventsRec>\n"
            if e1 =="FREE_LIVING":
                s+="<loss speciesLocation=\"" + e1 +"\"/>"
            else:
                s+="<loss speciesLocation=\"" + e.upper_right_or_loser_or_donor +"\"/>"
            s+="</eventsRec>\n"
            s+="</clade>\n"
            s+="<clade>\n"
            innername=innername+".l"
            s+=" <name>"+innername+"</name>\n"
            s+="<eventsRec>\n"

        if e0 =="TL":
            s+="<branchingOut speciesLocation=\"" + e1 +"\"/>"
            s+="\n"
            s+="</eventsRec>\n"
            s+="<clade>\n"
            n_loss+=1
            s+=" <name>"+innername+"LOSS"+"</name>\n"
            s+="<eventsRec>\n"

            if e1 =="FREE_LIVING":
                s+="<loss speciesLocation=\"" + e1 +"\"/>"
            else:
                s+="<loss speciesLocation=\"" + e.upper_right_or_loser_or_donor +"\"/>"


            s+="</eventsRec>\n"
            s+="</clade>\n"
            s+="<clade>\n"
            innername=innername+".l"
            s+=" <name>"+innername+"</name>\n"
            s+="<eventsRec>\n"
            s+="<transferBack destinationSpecies=\"" + e.upper_left_or_keeper_or_receiver +"\"/>"
        s+="\n"
    s+="</eventsRec>\n"
    if not tree_node.isLeaf():
        s+=rec_to_string(tree_node.left,transfer_back=transfer_back, transfer_back_specie=e.upper_left_or_keeper_or_receiver)
        #avant :
        #s+=rec_to_string(tree_node.left,transfer_back=transfer_back, transfer_back_specie=e.upper)

        s+=rec_to_string(tree_node.right)
    s+="</clade>\n"
    for i in range(n_loss):
        s+="</clade>\n"
    return s

def save_recphyloxml(rec,rec_scenario_by_fam,file):
    f=open(file, "w")
    f.write("<recPhylo>\n")
    for host_tree in rec.upper.tree_list:
        if not host_tree.added_for_free_living:
            f.write("<spTree>\n<phylogeny>\n")
            f.write(tree_to_string(host_tree))
            f.write("</phylogeny>\n</spTree>\n\n")
    id_counter=1
    for rec_scenario in rec_scenario_by_fam:
        f.write("<recGeneTree>\n<phylogeny rooted=\"true\"><id>"+str(id_counter)+"</id>")
        f.write(rec_to_string(rec_scenario.reconstructed_lower))
        f.write("</phylogeny>\n</recGeneTree>\n")
        id_counter+=1
    f.write("</recPhylo>\n")
    f.close()



#additional info not found in the recphyloxml, such as likelihood
def output_recphyloxml_additional_info(rec,rec_sol, file,rec_scenario_by_fam):
    s="log likelihood symbiont gene:" + str(rec_sol.log_likelihood) + "\n"
    if not rec_sol.upper_log_likelihood is None:
        s+="log likelihood host symbiont:" + str(rec_sol.upper_log_likelihood) + "\n"
    lower_scenario_likelihood=sum([scen.log_likelihood for scen in rec_scenario_by_fam])
    s+="log likelihood symbiont gene scenario:" + str(lower_scenario_likelihood) + "\n"
    s+=rec.rates.pp()+"\n"
    #if rec.third_level:
    #    upper_scenario_likelihood=np.log(sum([np.exp(u.log_likelihood) for u in rec_sol.upper_scenario]))

    #    s="log likelihood host symbiont:" + str(rec_sol.log_likelihood) + "\n"
    #   s+="log likelihood host symbiont scenario:" + str(upper_scenario_likelihood) + "\n"
    #    s+=rec.upper_rec.rates.pp()+"\n"
    #    s+="heuristic:"+rec.heuristic

    s+="time (in s):"+str(rec.lower_tree_computation.time)

    f=open(file,"w")
    f.write(s)



def save_recphyloxml_from_rec(rec,rec_scenario_by_fam,file, rec_sol):
    save_recphyloxml(rec, rec_scenario_by_fam,file)
    output_recphyloxml_additional_info(rec,rec_sol,file+".additional_info",rec_scenario_by_fam)
    for host_tree in rec.upper.tree_list:
        arbre.save_tree(host_tree,file+"upper_tree.nwk")
    for rec_scenario in rec_scenario_by_fam:
        arbre.save_tree(rec_scenario.reconstructed_lower, file+"lower_tree.nwk")
