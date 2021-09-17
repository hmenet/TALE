#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:17:16 2020

@author: hmenet
"""

#######################
#sauvegarde recphyloxml
#######################

import arbre

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

def rec_to_string(tree_node,transfer_back=False):
    n_loss=0
    s=""
    s+="<clade>\n"
    innername=tree_node.name
    s+=" <name>"+innername+"</name>\n"
    s+="<eventsRec>\n"
    if transfer_back:
        s+="<transferBack destinationSpecies=\"" + e.upper +"\"/>"
    transfer_back=False
    le=tree_node.event_list
    for e in le:
        e1=e.upper
        e0=e.name
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
            s+="<loss speciesLocation=\"" + e.right +"\"/>"
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
            s+="<loss speciesLocation=\"" + e.right +"\"/>"
            s+="</eventsRec>\n"
            s+="</clade>\n"
            s+="<clade>\n"
            innername=innername+".l"
            s+=" <name>"+innername+"</name>\n"
            s+="<eventsRec>\n"
            s+="<transferBack destinationSpecies=\"" + e.left +"\"/>"
        s+="\n"
    s+="</eventsRec>\n"
    if not inner.isLeaf():
        s+=rec_to_string(tree_node.left,transfer_back=transfer_back)
        s+=rec_to_string(tree_node.right)
    s+="</clade>\n"
    for i in range(n_loss):
        s+="</clade>\n"
    return s

def save_recphyloxml(rec,rec_scenario,file):
    f=open(file, "w")
    f.write("<recPhylo>\n")
    for host_tree in rec.upper:
        f.write("<spTree>\n<phylogeny>\n")
        f.write(tree_to_string(host))
        f.write("</phylogeny>\n</spTree>\n\n")
    id_counter=1
    for reconstructed_tree in rec_scenario:
        f.write("<recGeneTree>\n<phylogeny rooted=\"true\"><id>"+str(id_counter)+"</id>")
        f.write(rec_to_string(reconstructed_tree)
        f.write("</phylogeny>\n</recGeneTree>\n")
        id_counter+=1
    f.write("</recPhylo>\n")
    f.close()

def free_living_check(reconstructed_tree):
    for u in reconstructed_tree.post_order_traversal():
        for e in u.event_list:
            if e.upper==e.lower:
                e.upper="FREE_LIVING"


#additional info not found in the recphyloxml, such as likelihood
def output_recphyloxml_additional_info(rec,upper_scenario, file):
    s="log likelihood symbiont gene:" + str(rec.lower_tree_computation.log_likelihood) + "\n"
    s+=rec.rates.pp()
    if rec.third_level:
        s="log likelihood host symbiont:" + str(upper_scenario.log_likelihood) + "\n"
        s+=rec.upper_rec.rates.pp()



def save_recphyloxml_from_rec(rec,rec_scenario,upper_scenario,file)
    for rec_tree in rec_scenario:
        free_living_check(rec_tree)
    save_recphyloxml(rec, rec_scenario,file)
    output_recphyloxml_additional_info(rec,upper_scenario,file+".additional_info")

