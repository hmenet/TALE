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




def add_event(event_name,e,u,f,v,g,w, d_events, name_to_tree, t, Loss=False, Transfer=False):
    t.birth(1)
    t.left.name=v
    t.right.name=w
    name_to_tree[v]=t.left
    name_to_tree[w]=t.right
    d_events[v]=[]
    d_events[w]=[]
    if Loss :
        d_events[w]=[["L",g]]
        name_to_tree[u]=name_to_tree[v]
        rename=u.split(sep=".")[0]
        name_to_tree[rename]=name_to_tree[v]
    if Transfer:
        d_events[v]=[["transferBack", f]]
    d_events[u].append([event_name,e])

def from_l_e_to_recphylo(l_e, leaf_matching, symbiont=None, clade=False,clade_elements=None, starting_dictionary=dict()):
    t_root=arbre.Tree()
    if clade:
        t_root.name="0"
    else:
        t_root.name=symbiont.name
    t_root.id=1
    d_events=starting_dictionary#associe aux gènes leurs evenements #anciennement d
    d_events[t_root.name]=[]
    name_to_tree=dict()#associe nom et arbre #anciennement a
    if clade:
        name_to_tree['0']=t_root
    else:
        name_to_tree[t_root.name]=t_root
    for event in l_e:
        Loss=False
        Transfer=False
        #selon si on travaille avec clade ou arbre, pas même noms
        event_name,e,u,f,v,g,w=event
        if clade:
            f=f.name
            g=g.name
            e=e.name
            if len(clade_elements[u])==1:
                u=str(clade_elements[u][0])
                u=name_to_tree[u].name
            else:
                #u=str(u)
                u=name_to_tree[str(u)].name
            if len(clade_elements[v])==1:
                v=str(clade_elements[v][0])
            else:
                v=str(v)
            if len(clade_elements[w])==1:
                w=str(clade_elements[w][0])
            else:
                w=str(w)

        else:
            u=name_to_tree[u.name].name

            v=v.name
            w=w.name
            f=f.name
            g=g.name
            e=e.name

        if "L" in event_name:
            Loss = True
            w=u+"LOSS"
            v=u+".l"
            if event_name=="SL":
                event_name="S"
        if "T" in event_name:
            Transfer=True
            event_name="branchingOut"
        add_event(event_name,e,u,f,v,g,w, d_events, name_to_tree, name_to_tree[u], Loss=Loss, Transfer=Transfer)
    for leaf in t_root.leaves():
            if not "LOSS" in leaf.name: #les loss sont des feuilles ajoutés a posteriori
                rename=""
                k=0
                continue_condition=True
                while k < len(leaf.name) and continue_condition:
                    c=leaf.name[k]
                    if c!="." or leaf.name[k+1]!="l":
                        rename+=c
                        if not clade:
                            t_root.name=symbiont.name

                    else:
                        continue_condition=False
                    k+=1
                d_events[leaf.name]
                d_events[leaf.name].append(["leaf",leaf_matching[rename]])
    return d_events,t_root

def tree_to_string(tree):
    s=""
    s+="<clade>\n"
    s+=" <name>"+tree.name+"</name>\n"
    if not tree.isLeaf():
        s+=tree_to_string(tree.left)
        s+=tree_to_string(tree.right)
    s+="</clade>\n"
    return s

def rec_to_string(inner, l_events):
    s=""
    s+="<clade>\n"


    s+=" <name>"+inner.name+"</name>\n"
    s+="<eventsRec>\n"
    le=l_events[inner.name]
    for e in le:
        if e[0] == "S":
            s+="<speciation speciesLocation=\"" + e[1] +"\"/>"
        if e[0]=="L":
            s+="<loss speciesLocation=\"" + e[1] +"\"/>"
        if e[0] == "D":
            s+="<duplication speciesLocation=\"" + e[1] +"\"/>"
        if e[0] == "transferBack":
            s+="<transferBack destinationSpecies=\"" + e[1] +"\"/>"
        if e[0] =="branchingOut":
            s+="<branchingOut speciesLocation=\"" + e[1] +"\"/>"
        if e[0] =="leaf":
            s+="<leaf speciesLocation=\"" + e[1] +"\"/>"
        s+="\n"
    s+="</eventsRec>\n"
    if not inner.isLeaf():
        s+=rec_to_string(inner.left,l_events)
        s+=rec_to_string(inner.right,l_events)
    s+="</clade>\n"
    return s

def save_recphyloxml(host, gene_list, d_events_list, file):
    f=open(file, "w")
    f.write("<recPhylo>\n<spTree>\n<phylogeny>\n")
    f.write(tree_to_string(host))
    f.write("</phylogeny>\n</spTree>\n\n")
    id_counter=1
    for i_gene in range(len(gene_list)):
        gene=gene_list[i_gene]
        f.write("<recGeneTree>\n<phylogeny rooted=\"true\"><id>"+str(id_counter)+"</id>")


        f.write(rec_to_string(gene, d_events_list[i_gene]))
        f.write("</phylogeny>\n</recGeneTree>\n")
        id_counter+=1
    f.write("</recPhylo>\n")
    f.close()

def big_host(host_list):
    t_root=arbre.Tree()
    t=t_root
    t_root.name="big_host"
    for i_host in range(len(host_list)-1):
        host=host_list[i_host]
        if len(host_list)-i_host>2:
            t.birth(1)
            t.left=host
            t=t.right
            t.name="big_host"+str(i_host)
        else:
            t.left=host
            t.right=host_list[i_host+1]
    return t_root



def save_recphyloxml_from_l_event(host_list, l_event_by_family, file,symbiont_list=None, c_match_list=None, clade_data_list=None, clade=False, leaf_matching_list=None):
    if len(host_list)>1:
        host=big_host(host_list)
    else:
        host=host_list[0]
    gene_list=[]
    d_list=[]
    for i_family in range(len(l_event_by_family)):
        l_event=l_event_by_family[i_family]
        if leaf_matching_list == None:
            leaf_matching=dict()
            if clade:
                clade_elements=clade_data_list[i_family][2]
                c_match=c_match_list[i_family]
                for c in c_match:
                    c_name=str(clade_elements[c][0])
                    leaf_matching[c_name]=c_match[c].name
            else:
                symbiont=symbiont_list[i_family]
                for leaf in symbiont.leaves():
                    leaf_matching[leaf.name]=leaf.match.name
        else:
            leaf_matching=leaf_matching_list[i_family]
            if clade:
                clade_elements=clade_data_list[i_family][2]
            else:
                symbiont=symbiont_list[i_family]

        if clade:
            d,t=from_l_e_to_recphylo(l_event, leaf_matching,symbiont=None, clade_elements=clade_elements, clade=True)
        else:
            d,t=from_l_e_to_recphylo(l_event, leaf_matching,symbiont=symbiont, clade=clade)
        d_list.append(dict(d))
        gene_list.append(t)
    save_recphyloxml(host, gene_list, d_list, file)