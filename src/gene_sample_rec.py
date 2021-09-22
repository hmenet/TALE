#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:18:17 2020

@author: hmenet
"""

import random as rd
import numpy as np

from rec_aux_func import log_add, log_minus, log_add_list, is_mult_match
from arbre import Tree
from rec_classes import Rec_event, Rec_scenario

def sample_gene_upper_rec(rec, best=False):
    am_tree=rec.single_lower
    mult_gene_match=is_mult_match(am_tree)
    r=dict()



    for u in am_tree.reverse_post_order:
        r[u]=[]

    #we create a scenario instance
    in_sampling_scenario=Rec_scenario()
    reconstructed_lower=in_sampling_scenario.reconstructed_lower
    am_tree_to_reconstructed=in_sampling_scenario.am_tree_to_reconstructed_tree
    am_tree_to_reconstructed[am_tree]=reconstructed_lower
    reconstructed_lower.corresponding_am_tree=am_tree
    reconstructed_lower.name=am_tree.name

    l_sampled_events=in_sampling_scenario.event_list#list of all the events that compose the scenario we will have sampled at the end of this function

    gene=am_tree.reverse_post_order[0] #root des clades
    parasite_post_order=rec.upper.post_order

    d_r= rec.rates.ldr
    l_r= rec.rates.llr
    t_r= rec.rates.ltr
    s_r= rec.rates.lsr


    s=rec.lower_tree_computation.log_l
    P=rec.lower_tree_computation.P
    P_TL=rec.lower_tree_computation.P_TL
    E=rec.upper_tree_computation.E
    P_transfer=rec.upper_tree_computation.P_transfer
    ancestrale_correction_size=rec.lower_tree_computation.corr_size

    k=0
    x=np.log(rd.random())+s


    x_current=P[parasite_post_order[0]][gene]
    while x_current<x:
        k+=1
        x_current=log_add(x_current, P[parasite_post_order[k]][gene])
    #il faudrait des proba d'origination
    r[gene].append(parasite_post_order[k])


    clade_to_look=[gene]
    while len(clade_to_look)>0:
        c=clade_to_look.pop()

        reconstructed_tree_node=am_tree_to_reconstructed[c]

        #for c in clade_post_order:
        TL_done=False #on autorise un seul TL
        if len(r[c])>0:
            queue_species=[r[c][0]]
        else:
            queue_species=[]
        while len(queue_species)>0:
            if TL_done :
                P_to_use=rec.lower_tree_computation.P_TL
            else:
                P_to_use=rec.lower_tree_computation.P
            e=queue_species.pop()

            if not mult_gene_match:
                bool_continue_mult_match=not (c.is_leaf() and c.match==e)
            else:
                bool_continue_mult_match=not (c.is_leaf() and e in c.match)

            if bool_continue_mult_match:
                #cas d'arrêt
                #on étend P[e][u] pour choisir ce qu'il s'est produit
                l_proba=[]#liste des proba des événéments
                l_event=[]#liste des événements, un evenement est un tuple, de taille 2 lorsqu'on fait se déplacer le parasite (lui et son nouveau match), 4 lorsqu'on passe à ses enfants (les deux enfants et leurs deux matchs), et 3 lors d'un TL, pour changer la valeur du compteur TL si choisi
                #liste des evenements possibles pour cet espece dans celle la
                if not e.isLeaf():
                    for (cL,cR) in c.log_child_frequencies:
                        l_proba.append(c.log_child_frequencies[(cL,cR)]+s_r+P[e.left][cL]+P[e.right][cR])
                        l_event.append(("S", e.left, cL, e.right, cR))
                        l_proba.append(c.log_child_frequencies[(cL,cR)]+s_r+P[e.left][cR]+P[e.right][cL])
                        l_event.append(("S", e.left, cR, e.right, cL))
                    if c in P_to_use[e.right]:
                        l_proba.append(s_r+E[e.left]+P_to_use[e.right][c])
                        l_event.append(("SL", e.right, c, e.left,c))
                    if c in P_to_use[e.left]:
                        l_proba.append(s_r+E[e.right]+P_to_use[e.left][c])
                        l_event.append(("SL", e.left, c, e.right, c))
                for (cL,cR) in c.log_child_frequencies:
                    l_proba.append(c.log_child_frequencies[(cL,cR)]+d_r+P[e][cL]+P[e][cR])
                    l_event.append(("D", e, cL, e,cR))
                if not P_transfer:
                    ancestor_e=e.ancestor_dict()
                for (cL,cR) in c.log_child_frequencies:
                    if P_transfer:
                        possible_transfer=[h for h in P_transfer[e].keys()]
                    else:
                        possible_transfer=[h for h in parasite_post_order if not h in ancestor_e]
                    for h in possible_transfer:
                        #l_proba.append(lclade_frequencies[c][(cL,cR)]+P_to_use[e][cR]+t_r+P_to_use[h][cL]+P_t[e][h])
                        if P_transfer:
                            l_proba.append(c.log_child_frequencies[(cL,cR)]+P[e][cR]+t_r+P[h][cL]+P_transfer[e][h])
                            l_proba.append(c.log_child_frequencies[(cL,cR)]+P[e][cL]+t_r+P[h][cR]+P_transfer[e][h])
                        else:
                            l_proba.append(c.log_child_frequencies[(cL,cR)]+P[e][cR]+t_r+P[h][cL]-np.log(len(parasite_post_order)-ancestrale_correction_size[e]))
                            l_proba.append(c.log_child_frequencies[(cL,cR)]+P[e][cL]+t_r+P[h][cR]-np.log(len(parasite_post_order)-ancestrale_correction_size[e]))
                        l_event.append(("T", h,cL,e,cR))
                        l_event.append(("T",h,cR,e,cL))
                if not TL_done:
                    if P_transfer:
                        possible_transfer=[h for h in P_transfer[e].keys()]
                    else:
                        possible_transfer=[h for h in parasite_post_order if not h in ancestor_e]
                    for h in possible_transfer:
                        if c in P_TL[h]:
                            if P_transfer:
                                l_proba.append(P_TL[h][c]+t_r+E[e]+P_transfer[e][h])
                            else:
                                l_proba.append(P_TL[h][c]+t_r+E[e]-np.log(len(parasite_post_order)-ancestrale_correction_size[e]))
                            l_event.append(("TL", h,c,e,c))
                if len(l_proba)==0:

                    print(len(P),len(P_TL),len(rec.single_lower.leaves),len(rec.upper.post_order), t_r,d_r,l_r,s_r)
                    for rec_event in l_sampled_events:
                        print(rec_event.name)

                    print("normalement pas besoin",TL_done, e.isLeaf(), c.match[0].name, e.name, P[e][c], c in P_TL[e])
                    event=("TL",e,c,c_match[c],c)
                else:
                    s=log_add_list(l_proba)
                    if best:
                        k=np.argmax(l_proba)
                    else:
                        x=s+np.log(rd.random())
                        x_current=l_proba[0]
                        k=0
                        while x_current<x:
                            k+=1
                            x_current=log_add(x_current, l_proba[k])
                    event=l_event[k]

                    in_sampling_scenario.log_likelihood+=l_proba[k]


                #creating the event
                event_name,f,v,g,w=event
                event=Rec_event()
                event.name=event_name
                event.upper=e
                event.lower=c
                if not event_name=="D":
                    event.upper_left_or_keeper_or_receiver=f
                    event.upper_right_or_loser_or_donor=g

                if event_name in ["S","T","D"]:
                    r[v].append(f)
                    r[w].append(g)
                    clade_to_look.append(v)
                    clade_to_look.append(w)
                    event.lower_left=v
                    event.lower_right=w

                    #adding new node to the reconstructed tree
                    reconstructed_tree_node.undated_birth()
                    reconstructed_tree_node.left.match=r[v]
                    reconstructed_tree_node.right.match=r[w]
                    reconstructed_tree_node.left.event_list=[]
                    reconstructed_tree_node.right.event_list=[]
                    reconstructed_tree_node.left.corresponding_am_tree=v
                    reconstructed_tree_node.right.corresponding_am_tree=w
                    reconstructed_tree_node.left.name=v.name
                    reconstructed_tree_node.right.name=w.name
                    am_tree_to_reconstructed[v]=reconstructed_tree_node.left
                    am_tree_to_reconstructed[w]=reconstructed_tree_node.right

                if "L" in event.name:
                    r[c].append(event.upper_left_or_keeper_or_receiver)
                    if event.name=="TL":
                        TL_done=True
                    queue_species.append(event.upper_left_or_keeper_or_receiver)

                event.init_name()
                l_sampled_events.append(event)
                if reconstructed_tree_node.event_list is None:
                    reconstructed_tree_node.event_list=[]
                else:
                    reconstructed_tree_node.event_list.append(event)

            else:
                #we got to a matching leaf
                event=Rec_event()
                event.name="E"
                event.upper=e
                event.lower=c
                reconstructed_tree_node.name=c.clade_leaves[0]
                l_sampled_events.append(event)

    return in_sampling_scenario
