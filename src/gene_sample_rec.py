#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:18:17 2020

@author: hmenet
"""

import random as rd
import numpy as np

from rec_aux_func import log_add, log_minus, log_add_list, is_mult_match


def sample_gene_upper_rec(P, P_TL,E, parasite_post_order, clades_data, rates, c_match, likelihood, ancestrale_correction_size, return_r=False, best=False, P_transfer=None):
    mult_gene_match=is_mult_match(c_match)
    clade_post_order, clade_frequencies, clade_elements, clade_keys=clades_data
    r=dict()
    for u in clade_post_order:
        r[u]=[]
    l_sampled_events=[]#list of all the events that compose the scenario we will have sampled at the end of this function
    #ce choix n'est pas tout à fait correct, mais représente une bonne approx
    gene =clade_post_order[0] #root des clades
    d_r= np.log(rates["D"])
    l_r= np.log(rates["L"])
    t_r= np.log(rates["T"])
    s_r=np.log(1-rates["D"]-rates["L"]-rates["T"])
    s=likelihood
    k=0
    x=np.log(rd.random())+s
    x_current=P[parasite_post_order[0]][gene]
    while x_current<x:
        k+=1
        x_current=log_add(x_current, P[parasite_post_order[k]][gene])
    #il faudrait des proba d'origination
    r[gene].append(parasite_post_order[k])

    #log
    lclade_frequencies=dict()
    for c in clade_frequencies:
        lclade_frequencies[c]=dict()
        for u in clade_frequencies[c]:
            lclade_frequencies[c][u]=np.log(clade_frequencies[c][u])


    clade_to_look=[gene]
    while len(clade_to_look)>0:
        c=clade_to_look.pop()

        #for c in clade_post_order:
        TL_done=False #on autorise un seul TL
        if len(r[c])>0:
            queue_species=[r[c][0]]
        else:
            queue_species=[]
        while len(queue_species)>0:
            if TL_done :
                P_to_use=P_TL
            else:
                P_to_use=P
            e=queue_species.pop()

            if not mult_gene_match:
                bool_continue_mult_match=not (len(clade_elements[c])==1 and c_match[c]==e)
            else:
                bool_continue_mult_match=not (len(clade_elements[c])==1 and e in c_match[c])

            if bool_continue_mult_match:
                #cas d'arrêt
                #on étend P[e][u] pour choisir ce qu'il s'est produit
                l_proba=[]#liste des proba des événéments
                l_event=[]#liste des événements, un evenement est un tuple, de taille 2 lorsqu'on fait se déplacer le parasite (lui et son nouveau match), 4 lorsqu'on passe à ses enfants (les deux enfants et leurs deux matchs), et 3 lors d'un TL, pour changer la valeur du compteur TL si choisi
                #liste des evenements possibles pour cet espece dans celle la
                if not e.isLeaf():
                    for (cL,cR) in clade_frequencies[c]:
                        l_proba.append(lclade_frequencies[c][(cL,cR)]+s_r+P[e.left][cL]+P[e.right][cR])
                        l_event.append(("S", e.left, cL, e.right, cR))
                        l_proba.append(lclade_frequencies[c][(cL,cR)]+s_r+P[e.left][cR]+P[e.right][cL])
                        l_event.append(("S", e.left, cR, e.right, cL))
                    if c in P_to_use[e.right]:
                        l_proba.append(s_r+E[e.left]+P_to_use[e.right][c])
                        l_event.append(("SL", e.right, c, e.left,c))
                    if c in P_to_use[e.left]:
                        l_proba.append(s_r+E[e.right]+P_to_use[e.left][c])
                        l_event.append(("SL", e.left, c, e.right, c))
                for (cL,cR) in clade_frequencies[c]:
                    l_proba.append(lclade_frequencies[c][(cL,cR)]+d_r+P[e][cL]+P[e][cR])
                    l_event.append(("D", e, cL, e,cR))
                if not P_transfer:
                    ancestor_e=e.ancestor_dict()
                for (cL,cR) in clade_frequencies[c]:
                    if P_transfer:
                        possible_transfer=[h for h in P_transfer[e].keys()]
                    else:
                        possible_transfer=[h for h in parasite_post_order if not h in ancestor_e]
                    for h in possible_transfer:
                        #l_proba.append(lclade_frequencies[c][(cL,cR)]+P_to_use[e][cR]+t_r+P_to_use[h][cL]+P_t[e][h])
                        if P_transfer:
                            l_proba.append(lclade_frequencies[c][(cL,cR)]+P[e][cR]+t_r+P[h][cL]+np.log(P_transfer[e][h]))
                            l_proba.append(lclade_frequencies[c][(cL,cR)]+P[e][cL]+t_r+P[h][cR]+np.log(P_transfer[e][h]))
                        else:
                            l_proba.append(lclade_frequencies[c][(cL,cR)]+P[e][cR]+t_r+P[h][cL]-np.log(len(parasite_post_order)-ancestrale_correction_size[e]))
                            l_proba.append(lclade_frequencies[c][(cL,cR)]+P[e][cL]+t_r+P[h][cR]-np.log(len(parasite_post_order)-ancestrale_correction_size[e]))
                        l_event.append(("T",e,cR, h,cL))
                        l_event.append(("T",e,cL, h,cR))
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
                            l_event.append(("TL", e,c,h,c))
                if len(l_proba)==0:
                    print("normalement pas besoin",TL_done, e.isLeaf(), c in c_match)
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

                if event[0]=="S":
                    event_name, f,v,g,w = event
                    r[v].append(f)
                    r[w].append(g)
                    clade_to_look.append(v)
                    clade_to_look.append(w)
                    event_to_append=(event_name,e,c,f,v,g,w)
                elif event[0]=="SL":
                    event_name, f, u,g,u1 = event
                    r[u].append(f)
                    queue_species.append(f)
                    event_to_append=(event_name,e,c,f,c,g,c)
                elif event[0]=="D":
                    event_name, e,v,e,w = event
                    r[v].append(e)
                    r[w].append(e)
                    clade_to_look.append(v)
                    clade_to_look.append(w)
                    event_to_append=(event_name,e,c,e,v,e,w)
                elif event[0]=="T":
                    event_name, e,v,h,w=event
                    r[v].append(e)
                    r[w].append(h)
                    clade_to_look.append(v)
                    clade_to_look.append(w)
                    event_to_append=(event_name,e,c,h,w,e,v)
                elif event[0]=="TL":
                    event_name,e,u,h,u=event
                    r[u].append(h)
                    TL_done = True
                    event_to_append=(event_name, e,c,h,c,e,c)
                    queue_species.append(h)
                l_sampled_events.append(event_to_append)
    #if return_r:
    return l_sampled_events, r
    #else:
    #    return l_sampled_events
