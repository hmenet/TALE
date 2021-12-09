#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:00:16 2020

@author: hmenet
"""


import numpy as np
import arbre


#####
### Compute transfer probabilities of genes between symbionts #####
#####




#calcul des proba de transferts de gène entre parasite

#rates_hp : taux evolution hote parasite
#match_hp: real match entre parasite et hôte, dict étant donné un parasite renvoie la liste de ses hôtes
#E : probabilité d'être perdu dans l'hôte
#heuristic is the name of the chosen heuristic, unaware, alea, sequential
#monte_carlo and decouple are sequential
#host_info depend on the heuristics
# alea : P_hp = host_info
# sequential: rates_hp,match_hp, E = host_info


def compute_transfer_prob(rec):
    heuristic=rec.heuristic
    #monte_carlo and decouple use the same process
    if heuristic in ["monte_carlo", "decouple", "sequential"]:
        P_transfer=prob_transfer_sequential(rec)
    if heuristic == "alea":
        P_transfer=prob_transfer_alea(rec)
    return P_transfer



### alea #####

def prob_transfer_alea(rec):
    P_hp = rec.upper_rec.P
    P_transfer=dict()
    #P_transfer[p1][p2] : proba de transférer entre les parasites p1 et p2
    for p1 in parasite_post_order:
        P_transfer[p1]=dict()
        for p2 in parasite_post_order:
            if (not p1 == p2) and (not p2.isAscendant(p1)):
                P_transfer[p1][p2]=sum([P_hp[e][p1] * P_hp[e][p2] for e in P_hp.keys()])
            else :
                P_transfer[p1][p2]=0
    return P_transfer

### sequential #####

#match_hp : clade format
# at the end we want tree format
#inter list is the reconstructed from amalgamated trees inter trees
#rec is the upper rec problem
def prob_transfer_sequential(rec):

    inter_list=rec.upper
    host_post_order=rec.upper_rec.upper.post_order

    E=rec.upper_rec.upper_tree_computation.E_no_log
    d_r= rec.upper_rec.rates.dr
    l_r= rec.upper_rec.rates.lr
    t_r= rec.upper_rec.rates.tr
    s_r=rec.upper_rec.rates.sr
    #on calcule d'abord la proba de transfert entre de gène de parasite sachant l'hôte de chacun des parasites
    P_transfer_h=dict()
    #P_transfer_h[h1][h2] : proba de transfert de gène entre un parasite placé dans l'hôte h1 et un dans l'hôte h2
    for h in host_post_order:
        P_transfer_h[h]=dict()
    #calcul du nombre de parasite dans chaque hôte
    match_hp_inv=dict()
    for p in inter_list.post_order:
        for h in p.match:
            if not h in match_hp_inv:
                match_hp_inv[h]=[p]
            else:
                if not p in match_hp_inv[h]:
                    match_hp_inv[h].append(p)
    #N_parasites[e] : nombre de parasite dans l'hôte e
    N_parasites=dict()
    for h in match_hp_inv.keys():
        N_parasites[h]=len(match_hp_inv[h])

    if rec.heuristic=="dec_no_ghost":
        for h in P_transfer_h:
            for h2 in host_post_order:
                P_transfer_h[h][h2]=0
        for h in match_hp_inv.keys():
            P_transfer_h[h][h]=1/N_parasites[h]
    else:
        if rec.heuristic=="dd_dec":
            P_transfer_upper=distance_dependent(rec.upper_rec)
            for h in P_transfer_h:
                for h2 in P_transfer_upper[h]:
                    if h2 in N_parasites:
                        P_transfer_h[h][h2]=np.exp(P_transfer_upper[h][h2])/N_parasites[h2]
            for h in match_hp_inv.keys():
                P_transfer_h[h][h]=1/N_parasites[h]
        else:
            for target_e in host_post_order:
                if target_e in match_hp_inv:#si on ne possède aucun parasite pas de calcul
                    P=dict()
                    P_TL=dict()
                    #target_e # where we want to transfer
                    host_queue=list(host_post_order)
                    while len(host_queue)> 0 :
                        e=host_queue.pop()
                        a=0
                        b=0
                        if e==target_e:
                            #N_parasites[e] : nombre de parasite dans l'hôte e
                            b+=1/N_parasites[e]
                        if not e.isLeaf():
                            b+=s_r*(P_TL[e.left]*E[e.right] + P_TL[e.right]*E[e.left])
                        a+=1-2*d_r*E[e]
                        P_TL[e]=b/a
                    host_queue=list(host_post_order)
                    while len(host_queue)> 0 :
                        e=host_queue.pop()
                        #if not target_e.isFreeLiving:
                        a=0
                        b=0
                        if e==target_e:
                            #N_parasites[e] : nombre de parasite dans l'hôte e
                            b+=1/N_parasites[e]
                        if not e.isLeaf():
                            b+=s_r*(P[e.left]*E[e.right] + P[e.right]*E[e.left])
                        N_h=0
                        tmp=0
                        for h in host_post_order:
                            if not h.isAscendant(e) and not h==e:
                                tmp+=P_TL[h]
                                N_h+=1
                        b+=tmp*t_r*E[e]/N_h
                        a+=1-2*d_r*E[e]
                        P[e]=b/a
                    for e in host_post_order:
                        P_transfer_h[e][target_e]=P[e]
    P_transfer=dict()
    #P_transfer[p1][p2] : proba de transférer entre les parasites p1 et p2

    for p1 in inter_list.post_order:
        P_transfer[p1]=dict()
        norm_factor=0
        for p2 in inter_list.post_order:
            if (not p1 == p2) and (not p2.isAscendant(p1)):
                if rec.heuristic=="dd_dec":
                    p_transfertmp=0
                    for h1 in p1.match:
                        for h2 in p2.match:
                            if h2 in P_transfer_h[h1]:
                                p_transfertmp+=P_transfer_h[h1][h2]
                else:
                    p_transfertmp=sum([sum([P_transfer_h[h1][h2] for h1 in p1.match])for h2 in p2.match])

                if p_transfertmp != 0:
                    P_transfer[p1][p2]=p_transfertmp
                    norm_factor+=p_transfertmp
        for p2 in P_transfer[p1]:
            P_transfer[p1][p2]=np.log(P_transfer[p1][p2]/norm_factor)

    return P_transfer

def distance_dependent(rec):

    alpha=1

    P_transfer=dict()
    host_post_order=rec.upper.post_order
    t_r=rec.rates.tr
    #dans ce cas on ne compte pas les i dans la liste des transferts possibles
    if rec.rates.ir != 0:
        ir=True
    else:
        ir=False
    for e1 in host_post_order:
        s=0
        P_transfer_tmp=dict()
        P_transfer[e1]=dict()
        for e2 in host_post_order:
            dse1e2=arbre.distance(e1,e2)
            if (not e1 == e2) and (not e2.isAscendant(e1)):
                if not ir or e2.isRoot() or not e2.parent==e1:
                    P_transfer_tmp[e2]=t_r**(dse1e2*alpha)
                    s+=P_transfer_tmp[e2]
        for e2 in P_transfer_tmp:
            P_transfer[e1][e2]=np.log(P_transfer_tmp[e2]/s)
    #for e1 in host_post_order:
    #    #print([np.exp(P_transfer[e1][e2]) for e2 in P_transfer[e1]])
    #   print(sum([np.exp(P_transfer[e1][e2]) for e2 in P_transfer[e1]]))
    return P_transfer


