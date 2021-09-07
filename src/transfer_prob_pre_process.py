#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:00:16 2020

@author: hmenet
"""


import numpy as np

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


def compute_transfer_prob(host_info, parasite_post_order, heuristic):

    #monte_carlo and decouple use the same process
    if heuristic in ["monte_carlo", "decouple", "sequential"]:
        P_transfer=prob_transfer_sequential(host_info, parasite_post_order)
    if heuristic == "unaware":
        P_transfer=prob_transfer_unaware(parasite_post_order)
    if heuristic == "alea":
        P_transfer=prob_transfer_alea(host_info, parasite_post_order)


    #for u in P_transfer:
    #    for v in P_transfer[u]:
    #        if P_transfer[u][v]==0:
    #            if u != v and (not v.isAscendant(u)):
    #               print(u.root.name, v.root.name, P_transfer[u][v])
    return P_transfer



### alea #####

def prob_transfer_alea(host_info, parasite_post_order):
    P_hp = host_info
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
def prob_transfer_sequential(host_info, parasite_post_order):
    host_post_order, rates_hp,match_hp, E,heuristic = host_info
    d_r= rates_hp["D"]
    l_r= rates_hp["L"]
    t_r= rates_hp["T"]
    s_r=1-d_r-l_r-t_r
    #on calcule d'abord la proba de transfert entre de gène de parasite sachant l'hôte de chacun des parasites
    P_transfer_h=dict()
    #P_transfer_h[h1][h2] : proba de transfert de gène entre un parasite placé dans l'hôte h1 et un dans l'hôte h2
    for h in host_post_order:
        P_transfer_h[h]=dict()
    #calcul du nombre de parasite dans chaque hôte
    match_hp_inv=dict()
    for p in match_hp.keys():
        for h in match_hp[p]:
            if not h in match_hp_inv:
                match_hp_inv[h]=[p]
            else:
                if not p in match_hp_inv[h]:
                    match_hp_inv[h].append(p)
    #N_parasites[e] : nombre de parasite dans l'hôte e
    N_parasites=dict()
    for h in match_hp_inv.keys():
        N_parasites[h]=len(match_hp_inv[h])

    if heuristic=="dec_no_ghost":
        for h in P_transfer_h:
            for h2 in host_post_order:
                P_transfer_h[h][h2]=0
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
    for p1 in parasite_post_order:
        P_transfer[p1]=dict()
        for p2 in parasite_post_order:
            if (not p1 == p2) and (not p2.isAscendant(p1)):
                p_transfertmp=sum([sum([P_transfer_h[h1][h2] for h1 in match_hp[p1]])for h2 in match_hp[p2]])
                if p_transfertmp != 0:
                    P_transfer[p1][p2]=np.log(p_transfertmp)
            #else:
            #    P_transfer[p1][p2]=0#pas de transfert vers soi même


    s=0
    n=0
    #for u in P_transfer:
    #    voisinage=[]
    #    for h in match_hp[u]:
    #        for p in match_hp_inv[h]:
    #            voisinage.append(p)
    #    for v in P_transfer[u]:
    #        if not v in voisinage:
    #            s+=np.exp(P_transfer[u][v])
    #            n+=1
    #print("ptransfer inter", s, s/n)
    return P_transfer
