#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:28:33 2020

@author: hmenet
"""

import numpy as np

from rec_aux_func import log_add, log_minus, log_add_list, is_mult_match, inter_list


################################
#vraisemblance rec symbiote gene
################################

#return the positive solution we're seeking
def solution_polynome_snd(a,b,c):
    delta=b**2-4*a*c
    return (-b - delta**0.5)/(2*a)


#pareil que TL_compute_parasite_E_5_order mais peut prendre une lsite de parasite en entrée en donnant la concatenation des post order de ces parasites
def compute_upper_gene_E(upper_post_order, rates, P_transfer=None):
    d_r= rates["D"]
    l_r= rates["L"]
    t_r= rates["T"]
    s_r=1-d_r-l_r-t_r
    E=dict()
    queue=list(upper_post_order)
    Eavg=dict()
    for h in upper_post_order:
        tmp=l_r
        tmp+=d_r*l_r*l_r
        if not h.isLeaf():
            tmp+=s_r*l_r*l_r
        if not h.isRoot():
            tmp+=l_r*l_r*t_r
        Eavg[h]=tmp
    Eavg1=sum([t_r*Eavg[h] for h in upper_post_order])
    while len(queue)>0:
        e=queue.pop()
        a=d_r
        #b= -1 + t_r*sum([P_t[e][h]*Eavg[h] for h in P_t[e].keys()])

        l_tmp=e.ancestor()
        correction_ancestral=sum([t_r*Eavg[h] for h in l_tmp])
        b=-1 + (Eavg1-correction_ancestral)/(len(upper_post_order)-len(l_tmp))

        c=l_r
        if not e.isLeaf():
            c+=s_r*E[e.left]*E[e.right]
        E[e]=solution_polynome_snd(a,b,c)

    Eavg_no_log=(-1)/len(upper_post_order)*sum([t_r*E[h] for h in upper_post_order])
    for h in E.keys():
        E[h]=np.log(E[h])

    return E, Eavg_no_log


#proba de transfert pre processed, dans P_transfer[e][e2] de e vers e2
#c_match = clade_to_matched_node
def compute_upper_gene_P(upper_post_order,E, Eavg_no_log, clades_data, rates, c_match, more_output=False, P_transfer=None):
    mult_gene_match=is_mult_match(c_match)
    d_r= np.log(rates["D"])
    l_r= np.log(rates["L"])
    t_r= np.log(rates["T"])
    s_r=np.log(1-rates["D"]-rates["L"]-rates["T"])
    clade_post_order, clade_frequencies, clade_elements, clade_keys=clades_data
    P=dict()
    P_TL=dict()
    P_avg=dict()
    P_avg_TL=dict()
    correction_ancestrale=dict()
    correction_ancestrale_TL=dict()
    correction_ancestrale_size=dict()

    for e in upper_post_order:
        P[e]=dict()
        P_TL[e]=dict()
        correction_ancestrale[e]=dict()
        correction_ancestrale_TL[e]=dict()
    if not P_transfer:
        for e in upper_post_order:
            if e.isRoot():
                correction_ancestrale_size[e]=1
            else:
                correction_ancestrale_size[e]=correction_ancestrale_size[e.parent]

    if not mult_gene_match:
        for c in c_match:
            #print(c, c_match[c].root.name, c_match[c].root)
            P[c_match[c]][c]=0
            P_TL[c_match[c]][c]=0
    else:
        for c in c_match:
            #c_match[c] si mult match c'est une liste
            for species_tmp in c_match[c]:
                P[species_tmp][c]=(-1)*np.log(len(c_match[c]))
                P_TL[species_tmp][c]=(-1)*np.log(len(c_match[c]))
    clade_queue=list(clade_post_order)


    #for c in clade_post_order:
    #    for e in P.keys():
    #        if not c in P[e]:
    #            P[e][c]=0
    #for c in clade_post_order:
    #    for u in clades_data[2][c]:
    #        if sum([P[e][clades_data[3][tuple([u])]] for e in P.keys()]) != 1 :
    #            print("hey",u,clades_data[3][tuple([u])], sum([P[e][clades_data[3][tuple([u])]] for e in P.keys()]))

    while len(clade_queue)>0:
        c=clade_queue.pop()
        upper_queue=list(upper_post_order)
        a_store=dict()
        b_store=dict()
        while len(upper_queue)>0:
            e=upper_queue.pop()
            if not mult_gene_match:
                bool_continue_mult_match=not (c in c_match.keys() and c_match[c]==e)
                bool_continue_mult_match2=not(c in c_match.keys()) or (c_match[c] in e.leaves())
            else:
                bool_continue_mult_match=not (c in c_match.keys() and e in c_match[c])
                bool_continue_mult_match2=not(c in c_match.keys()) or (inter_list(c_match[c],e.leaves()))
            if bool_continue_mult_match:
                if bool_continue_mult_match2:
                    #a=1
                    #a+=(-2)*d_r*E[e]
                    #a+=Eavg
                    #b=0
                    a=1
                    a+=(-2)*np.exp(d_r)*np.exp(E[e])
                    a+=Eavg_no_log
                    b_l=[]
                    if not e.isLeaf():
                        for cL,cR in clade_frequencies[c]:

                            b1=log_add(P[e.left][cL]+P[e.right][cR],P[e.left][cR]+P[e.right][cL])
                            b1+=np.log(clade_frequencies[c][(cL,cR)])


                            #b+=(P[e.left][cL]*P[e.right][cR]+P[e.left][cR]*P[e.right][cL])*clade_frequencies[c][(cL,cR)]
                            b_l.append(b1)
                        if c in c_match.keys():
                            if mult_gene_match:
                                c_match_c=c_match[c]
                            else:
                                c_match_c=[c_match[c]]
                            if inter_list(c_match_c, e.right.leaves()):
                                b21=E[e.left]+P_TL[e.right][c]+s_r
                                b_l.append(b21)
                            if inter_list(c_match_c, e.left.leaves()):
                                b22=E[e.right]+P_TL[e.left][c]+s_r
                                b_l.append(b22)

                        else:
                            b2=log_add(E[e.left]+P_TL[e.right][c],E[e.right]+P_TL[e.left][c])
                            b2+=s_r
                            b_l.append(b2)

                        #b+=E[e.left]*P_TL[e.right][c]
                        #b+=E[e.right]*P_TL[e.left][c]
                        #b*=s_r

                    for cL,cR in clade_frequencies[c]:
                        if P_transfer:
                            for h in P_transfer[e].keys():
                                b3=np.log(clade_frequencies[c][(cL,cR)])+t_r+np.log(P_transfer[e][h])+P[h][cL]+P[e][cR]
                                b4=np.log(clade_frequencies[c][(cL,cR)])+t_r+np.log(P_transfer[e][h])+P[h][cR]+P[e][cL]
                                b_l.append(b3)
                                b_l.append(b4)
                        else:
                            b3=np.log(clade_frequencies[c][(cL,cR)])+t_r-np.log(len(upper_post_order)-correction_ancestrale_size[e])+log_minus(P_avg[cL],correction_ancestrale[e][cL])+P[e][cR]
                            b4=np.log(clade_frequencies[c][(cL,cR)])+t_r-np.log(len(upper_post_order)-correction_ancestrale_size[e])+log_minus(P_avg[cR],correction_ancestrale[e][cR])+P[e][cL]
                            b_l.append(b3)
                            b_l.append(b4)

                        #b+=clade_frequencies[c][(cL,cR)]*t_r*sum([P_transfer[e][h]*P[h][cL] for h in P_transfer[e].keys()])*P[e][cR]
                        #b+=clade_frequencies[c][(cL,cR)]*t_r*sum([P_transfer[e][h]*P[h][cR] for h in P_transfer[e].keys()])*P[e][cL]
                        #duplication

                        b6=d_r+P[e][cL]+P[e][cR]+np.log(clade_frequencies[c][(cL,cR)])
                        b_l.append(b6)
                        #b+=d_r*P[e][cL]*P[e][cR]*clade_frequencies[c][(cL,cR)]
                    b=log_add_list(b_l)

                    resultat=b-np.log(a)
                    P_TL[e][c]=resultat
                    #P_avg_TL[e]=log_add(P_avg_TL[e],resultat)
                    a_store[e]=a
                    b_store[e]=b
        if c in c_match.keys():
            transfer_list_e=[]
            for e in upper_post_order:
                if c in P_TL[e]:
                    transfer_list_e.append(e)
            P_avg_TL[c]=log_add_list([P_TL[e][c] for e in transfer_list_e])
        else:
            P_avg_TL[c]=log_add_list([P_TL[e][c] for e in upper_post_order])
        if c in c_match.keys():
            list_to_traverse=transfer_list_e
        else:
            list_to_traverse=upper_post_order
        if not P_transfer:
            for e in list_to_traverse:
                if e.isRoot():
                    correction_ancestrale_TL[e][c]=P_TL[e][c]
                else:
                    correction_ancestrale_TL[e][c]=log_add(P_TL[e][c],correction_ancestrale_TL[e.parent][c])
        upper_queue=list(upper_post_order)
        while len(upper_queue)>0:
            e=upper_queue.pop()

            if not mult_gene_match:
                bool_continue_mult_match=not (c in c_match.keys() and c_match[c]==e)
            else:
                bool_continue_mult_match=not (c in c_match.keys() and e in c_match[c])
            if bool_continue_mult_match:
                if e in b_store:
                    b_l=[b_store[e]]
                    a=a_store[e]
                    #b+=t_r*sum([P_transfer[e][h]*P_TL[h][c] for h in P_transfer[e].keys()])*E[e]
                    if P_transfer:
                        for h in P_transfer[e].keys():
                            b1=t_r+np.log(P_transfer[e][h])+P_TL[h][c]+E[e]
                            b_l.append(b1)
                    else:
                        b1=t_r-np.log(len(upper_post_order)-correction_ancestrale_size[e])+log_minus(P_avg_TL[c],correction_ancestrale_TL[e][c])+E[e]
                        b_l.append(b1)


                    if not e.isLeaf():
                        #b+=s_r*(E[e.left]*(P[e.right][c]-P_TL[e.right][c] )  + E[e.right]*(P[e.left][c]-P_TL[e.left][c]))
                        if c in P_TL[e.right]:
                            if not P[e.right][c]==P_TL[e.right][c]:
                                b2=s_r+E[e.left]+log_minus(P[e.right][c],P_TL[e.right][c])
                                b_l.append(b2)
                        else:
                            b2=s_r+E[e.left]+P[e.right][c]
                            b_l.append(b2)
                        if c in P_TL[e.left]:
                            if not P[e.left][c]==P_TL[e.left][c]:
                                b3=s_r+E[e.right]+log_minus(P[e.left][c],P_TL[e.left][c])
                                b_l.append(b3)
                        else:
                            b3=s_r+E[e.right]+P[e.left][c]
                            b_l.append(b3)
                    b=log_add_list(b_l)
                    resultat=b-np.log(a)
                    P[e][c]=resultat
                else:
                    a=1
                    a+=(-2)*np.exp(d_r)*np.exp(E[e])
                    a+=Eavg_no_log
                    b_l=[]
                    if P_transfer:
                        for h in P_transfer[e].keys():
                            b1=t_r+np.log(P_transfer[e][h])+P_TL[h][c]+E[e]
                            b_l.append(b1)
                    else:
                        b1=t_r-np.log(len(upper_post_order)-correction_ancestrale_size[e])+P_avg_TL[c]+E[e]
                        b_l.append(b1)
                    if not e.isLeaf():
                        b2=s_r+E[e.left]+P[e.right][c]
                        b3=s_r+E[e.right]+P[e.left][c]
                        b_l.append(b2)
                        b_l.append(b3)
                    b=log_add_list(b_l)
                    resultat=b-np.log(a)
                    P[e][c]=resultat
        P_avg[c]=log_add_list([P[e][c] for e in upper_post_order])
        if not P_transfer:
            for e in upper_post_order:
                if e.isRoot():
                    correction_ancestrale[e][c]=P[e][c]
                else:
                    correction_ancestrale[e][c]=log_add(P[e][c],correction_ancestrale[e.parent][c])


    if more_output:
        likelihood=log_add_list([P[k][0] for k in P.keys()])
        return P, P_TL,likelihood, correction_ancestrale_size
    else:
        return P