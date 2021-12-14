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
def compute_upper_gene_E(rec, P_transfer=None):
    d_r= rec.rates.dr
    l_r= rec.rates.lr
    t_r= rec.rates.tr
    rec.rates.reinit()
    s_r=rec.rates.sr
    E=dict()
    queue=list(rec.upper.post_order)
    Eavg=dict()
    upper_post_order=rec.upper.post_order
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
        l_tmp=e.ancestor()
        correction_ancestral=sum([t_r*Eavg[h] for h in l_tmp])
        b=-1 + (Eavg1-correction_ancestral)/(len(upper_post_order)-len(l_tmp))
        c=l_r
        if not e.isLeaf():
            c+=s_r*E[e.left]*E[e.right]
        E[e]=solution_polynome_snd(a,b,c)
    Eavg_no_log=(-1)/len(upper_post_order)*sum([t_r*E[h] for h in upper_post_order])
    E_log=dict()
    for h in E.keys():
        E_log[h]=np.log(E[h])
    #return
    rec.upper_tree_computation.Eavg_no_log=Eavg_no_log
    rec.upper_tree_computation.E_no_log=E
    rec.upper_tree_computation.E=E_log


def compute_upper_gene_P(rec_problem):
    if rec_problem.slm=="joint_ml":
        return compute_upper_gene_P_joint_ml(rec_problem)
    if rec_problem.slm=="tree_ml":
        return compute_upper_gene_P_ml_tree(rec_problem)
    else:
        return compute_upper_gene_P_l(rec_problem)

#proba de transfert pre processed, dans P_transfer[e][e2] de e vers e2
#c_match = clade_to_matched_node
def compute_upper_gene_P_l(rec_problem):

    E=rec_problem.upper_tree_computation.E
    Eavg_no_log=rec_problem.upper_tree_computation.Eavg_no_log
    am_tree=rec_problem.single_lower
    P_transfer=rec_problem.upper_tree_computation.P_transfer

    mult_gene_match=is_mult_match(am_tree)
    rec_problem.rates.reinit()
    d_r= rec_problem.rates.ldr
    l_r= rec_problem.rates.llr
    t_r= rec_problem.rates.ltr
    s_r= rec_problem.rates.lsr
    i_r=rec_problem.rates.lir
    upper_post_order=rec_problem.upper.post_order

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
                correction_ancestrale_size[e]=correction_ancestrale_size[e.parent]+1

    if not mult_gene_match:
        for c in am_tree.leaves:
            P[c.match][c]=0
            P_TL[c.match][c]=0
    else:
        for c in am_tree.leaves:
            #c_match[c] si mult match c'est une liste
            for species_tmp in c.match:
                P[species_tmp][c]=(-1)*np.log(len(c.match))
                P_TL[species_tmp][c]=(-1)*np.log(len(c.match))

    clade_queue=list(am_tree.reverse_post_order)

    while len(clade_queue)>0:
        c=clade_queue.pop()
        upper_queue=list(upper_post_order)
        a_store=dict()
        b_store=dict()
        #reachable_notl=notl_reachable(c,mult_gene_match)
        while len(upper_queue)>0:
            e=upper_queue.pop()
            if not mult_gene_match:
                bool_continue_mult_match=not (c.match==e)
                bool_continue_mult_match2=not(c.is_leaf()) or (c.match in e.leaves())
            else:
                bool_continue_mult_match=not (c.is_leaf() and e in c.match)
                bool_continue_mult_match2=not(c.is_leaf()) or (inter_list(c.match,e.leaves()))
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
                        for cL,cR in c.log_child_frequencies:

                            b1=log_add(P[e.left][cL]+P[e.right][cR],P[e.left][cR]+P[e.right][cL])
                            b1+=c.log_child_frequencies[(cL,cR)]+s_r


                            #b+=(P[e.left][cL]*P[e.right][cR]+P[e.left][cR]*P[e.right][cL])*clade_frequencies[c][(cL,cR)]
                            b_l.append(b1)
                        if c.is_leaf():
                            if mult_gene_match:
                                c_match_c=c.match
                            else:
                                c_match_c=[c.match]
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

                    #incomplete sorting event
                    if not i_r is None:
                        if not e.isLeaf():
                            for cL,cR in c.log_child_frequencies:

                                b7=P[e.left][cL]+P[e][cR]+c.log_child_frequencies[(cL,cR)]+i_r
                                b8=P[e.left][cR]+P[e][cL]+c.log_child_frequencies[(cL,cR)]+i_r
                                b9=P[e.right][cL]+P[e][cR]+c.log_child_frequencies[(cL,cR)]+i_r
                                b10=P[e.right][cR]+P[e][cL]+c.log_child_frequencies[(cL,cR)]+i_r

                                b_l.append(b7)
                                b_l.append(b8)
                                b_l.append(b9)
                                b_l.append(b10)


                    for cL,cR in c.log_child_frequencies:
                        if P_transfer:
                            for h in P_transfer[e].keys():
                                b3=c.log_child_frequencies[(cL,cR)]+t_r+P_transfer[e][h]+P[h][cL]+P[e][cR]
                                b4=c.log_child_frequencies[(cL,cR)]+t_r+P_transfer[e][h]+P[h][cR]+P[e][cL]
                                b_l.append(b3)
                                b_l.append(b4)
                        else:
                            if P_avg[cL]>correction_ancestrale[e][cL]:
                                b3=c.log_child_frequencies[(cL,cR)]+t_r-np.log(len(upper_post_order)-correction_ancestrale_size[e])+log_minus(P_avg[cL],correction_ancestrale[e][cL])+P[e][cR]
                                b_l.append(b3)
                            if P_avg[cR]>correction_ancestrale[e][cR]:
                                b4=c.log_child_frequencies[(cL,cR)]+t_r-np.log(len(upper_post_order)-correction_ancestrale_size[e])+log_minus(P_avg[cR],correction_ancestrale[e][cR])+P[e][cL]
                                b_l.append(b4)

                        #b+=clade_frequencies[c][(cL,cR)]*t_r*sum([P_transfer[e][h]*P[h][cL] for h in P_transfer[e].keys()])*P[e][cR]
                        #b+=clade_frequencies[c][(cL,cR)]*t_r*sum([P_transfer[e][h]*P[h][cR] for h in P_transfer[e].keys()])*P[e][cL]
                        #duplication

                        b6=d_r+P[e][cL]+P[e][cR]+c.log_child_frequencies[(cL,cR)]
                        b_l.append(b6)
                        #b+=d_r*P[e][cL]*P[e][cR]*clade_frequencies[c][(cL,cR)]

                    #print(len(b_l),e.name,c.name)
                    b=log_add_list(b_l)

                    resultat=b-np.log(a)
                    P_TL[e][c]=resultat
                    #P_avg_TL[e]=log_add(P_avg_TL[e],resultat)
                    a_store[e]=a
                    b_store[e]=b
        if c.is_leaf():
            transfer_list_e=[]
            for e in upper_post_order:
                if c in P_TL[e]:
                    transfer_list_e.append(e)
            P_avg_TL[c]=log_add_list([P_TL[e][c] for e in transfer_list_e])
        else:
            P_avg_TL[c]=log_add_list([P_TL[e][c] for e in upper_post_order])
        if c.is_leaf():
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
                bool_continue_mult_match=not (c.is_leaf() and c.match==e)
            else:
                bool_continue_mult_match=not (c.is_leaf() and e in c.match)
            if bool_continue_mult_match:
                if e in b_store:
                    b_l=[b_store[e]]
                    a=a_store[e]
                    #b+=t_r*sum([P_transfer[e][h]*P_TL[h][c] for h in P_transfer[e].keys()])*E[e]
                    if P_transfer:
                        for h in P_transfer[e].keys():
                            if c in P_TL[h]:
                                b1=t_r+P_transfer[e][h]+P_TL[h][c]+E[e]
                                b_l.append(b1)
                    else:
                        #if P_avg_TL[c]==correction_ancestrale_TL[e][c]:
                        #    #print(P_avg_TL[c],correction_ancestrale_TL[e][c])
                        #    #print(c.name,e.name,e.parent.name,len(c.clade_leaves),correction_ancestrale_TL[e.parent][c],P_TL[e.parent][c],P_TL[e][c],P_avg_TL[c],log_add_list([P_TL[e][c] for e in upper_post_order]))
                        #    #for e in upper_post_order:
                        #    #    print(e.name, P_TL[e][c])correction_ancestrale_size
                        #    #print(log_add(P_TL[upper_post_order[0]][c],P_TL[upper_post_order[-1]][c]),P_TL[upper_post_order[0]][c],P_TL[upper_post_order[-1]][c])
                        if P_avg_TL[c]>correction_ancestrale_TL[e][c]:
                            if None in [t_r,np.log(len(upper_post_order)-correction_ancestrale_size[e]),log_minus(P_avg_TL[c],correction_ancestrale_TL[e][c]),E[e]]:
                                print([t_r,np.log(len(upper_post_order)-correction_ancestrale_size[e]),P_avg_TL[c],correction_ancestrale_TL[e][c],E[e]])
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
                            if c in P_TL[h]:
                                b1=t_r+P_transfer[e][h]+P_TL[h][c]+E[e]
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


    likelihood=log_add_list([P[k][am_tree] for k in P.keys()])
    return P, P_TL,likelihood, correction_ancestrale_size



#proba de transfert pre processed, dans P_transfer[e][e2] de e vers e2
#new version to account for sampling ml scenario and ml amalgamated tree
def compute_upper_gene_P_joint_ml(rec_problem):

    E=rec_problem.upper_tree_computation.E
    Eavg_no_log=rec_problem.upper_tree_computation.Eavg_no_log
    am_tree=rec_problem.single_lower
    P_transfer=rec_problem.upper_tree_computation.P_transfer

    mult_gene_match=is_mult_match(am_tree)
    rec_problem.rates.reinit()
    d_r= rec_problem.rates.ldr
    l_r= rec_problem.rates.llr
    t_r= rec_problem.rates.ltr
    s_r= rec_problem.rates.lsr
    i_r=rec_problem.rates.lir
    upper_post_order=rec_problem.upper.post_order

    P=dict()
    P_TL=dict()
    P_max=dict()
    P_max_TL=dict()
    correction_ancestrale=dict()
    correction_ancestrale_TL=dict()
    correction_ancestrale_size=dict()

    if not P_transfer:
        for e in upper_post_order:
            if e.isRoot():
                correction_ancestrale_size[e]=1
            else:
                correction_ancestrale_size[e]=correction_ancestrale_size[e.parent]


    for e in upper_post_order:
        P[e]=dict()
        P_TL[e]=dict()

    if not mult_gene_match:
        for c in am_tree.leaves:
            P[c.match][c]=0
            P_TL[c.match][c]=0
    else:
        for c in am_tree.leaves:
            #c_match[c] si mult match c'est une liste
            for species_tmp in c.match:
                P[species_tmp][c]=(-1)*np.log(len(c.match))
                P_TL[species_tmp][c]=(-1)*np.log(len(c.match))

    clade_queue=list(am_tree.reverse_post_order)

    while len(clade_queue)>0:
        c=clade_queue.pop()
        upper_queue=list(upper_post_order)
        b_store=dict()
        #reachable_notl=notl_reachable(c,mult_gene_match)
        while len(upper_queue)>0:
            e=upper_queue.pop()
            if not mult_gene_match:
                bool_continue_mult_match=not (c.match==e)
                bool_continue_mult_match2=not(c.is_leaf()) or (c.match in e.leaves())
            else:
                bool_continue_mult_match=not (c.is_leaf() and e in c.match)
                bool_continue_mult_match2=not(c.is_leaf()) or (inter_list(c.match,e.leaves()))
            if bool_continue_mult_match:
                if bool_continue_mult_match2:
                    b_l=[]
                    if not e.isLeaf():
                        for cL,cR in c.log_child_frequencies:
                            b_l.append(P[e.left][cL]+P[e.right][cR]+c.log_child_frequencies[(cL,cR)]+s_r)
                            b_l.append(P[e.left][cR]+P[e.right][cL]+c.log_child_frequencies[(cL,cR)]+s_r)
                        if c.is_leaf():
                            if mult_gene_match:
                                c_match_c=c.match
                            else:
                                c_match_c=[c.match]
                            if inter_list(c_match_c, e.right.leaves()):
                                b_l.append(E[e.left]+P_TL[e.right][c]+s_r)
                            if inter_list(c_match_c, e.left.leaves()):
                                b_l.append(E[e.right]+P_TL[e.left][c]+s_r)
                        else:
                            b_l.append(E[e.left]+P_TL[e.right][c]+s_r)
                            b_l.append(E[e.right]+P_TL[e.left][c]+s_r)
                    #incomplete sorting event
                    if not i_r is None:
                        if not e.isLeaf():
                            for cL,cR in c.log_child_frequencies:

                                b7=P[e.left][cL]+P[e][cR]+c.log_child_frequencies[(cL,cR)]+i_r
                                b8=P[e.left][cR]+P[e][cL]+c.log_child_frequencies[(cL,cR)]+i_r
                                b9=P[e.right][cL]+P[e][cR]+c.log_child_frequencies[(cL,cR)]+i_r
                                b10=P[e.right][cR]+P[e][cL]+c.log_child_frequencies[(cL,cR)]+i_r

                                b_l.append(b7)
                                b_l.append(b8)
                                b_l.append(b9)
                                b_l.append(b10)


                    for cL,cR in c.log_child_frequencies:
                        if P_transfer:
                            for h in P_transfer[e].keys():
                                b3=c.log_child_frequencies[(cL,cR)]+t_r+P_transfer[e][h]+P[h][cL]+P[e][cR]
                                b4=c.log_child_frequencies[(cL,cR)]+t_r+P_transfer[e][h]+P[h][cR]+P[e][cL]
                                b_l.append(b3)
                                b_l.append(b4)
                        else:


                            b3=c.log_child_frequencies[(cL,cR)]+t_r-np.log(len(upper_post_order))+P_max[cL]+P[e][cR]
                            b_l.append(b3)
                            b4=c.log_child_frequencies[(cL,cR)]+t_r-np.log(len(upper_post_order))+P_max[cR]+P[e][cL]
                            b_l.append(b4)




                        #duplication

                        b6=d_r+P[e][cL]+P[e][cR]+c.log_child_frequencies[(cL,cR)]
                        b_l.append(b6)

                    b=max(b_l)
                    P_TL[e][c]=b
        if c.is_leaf():
            transfer_list_e=[]
            q=list(upper_post_order)
            for e in upper_post_order:
                if c in P_TL[e]:
                    transfer_list_e.append(e)
                    if not c in P_max_TL:
                        P_max_TL[c]=P_TL[e][c]
                    else:
                        P_max_TL[c]=max(P_max_TL[c],P_TL[e][c])
        else:
            P_max_TL[c]=max([P_TL[e][c] for e in upper_post_order])
        upper_queue=list(upper_post_order)
        while len(upper_queue)>0:
            e=upper_queue.pop()
            if not mult_gene_match:
                bool_continue_mult_match=not (c.is_leaf() and c.match==e)
            else:
                bool_continue_mult_match=not (c.is_leaf() and e in c.match)
            if bool_continue_mult_match:
                if e in b_store:
                    b_l=[b_store[e]]
                    if P_transfer:
                        for h in P_transfer[e].keys():
                            if c in P_TL[h]:
                                b1=t_r+P_transfer[e][h]+P_TL[h][c]+E[e]
                                b_l.append(b1)
                    else:
                        b1=t_r-np.log(len(upper_post_order))+(P_max_TL[c])+E[e]
                        b_l.append(b1)

                    if not e.isLeaf():
                        if c in P_TL[e.right]:
                            if not P[e.right][c]==P_TL[e.right][c]:
                                b2=s_r+E[e.left]+P[e.right][c]
                                b_l.append(b2)
                        else:
                            b2=s_r+E[e.left]+P[e.right][c]
                            b_l.append(b2)
                        if c in P_TL[e.left]:
                            if not P[e.left][c]==P_TL[e.left][c]:
                                b3=s_r+E[e.right]+P[e.left][c]
                                b_l.append(b3)
                        else:
                            b3=s_r+E[e.right]+P[e.left][c]
                            b_l.append(b3)
                    b=max(b_l)
                    P[e][c]=b
                else:
                    b_l=[]
                    if P_transfer:
                        for h in P_transfer[e].keys():
                            if c in P_TL[h]:
                                b1=t_r+P_transfer[e][h]+P_TL[h][c]+E[e]
                                b_l.append(b1)
                    else:
                        b1=t_r-np.log(len(upper_post_order))+P_max_TL[c]+E[e]
                        b_l.append(b1)
                    if not e.isLeaf():
                        b2=s_r+E[e.left]+P[e.right][c]
                        b3=s_r+E[e.right]+P[e.left][c]
                        b_l.append(b2)
                        b_l.append(b3)
                    b=max(b_l)
                    P[e][c]=b
        P_max[c]=max([P[e][c] for e in upper_post_order])
    likelihood=max([P[k][am_tree] for k in P.keys()])
    return P, P_TL,likelihood, correction_ancestrale_size

#proba de transfert pre processed, dans P_transfer[e][e2] de e vers e2
#new version to account for sampling ml scenario and ml amalgamated tree
def compute_upper_gene_P_ml_tree(rec_problem):

    E=rec_problem.upper_tree_computation.E
    Eavg_no_log=rec_problem.upper_tree_computation.Eavg_no_log
    am_tree=rec_problem.single_lower
    P_transfer=rec_problem.upper_tree_computation.P_transfer

    mult_gene_match=is_mult_match(am_tree)
    rec_problem.rates.reinit()
    d_r= rec_problem.rates.ldr
    l_r= rec_problem.rates.llr
    t_r= rec_problem.rates.ltr
    s_r= rec_problem.rates.lsr
    i_r=rec_problem.rates.lir
    upper_post_order=rec_problem.upper.post_order

    P=dict()
    P_TL=dict()
    P_max=dict()
    P_max_TL=dict()
    correction_ancestrale=dict()
    correction_ancestrale_TL=dict()
    correction_ancestrale_size=dict()

    if not P_transfer:
        for e in upper_post_order:
            if e.isRoot():
                correction_ancestrale_size[e]=1
            else:
                correction_ancestrale_size[e]=correction_ancestrale_size[e.parent]


    for e in upper_post_order:
        P_TL[e]=dict()
        P[e]=dict()


    if not mult_gene_match:
        for c in am_tree.leaves:
            P[c.match][c]=0
            P_TL[c.match][c]=0
    else:
        for c in am_tree.leaves:
            #c_match[c] si mult match c'est une liste
            for species_tmp in c.match:
                P[species_tmp][c]=(-1)*np.log(len(c.match))
                P_TL[species_tmp][c]=(-1)*np.log(len(c.match))

    clade_queue=list(am_tree.reverse_post_order)

    while len(clade_queue)>0:
        c=clade_queue.pop()
        upper_queue=list(upper_post_order)
        b_store=dict()
        #reachable_notl=notl_reachable(c,mult_gene_match)
        while len(upper_queue)>0:
            e=upper_queue.pop()
            if not mult_gene_match:
                bool_continue_mult_match=not (c.match==e)
                bool_continue_mult_match2=not(c.is_leaf()) or (c.match in e.leaves())
            else:
                bool_continue_mult_match=not (c.is_leaf() and e in c.match)
                bool_continue_mult_match2=not(c.is_leaf()) or (inter_list(c.match,e.leaves()))
            if bool_continue_mult_match:
                if bool_continue_mult_match2:
                    b_l=[]
                    if not e.isLeaf():
                        if c.is_leaf():
                            if mult_gene_match:
                                c_match_c=c.match
                            else:
                                c_match_c=[c.match]
                            if inter_list(c_match_c, e.right.leaves()):
                                b_l.append(E[e.left]+P_TL[e.right][c]+s_r)
                            if inter_list(c_match_c, e.left.leaves()):
                                b_l.append(E[e.right]+P_TL[e.left][c]+s_r)
                        else:
                            bl.append(E[e.left]+P_TL[e.right][c]+s_r)
                            bl.append(E[e.right]+P_TL[e.left][c]+s_r)

                    bchild=[]
                    for cL,cR in c.log_child_frequencies:
                        bclcr=[]
                        if not e.isLeaf():
                            bclcr.append(P[e.left][cL]+P[e.right][cR]+c.log_child_frequencies[(cL,cR)]+s_r)
                            bclcr.append(P[e.left][cR]+P[e.right][cL]+c.log_child_frequencies[(cL,cR)]+s_r)


                        #incomplete sorting event
                        if not i_r is None:
                            if not e.isLeaf():

                                b7=P[e.left][cL]+P[e][cR]+c.log_child_frequencies[(cL,cR)]+i_r
                                b8=P[e.left][cR]+P[e][cL]+c.log_child_frequencies[(cL,cR)]+i_r
                                b9=P[e.right][cL]+P[e][cR]+c.log_child_frequencies[(cL,cR)]+i_r
                                b10=P[e.right][cR]+P[e][cL]+c.log_child_frequencies[(cL,cR)]+i_r

                                bclcr.append(b7)
                                bclcr.append(b8)
                                bclcr.append(b9)
                                bclcr.append(b10)


                        if P_transfer:
                            for h in P_transfer[e].keys():
                                b3=c.log_child_frequencies[(cL,cR)]+t_r+P_transfer[e][h]+P[h][cL]+P[e][cR]
                                bclcr.append(b3)
                                b4=c.log_child_frequencies[(cL,cR)]+t_r+P_transfer[e][h]+P[h][cR]+P[e][cL]
                                bclcr.append(b4)
                        else:

                            b3=c.log_child_frequencies[(cL,cR)]+t_r-np.log(len(upper_post_order))+P_max[cL]+P[e][cR]
                            bclcr.append(b3)
                            b4=c.log_child_frequencies[(cL,cR)]+t_r-np.log(len(upper_post_order))+P_max[cR]+P[e][cL]
                            bclcr.append(b4)




                        #duplication

                        b6=d_r+P[e][cL]+P[e][cR]+c.log_child_frequencies[(cL,cR)]
                        bclcr.append(b6)

                        bchild.append(log_add_list(bclcr))
                    if len(bchild)>0:
                        b_l.append(max(bchild))

                    b=max(b_l)
                    P_TL[e][c]=b
        if c.is_leaf():
            transfer_list_e=[]
            q=list(upper_post_order)
            for e in upper_post_order:
                if c in P_TL[e]:
                    transfer_list_e.append(e)
                    if not c in P_max_TL:
                        P_max_TL[c]=P_TL[e][c]
                    else:
                        P_max_TL[c]=max(P_max_TL[c],P_TL[e][c])


        else:
            P_max_TL[c]=max([P_TL[e][c] for e in upper_post_order])
        upper_queue=list(upper_post_order)
        while len(upper_queue)>0:
            e=upper_queue.pop()
            if not mult_gene_match:
                bool_continue_mult_match=not (c.is_leaf() and c.match==e)
            else:
                bool_continue_mult_match=not (c.is_leaf() and e in c.match)
            if bool_continue_mult_match:
                if e in b_store:
                    b_l=[b_store[e]]
                    if P_transfer:
                        for h in P_transfer[e].keys():
                            if c in P_TL[h]:
                                b1=t_r+P_transfer[e][h]+P_TL[h][c]+E[e]
                                b_l.append(b1)
                    else:
                        b1=t_r-np.log(len(upper_post_order))+(P_max_TL[c])+E[e]
                        b_l.append(b1)

                    if not e.isLeaf():
                        if c in P_TL[e.right]:
                            if not P[e.right][c]==P_TL[e.right][c]:
                                b2=s_r+E[e.left]+P[e.right][c]
                                b_l.append(b2)
                        else:
                            b2=s_r+E[e.left]+P[e.right][c]
                            b_l.append(b2)
                        if c in P_TL[e.left]:
                            if not P[e.left][c]==P_TL[e.left][c]:
                                b3=s_r+E[e.right]+P[e.left][c]
                                b_l.append(b3)
                        else:
                            b3=s_r+E[e.right]+P[e.left][c]
                            b_l.append(b3)
                    b=max(b_l)
                    P[e][c]=b
                else:
                    b_l=[]
                    if P_transfer:
                        for h in P_transfer[e].keys():
                            if c in P_TL[h]:
                                b1=t_r+P_transfer[e][h]+P_TL[h][c]+E[e]
                                b_l.append(b1)
                    else:
                        b1=t_r-np.log(len(upper_post_order))+P_max_TL[c]+E[e]
                        b_l.append(b1)
                    if not e.isLeaf():
                        b2=s_r+E[e.left]+P[e.right][c]
                        b3=s_r+E[e.right]+P[e.left][c]
                        b_l.append(b2)
                        b_l.append(b3)
                    b=max(b_l)
                    P[e][c]=b
        P_max[c]=max([P[e][c] for e in upper_post_order])
    likelihood=max([P[k][am_tree] for k in P.keys()])



    likelihood=max([P[k][am_tree] for k in P.keys()])
    return P, P_TL,likelihood, correction_ancestrale_size


