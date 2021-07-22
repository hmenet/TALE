#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:16:35 2020

@author: hmenet
"""


import random as rd
import numpy as np

from reconciliation import reconciliation


def count_events(l_event, min_rate=0.0001):
    d=dict()
    d["S"]=0
    d["L"]=0
    d["T"]=0
    d["D"]=0
    total=0
    for l_e in l_event:
        for u in l_e:
            for c in u[0]:
                d[c]+=l_e[u]
                total+=l_e[u]
    for c in d.keys():
        d[c]=d[c]/total
        if d[c]==0:
            d[c]=min_rate
    return d


#version pour chercher maximum likelihood comme dans ALE ou Corepa
def gene_rates_ml(parasite_post_order,clades_data_list,c_match_list, n_steps, init_rates_g=[0.01,0.01,0.01], n_rec_sample=100, multi_process=False, multi_process_family=False, upper_input=None):

    rates_g=dict()
    rates_g["D"]=init_rates_g[0]
    rates_g["L"]=init_rates_g[1]
    rates_g["T"]=init_rates_g[2]

    i_steps=0
    while i_steps < n_steps :
        output_table=reconciliation(parasite_post_order, clades_data_list, c_match_list, rates_g,sample=True, n_sample=n_rec_sample, multi_process=multi_process, multi_process_family=multi_process_family,upper_input=upper_input)
        likelihood=output_table[0]
        l_event_gene=output_table[1]
        rates_g=count_events(l_event_gene)
        i_steps+=1
        print("rates estimation step ", i_steps, "/",n_steps, " likelihood: ", likelihood, "rates :", rates_g)
    return rates_g
