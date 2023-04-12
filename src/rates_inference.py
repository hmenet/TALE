#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:16:35 2020

@author: hmenet
"""


import random as rd
import numpy as np

from reconciliation import reconciliation

####### to do ######

def count_events(l_event, min_rate=0.0001,non_used_events=[],geo=False):
    #geo=True
    d=dict()
    d["S"]=0
    d["L"]=0
    d["T"]=0
    d["D"]=0
    d["I"]=0
    total=0
    for l_e in l_event:
        for u in l_e:
            if u.name != "E":
                for c in u.name:
                    d[c]+=l_e[u]
                    total+=l_e[u]

    if geo:
        non_disp_ext_events=["S","D","I"]
        d_tmp=dict()
        for event in non_disp_ext_events:
            total+=sum([d[e] for e in non_disp_ext_events])-d[event]
            d_tmp[event]=sum([d[e] for e in non_disp_ext_events])
        for event in non_disp_ext_events:
            d[event]=d_tmp[event]
    added=0
    for c in d.keys():
        d[c]=d[c]/total
        if d[c]<=min_rate and not c in non_used_events:
            d[c]=min_rate
            added+=min_rate
    to_retrieve=added/len([d[c] for c in d.keys() if d[c]>2*min_rate])
    for c in d.keys():
        if d[c]>2*min_rate:
            d[c]-=to_retrieve

    return d

#####

### done #########


#version to search for maximum likelihood as in ALE or Corepa
def gene_rates_ml(rec):
    rec.rates_inference=True
    print(rec.rates.pp())
    i_steps=0
    while i_steps < rec.n_steps :
        rec_sol=reconciliation(rec)
        if rec.rates.ir ==0 :
            non_used_events=["I"]
        else:
            non_used_events=[]
        rates=count_events(rec_sol.event_list_by_fam,non_used_events=non_used_events,geo=rec.geo)
        rec.rates.dr=rates["D"]
        rec.rates.tr=rates["T"]
        rec.rates.lr=rates["L"]
        rec.rates.ir=rates["I"]
        rec.rates.reinit()
        i_steps+=1
        print("rates estimation step ", i_steps, "/",rec.n_steps, " likelihood: ", rec_sol.log_likelihood, "rates :", rec.rates.pp())
    rec.rates_inference=False

