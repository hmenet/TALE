#to read Sagephy simulated files
#for comparing 3 level reconciliation

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:33:46 2020

@author: hmenet
"""
from os import walk, path, makedirs


import random as rd
import numpy as np
import matplotlib.pyplot as plt

import arbre
from read_sagephy_aux import compute_clade_frequencies_multiple_families


### read sagephy file  #####


def read_from_to(s, start_index, stop_caracters_list):
    c=s[start_index]
    word=""
    while not c in stop_caracters_list:
        word+=c
        start_index+=1
        c=s[start_index]
    start_index+=1
    return word, start_index


def read_tree_string_sagephy(s, more_output=False):
    nodes_info=dict()
    tree=arbre.Tree()
    tree.id=1
    tree.root=tree
    node=tree
    k=0
    while k < len(s):
        c=s[k]
        if c=="(":
            node.left=arbre.Tree()
            node.left.parent=node
            node.left.root=node.root
            node.left.id=node.id*2
            node=node.left
        elif c==",":
            node=node.parent

            if node.right== None:
                node.right=arbre.Tree()
                node.right.parent=node
                node.right.root=node.root
                node.right.id=node.id*2+1
                node=node.right
            else:
                if node.right2==None:
                    node.right2=arbre.Tree()
                    node.right2.parent=node
                    node.right2.root=node.root
                    node.right2.id=node.id*2+1.5#noeud qui apparait quand l'arbre n'est pas raciné
                    node=node.right2
                else:
                    print("Wasn't expecting multifurcating tree")
        elif c==")":
            node=node.parent
        elif c==";":
            if more_output:
                return tree, nodes_info
            else:
                return tree
        else:
            name=""
            while not c in ["(",")",",",":", ";", "[","]"]:
                name+=c
                k+=1
                c=s[k]
            if c==":":
                while not c in ["(",")",",",";","[","]"]:
                    k+=1
                    c=s[k]
            if c=="[":
                words,k = read_from_to(s,k+1,["]"])
                l_word=words.split(sep=" ")
                l_word2=[]
                for u in l_word:
                    l_word2.append(u.split(sep="="))
                if not node in nodes_info:
                    #voir ce qu'on fait avec les infos dans l'autre cas
                    nodes_info[node]=l_word2


            if len(name)>0:
                node.name=name
            k-=1
        k+=1

#return a tree
def read_tree_sagephy(file, more_output=False):
    f=open(file, "r")
    s=f.read()
    return read_tree_string_sagephy(s, more_output=more_output)


########


#rd.seed(a=2)


def construct_file_list(directory):
    liste_fichiers = []
    for (repertoire, sousRepertoires, fichiers) in walk(directory):
        liste_fichiers.extend(fichiers)
    return liste_fichiers

#selon si fonction considère liste d'arbre ou arbre (pour amalgamation)
def construct_tree_list(directory, amalgamation=False):
    nodes_info_list=[]
    tree_list=[]
    tree_file_list=construct_file_list(directory)
    for tree_file in tree_file_list:
        #normalement il faudra choisir une autre fonction read tree, ou passer amalgamation à celle ci
        if amalgamation:
            tree,nodes_info=read_tree_sagephy(directory+tree_file, more_output=True)
            tree_list2=[tree]
            nodes_info_list.append(nodes_info)
            tree_list.append(tree_list2)
            id_counter=0
            for tree in tree_list2:
                tree.tree_name=tree_file+str(id_counter)
                id_counter+=1
        else:
            tree,nodes_info=read_tree_sagephy(directory+tree_file, more_output=True)
            nodes_info_list.append(nodes_info)
            tree.tree_name=tree_file
            tree_list.append(tree)
    return tree_list, tree_file_list, nodes_info_list



#return a  dict matching inner leaves to host leaves
# selon le format dans sagephy

#voir selon ce qu'on veut en faire apres
def construct_leaves_matching(matching_file, lower=False):
    d=dict()
    f=open(matching_file, "r")
    s=f.read()
    s1=s.split(sep="\n")
    for  u in s1:
        if len(u)>0:
            s2=u.split(sep="\t")
            if lower:
                d[s2[0]]=[s2[1],s2[2]]
            else:
                d[s2[0]]=s2[1]
    return d

#file_list contain the name of the files, as explored for the trees, so that the gene list and matching lists correspond to one another
def construct_leaves_matching_dir(matching_directory, file_list, lower=False):
    new_file_list=[]#some minor changes between leaf mapping and gene trees names
    for s in file_list:
        new_file_list.append(s[:-4]+"leafmap")
    leaf_matching_list=[]
    for i_file in range(len(new_file_list)):
        #pour correspondre aux noms donnés par sagephy
        d=construct_leaves_matching(matching_directory+new_file_list[i_file], lower=lower)
        leaf_matching_list.append(d)
    return leaf_matching_list

#à chaque nom de noeud renvoie le noeud correspondant
def construct_name_to_leaves(tree):
    d=dict()
    for u in tree.leaves():
        d[u.name]=u
    return d

def construct_name_to_clade(clade_keys):
    d=dict()
    for u in clade_keys:
        if len(u)==1:
            d[u[0]]=clade_keys[u]
    return d

#à chaque nom d'arbre renvoie l'arbre corrspondant
def construct_name_to_tree(tree_list):
    d=dict()
    for t in tree_list:
        d[t.tree_name]=t
    return d

def name_matching_to_tree_matching(name_to_node_host, name_to_node_inner,leaf_matching,lower=False, tree=False):
    if lower:
        d=dict()
    for inner_node_name in leaf_matching:
        if inner_node_name in name_to_node_inner : #sinon la feuille a été pruned de l'arbre
            inner_node=name_to_node_inner[inner_node_name]
            if lower:
                host_node_name,host_tree_name=leaf_matching[inner_node_name]
                host_node=name_to_node_host[host_tree_name][host_node_name]
                if tree:
                    inner_node.match=host_node
                else:
                    d[inner_node]=host_node
            else:
                host_node_name=leaf_matching[inner_node_name]
                host_node=name_to_node_host[host_node_name]
                inner_node.match=host_node
    if lower:
        return d

#to prune
def pre_match(host_list, symbiont_list, gene_list, upper_leaf_matching_list, lower_leaf_matching_list):
    name_to_node_host=dict()
    host=host_list[0]
    name_to_node_host=construct_name_to_leaves(host)
    for i_symbiont in range(len(symbiont_list)):
        symbiont=symbiont_list[i_symbiont]
        upper_leaf_matching=upper_leaf_matching_list[i_symbiont]
        name_to_node_inner=construct_name_to_leaves(symbiont)
        name_matching_to_tree_matching(name_to_node_host, name_to_node_inner, upper_leaf_matching, lower=False, tree=True)
    name_to_node_symbionts=dict()
    for symbiont in symbiont_list:
        name_to_node_symbionts[symbiont.tree_name]=construct_name_to_leaves(symbiont)
    for i_gene in range(len(gene_list)):
        gene_l=gene_list[i_gene]
        for gene in gene_l:
            lower_leaf_matching=lower_leaf_matching_list[i_gene]
            name_to_node_inner=construct_name_to_leaves(gene)
            name_matching_to_tree_matching(name_to_node_symbionts, name_to_node_inner, lower_leaf_matching, lower=True, tree=True)


def all_name_matching_to_tree_matching(host_list, symbiont_list, clades_data_list, upper_leaf_matching_list, lower_leaf_matching_list):
    name_to_node_host=dict()
    host=host_list[0]
    name_to_node_host=construct_name_to_leaves(host)
    for i_symbiont in range(len(symbiont_list)):
        symbiont=symbiont_list[i_symbiont]
        upper_leaf_matching=upper_leaf_matching_list[i_symbiont]
        name_to_node_inner=construct_name_to_leaves(symbiont)
        name_matching_to_tree_matching(name_to_node_host, name_to_node_inner, upper_leaf_matching, lower=False)
    lower_leaf_matching_nodes_list=[]
    name_to_node_symbionts=dict()
    for symbiont in symbiont_list:
        name_to_node_symbionts[symbiont.tree_name]=construct_name_to_leaves(symbiont)
    for clade_id in range(len(clades_data_list)):
        clade_keys=clades_data_list[clade_id][3]
        lower_leaf_matching=lower_leaf_matching_list[clade_id]
        name_to_node_inner=construct_name_to_clade(clade_keys)
        d=name_matching_to_tree_matching(name_to_node_symbionts, name_to_node_inner, lower_leaf_matching, lower=True)
        lower_leaf_matching_nodes_list.append(d)
    return lower_leaf_matching_nodes_list
















#################################################
### pour faire des tests sur les transferts######
#################################################

def tree_list_post_order_traversal(tree_list):
    list_post_order=[]
    for tree in tree_list:
        list_post_order+=tree.post_order_traversal()
    return list_post_order



#put to 0.5 the time of all unsampled and lost leaves
def nodes_lost(nodes_info):
    for node in nodes_info:
        for u in node:
            v=node[u]
            if v[0]=="VERTEXTYPE" and v[1] in ["UnsampledLeaf","Loss"]:
                node.time=0.5
            if v[0]=="VERTEXTYPE" and v[1]=="Leaf":
                node.time=1


def in_tree_descendant(tree):
    return len([u for u in tree.leaves() if u.time==1])>0

def first_in_tree_ancestor(tree):
    if in_tree_descendant(tree):
        return tree
    else:
        return first_in_tree_ancestor(tree.parent)

#if a node is collapsed with its parent, return the node in the tree, where it is collapsed
# a faire avec l'arbre pruned
def renamed_after_pruned(name, tree):
    for u in tree.liste():
        if name in u.pruned_name_list:
            return u

def name_to_tree(name, tree):
    for u in tree.liste():
        if u.name==name:
            return u



def transfer_list(nodes_info, tree_and_id_to_node_host, pruned_tree_from_name):
    l_t=[]
    for node in nodes_info:
        if in_tree_descendant(node):#on ajoute le transfert seulement si un descendant du gène transféré a atteint les feuilles
            if ['Transfer'] in nodes_info[node]:
                for u in nodes_info[node]:
                    if u[0]=="FROMTOLINEAGE":
                        transfer_info=u[1][1:-1]
                        s=transfer_info.split(sep=";")
                        s0=s[0]
                        s1=s[1]
                        gvr_node, rcv_node = s0.split(sep=",")
                        gvr_tree, rcv_tree = s1.split(sep=",")

                        gvr=tree_and_id_to_node_host[gvr_tree][gvr_node]
                        rcv=tree_and_id_to_node_host[rcv_tree][rcv_node]

                        #on transforme directement en transfer pruned

                        pruned_gvr_tree=pruned_tree_from_name[gvr_tree]
                        pruned_rcv_tree=pruned_tree_from_name[rcv_tree]
                        if in_tree_descendant(rcv):#si pas de descendant du receveur, on ajoute pas
                            gvr_unpruned_ancestor=first_in_tree_ancestor(gvr)
                            pruned_gvr=renamed_after_pruned(gvr_unpruned_ancestor.name, pruned_gvr_tree)
                            pruned_rcv=renamed_after_pruned(rcv.name, pruned_rcv_tree)
                            l_t.append([pruned_gvr, pruned_rcv])
    return l_t



def construct_tree_and_id_to_node_host(nodes_info_list, symbiont_list):
    d=dict()
    for i_symbiont in range(len(nodes_info_list)):
        symbiont=symbiont_list[i_symbiont]
        d[symbiont.tree_name]=dict()
        nodes_info=nodes_info_list[i_symbiont]
        for node in nodes_info:
            for u in nodes_info[node]:
                if u[0]=="ID":
                    d[symbiont.tree_name][u[1]]=node
    return d


def from_unpruned_to_pruned_with_info(host,host_nodes_info_list,symbiont_list, symbiont_nodes_info_list, gene_list, gene_nodes_info_list, observed_proba=1):

    nodes_lost(host_nodes_info_list)

    for leaf in host.leaves():
        x=rd.random()
        if x > observed_proba:
            leaf.time=0.5
        else :
            leaf.time=1

    #en fait sagephy renvoit seulement des arbres pruned
    #set the time to 0 for all lost or unsampled leaves in all symbiont tree
    #nodes_lost(symbiont_nodes_info)
    #nodes_lost(gene_nodes_info_list)


    for symbiont in symbiont_list:
        for leaf in symbiont.leaves():
            leaf.time=leaf.match.time
    for gene_l in gene_list:
        for gene in gene_l:
            for leaf in gene.leaves():
                leaf.time=leaf.match.time

    pruned_tree_from_name=dict()
    symbiont_pruned_list=[]

    for i_symbiont in range(len(symbiont_list)):
        symbiont=symbiont_list[i_symbiont]
        for u in symbiont.liste():
            u.pruned_name_list=[u.name]
        symbiont_pruned=symbiont.copy()
        symbiont_pruned.prune_tree()
        pruned_tree_from_name[symbiont_pruned.tree_name]=symbiont_pruned
        symbiont_pruned_list.append(symbiont_pruned)

    #construction de la liste des tupper_leaf_matching_directoryransfers
    tree_and_id_to_node_host=construct_tree_and_id_to_node_host(symbiont_nodes_info_list, symbiont_list)
    l_t_list=[]
    for gene_nodes_info in gene_nodes_info_list:
        l_t=transfer_list(gene_nodes_info, tree_and_id_to_node_host, pruned_tree_from_name)
        l_t_list.append(l_t)

    for gene_l in gene_list:
        for gene in gene_l:
            if len([u for u in gene.leaves() if u.time == 1]) < 5 :
                raise ValueError("pas assez de feuilles")
            print(len([u for u in gene.leaves() if u.time == 1]))
            gene.prune_tree()
    host.prune_tree()
    return symbiont_pruned_list, l_t_list


#on prend entrée les arbres non coupés, sauf l'arbre le plus bas (gene), qu'on prend deja pruned

#on les coupe en gardant les infos



#################################


def save_trees(tree_list, tree_file_list,output_path):
    if not path.isdir(output_path):
        makedirs(output_path)
    for i_tree in range(len(tree_list)):
        tree_file=tree_file_list[i_tree]
        arbre.save_tree(tree_list[i_tree],output_path+tree_file+".nwk")

def save_transfer(l_t,output_file):
    s=""
    for e in l_t:
        s+=e[0].name+"\t"+e[1].name+"\n"
    f=open(output_file, "w")
    f.write(s)
    f.close()


def save_transfer_list(l_t_list,tree_file_list,output_path):
    if not path.isdir(output_path):
        makedirs(output_path)
    for i_tree in range(len(l_t_list)):
        save_transfer(l_t_list[i_tree],output_path+tree_file_list[i_tree])


def save_matching(tree, output_file):
    s=""
    for u in tree.leaves():
        s+=u.name+"\t"+u.match.name+"\t"+u.match.root.tree_name+".nwk"+"\n"
    f=open(output_file, "w")
    f.write(s)
    f.close()

def save_matchings(tree_list,tree_file_list,output_path):
    if not path.isdir(output_path):
        makedirs(output_path)
    for i_tree in range(len(tree_list)):
        save_matching(tree_list[i_tree],output_path+tree_file_list[i_tree]+".leafmap")


#sagephy selon si on veut l'info des simulations
def sagephy_to_3_level_input(data_dir,host_file, symbiont_directory, gene_directory, upper_leaf_matching_directory, lower_leaf_matching_directory, output_path, observed_proba=1):



    host, host_nodes_info=read_tree_sagephy(host_file,more_output=True)

    host_file=host_file[host_file.rfind("/")+1:]
    host.tree_name=host_file

    symbiont_list, symbiont_file_list, symbiont_nodes_info_list=construct_tree_list(symbiont_directory, amalgamation=False)
    gene_list, gene_file_list, gene_nodes_info_list=construct_tree_list(gene_directory, amalgamation=True)

    host_list=[host]

    ### leaf matching #####

    construct_file_list(upper_leaf_matching_directory)
    upper_leaf_matching_list=construct_leaves_matching_dir(upper_leaf_matching_directory,symbiont_file_list, lower=False)
    lower_leaf_matching_list=construct_leaves_matching_dir(lower_leaf_matching_directory,gene_file_list, lower=True)




    ## prune ####
    #on a besoin d'un match dès maintenant pour couper les feuilles des arbres du dessous si leur match sont coupés
    #et d'un autre plus loin, puisque on travaille ensuite sur une copie des arbres, et pour les avoir tous pruned
    print([len(symbiont.liste()) for symbiont in symbiont_list])
    pre_match(host_list,symbiont_list,gene_list,upper_leaf_matching_list, lower_leaf_matching_list)
    while True:
        try:
            symbiont_list, l_t_list=from_unpruned_to_pruned_with_info(host, [host_nodes_info], symbiont_list, symbiont_nodes_info_list, gene_list, gene_nodes_info_list, observed_proba=observed_proba)
            break
        except ValueError:
            for u in gene_list:
                for v in u :
                    print(len(v.leaves()))
                    if len(v.leaves()) < 6 :
                        raise ValueError("Pas assez de feuille")
            print(data_dir)
    print([len(symbiont.liste()) for symbiont in symbiont_list])

    ### genes clade frequencies computation #####


    if not path.isdir(output_path):
        makedirs(output_path)

    new_gene_list=[]
    for u in gene_list:
        new_gene_list.append(u[0])

    save_trees([host],[host_file],output_path+"/species/")
    save_trees(symbiont_list,symbiont_file_list,output_path+"/symbiont/")

    save_transfer_list(l_t_list,gene_file_list,output_path+"/transfer_list/")

    save_matchings(new_gene_list,gene_file_list,output_path+"/lower_matching/")

    clades_data_list= compute_clade_frequencies_multiple_families(gene_list)


    c_match_list=all_name_matching_to_tree_matching(host_list,symbiont_list,clades_data_list,upper_leaf_matching_list, lower_leaf_matching_list)

    save_matchings(symbiont_list,symbiont_file_list,output_path+"/upper_matching/")

    save_trees(new_gene_list,gene_file_list,output_path+"/genes/")


    return host_list, symbiont_list, clades_data_list, upper_leaf_matching_list, c_match_list, l_t_list


"""

sim_list_list=[]

sim_list_list=[["sim_290621_med_"+j+"_coevol"+str(i)+"/" for i in range(1,11)] for j in ["02","05","08","10"]]

#sim_list_list.append(["sim_290621_highmed"+str(i)+"/" for i in range(1,6)])

#sim_list_list.append(["sim_290621_highlow"+str(i)+"/" for i in range(1,6)])
#sim_list_list.append(["sim_290621_high"+str(i)+"/" for i in range(1,6)])
#sim_list_list.append(["sim_290621_med"+str(i)+"/" for i in range(1,6)])
#sim_list_list.append(["sim_290621_low"+str(i)+"/" for i in range(1,6)])
sim_list_list=[]
for i_rate in ["1.0"]:#["0.0","0.2","0.4","0.6","0.8","1.0"]:
    sim_list_list.append(["med_coevol_"+i_rate+"_"+str(i)+"/" for i in range(1,50)])

for sim_list in sim_list_list:
    for sim in sim_list:
        if not sim in ["med_coevol_0.2_1/","med_coevol_0.2_5/","med_coevol_0.2_44/","med_coevol_0.4_3/","med_coevol_0.4_46/","med_coevol_0.4_31/","med_coevol_0.4_48/","med_coevol_1.0_15/","med_coevol_0.6_7/","med_coevol_0.6_9/","med_coevol_1.0_23/"]:
            output_path="/home/hmenet/Documents/Stage_M2/These/script/simulation_101121/three_level_input/coevol"
            #print(output_path+"/" +sim[:-1])
            if not path.exists(output_path+"/" +sim[:-1]):
                nsim=""
                if sim[-3]=="_":
                    nsim=int(sim[-2])
                else:
                    nsim=int(sim[-3:-1])
                #data_dir="/home/hmenet/Documents/Stage_M2/These/script/simulation/"+sim
                data_dir="/home/hmenet/Documents/Stage_M2/These/script/simulation_101121/sagephy_output/coevol/"+sim
                host_file="species.pruned.tree"
                symbiont_dir="symbiont_tree/"
                gene_dir="gene_tree/"
                upper_match="upper_matching/"
                lower_match="lower_matching/"
                file_links=[host_file, symbiont_dir,gene_dir, upper_match, lower_match]
                for i in range(len(file_links)):
                    file_links[i]=data_dir+file_links[i]
                host_file, symbiont_dir,gene_dir, upper_match, lower_match=file_links


                #output_path="/home/hmenet/Documents/Stage_M2/These/script/simulation/test_270721_sagephy"
                #if not path.isdir(output_path):
                #    makedirs(output_path)


                try:
                    host_list, symbiont_list, clades_data_list, upper_leaf_matching, c_match_list,l_t_list= sagephy_to_3_level_input(data_dir,host_file,symbiont_dir, gene_dir, upper_match, lower_match,output_path+"/" +sim[:-1],observed_proba=0.2)#0.15#0.08)
                except ValueError:
                    print(sim)

                print("\n lecture données ",sim," terminée")

"""

sim_list_list=[]
for i_rate in ["1.0"]:
    sim_list_list.append(["rep_tr"+str(i)+"/" for i in range(1,2)])

for sim_list in sim_list_list:
    for sim in sim_list:
        if not sim in []:
            output_path="/home/hmenet/Documents/Stage_M2/These/script/sim_pap2/three_level_input"
            #print(output_path+"/" +sim[:-1])
            if not path.exists(output_path+"/" +sim[:-1]):
                #data_dir="/home/hmenet/Documents/Stage_M2/These/script/simulation/"+sim
                data_dir="/home/hmenet/Documents/Stage_M2/These/script/sim_pap2/"+sim
                host_file="species.pruned.tree"
                symbiont_dir="symbiont_tree/"
                gene_dir="gene_tree/"
                upper_match="upper_matching/"
                lower_match="lower_matching/"
                file_links=[host_file, symbiont_dir,gene_dir, upper_match, lower_match]
                for i in range(len(file_links)):
                    file_links[i]=data_dir+file_links[i]
                host_file, symbiont_dir,gene_dir, upper_match, lower_match=file_links


                #output_path="/home/hmenet/Documents/Stage_M2/These/script/simulation/test_270721_sagephy"
                #if not path.isdir(output_path):
                #    makedirs(output_path)


                try:
                    host_list, symbiont_list, clades_data_list, upper_leaf_matching, c_match_list,l_t_list= sagephy_to_3_level_input(data_dir,host_file,symbiont_dir, gene_dir, upper_match, lower_match,output_path+"/" +sim[:-1],observed_proba=0.2)#0.15#0.08)
                except ValueError:
                    print(sim)

                print("\n lecture données ",sim," terminée")
