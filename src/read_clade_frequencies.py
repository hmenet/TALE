#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:19:06 2020

@author: hmenet
"""







import arbre

#il faut quelque chose de hashable dans un dictionnaire donc on
#prend les tuples obtenu a partir les listes triés de feuilles communes
#et on va utiliser un dictionnaire qui prend en entrée ces tuples de feuilles

#compter les bipartitions et tripartitions

#l'idée, si j'ai bien compris, c'est 3 choses
#1 dans l'arbre de gene on ne regarde plus les sommets mais les clades (ensembles de feuilles) (puisqu'il n'y a plus de sommet car plus d'arbre, mais un ensemble d'arbre)
#2 un clade n'a pas 2 fils comme un sommet d'un arbre mais plein de couples possibles qui correspondent à ses bipartitions
#3 on associe à chaque clade sa fréquence dans l'échantillon pour le pondérer

#il faut donc connaitre les fréquences des partitions de clades
#un clade : liste de feuilles, rangé par indice ? ou dans l'ordre donné par post order ?
#liste de clade

def sorted_common_leaves_from_tree(tree):
    l=tree.leaves_unrooted()
    l_leaves=[i.name for i in l]#maintenant les feuilles de la liste ont des indices qui correspondent bien entre arbre
    l_leaves.sort()
    return l_leaves

def sorted_common_leaves_from_leaves(leaves):
    l_leaves=[i.name for i in leaves]#maintenant les feuilles de la liste ont des indices qui correspondent bien entre arbre
    l_leaves.sort()
    return l_leaves


#les clades sont des indices entiers qui commence à 0

#try to add a clade, if not present
#tree_clade : a tree, the clade will be its leaves
#a clade is just an index, in clade_elements[clade] : sorted list of leaves corresponding to the clade
#and the leaves are the leaves names
def clade_from_tree(tree_clade, clade_elements, clade_keys):
    n_clade=len(clade_elements)
    leaves_clade=tuple(sorted_common_leaves_from_tree(tree_clade))
    clade=-1
    if leaves_clade in clade_keys:
        clade=clade_keys[leaves_clade]
    else:
        clade=n_clade
        #on met à jour clade_elements avec les elements de ce clade
        clade_elements[clade]=leaves_clade
        clade_keys[leaves_clade]=clade
    return clade

def clade_from_leaves(leaves_clade, clade_elements, clade_keys):
    n_clade=len(clade_elements)
    leaves_clade=tuple(sorted_common_leaves_from_leaves(leaves_clade))
    clade=-1
    if leaves_clade in clade_keys:
        clade=clade_keys[leaves_clade]
    else:
        clade=n_clade
        #on met à jour clade_elements avec les elements de ce clade
        clade_elements[clade]=leaves_clade
        clade_keys[leaves_clade]=clade
    return clade



#a la fin on veut travailler seulement avec des clades
# et être simplement capable de retrouver les feuilles correspondant à chaque clade

#clade_frequencies[u][(v,w)]=0# bipartition v w sachant bipartition u u barre
#l_tree : liste d'arbre, dont on veut obtenir les fréquences des clades
#il faut faire attention à prendre les trois racines possibles pour chaque bipartition


def partition_counter_one_tree(tree, clade_elements, clade_keys, bipartition_number, tripartition_number):
    #if rooted, we unroot
    if tree.right2==None:
        arbre.from_rooted_to_un(tree)
    tree_post_order = tree.post_order_traversal_unrooted()
    all_leaves=tree.leaves_unrooted()
    clade1=clade_from_leaves(all_leaves, clade_elements, clade_keys)#on definit le premier clade comme celui de toutes les feuilles
    for e in tree_post_order:
        if not e.isRoot():
            #on va couper l'arete qui relie e a son pere, creant un bipartition
            c1=e.leaves_unrooted()#la moitié de la bipartition
            c2=[leaf for leaf in all_leaves if not leaf in c1]
            clade1=clade_from_leaves(c1, clade_elements, clade_keys)
            clade2=clade_from_leaves(c2, clade_elements, clade_keys)
            bipartition=(min(clade1, clade2), max(clade1, clade2))
            tmp_compte_a_prendre=1
            #en unrooted la racine est un vrai noeud
            #if not e.isRoot() and e.parent==tree:
            #    tmp_compte_a_prendre=0.5 #la racine n'est pas un vrai noeud mais une branche
            if bipartition in bipartition_number:
                bipartition_number[bipartition]+=tmp_compte_a_prendre
            else:
                bipartition_number[bipartition]=tmp_compte_a_prendre

        #maintenant on compte la tripartition generé en supprimant le noeud considéré
        #ici pas besoin de faire ni les feuilles, ni la racine puisque ce n'est pas un vrai noeud
        if not e.isLeaf():
            #on a deja calculé un des clades de la tripartition, c2
            #on calcule maintenant les deux autres
            if not e.isRoot():
                child_clade1=clade_from_tree(e.left, clade_elements, clade_keys)
                child_clade2=clade_from_tree(e.right, clade_elements, clade_keys)
                child_clade3=clade2
            else:
                #si c'est la racine on a encore calculé aucun des clades
                child_clade1=clade_from_tree(e.left, clade_elements, clade_keys)
                child_clade2=clade_from_tree(e.right, clade_elements, clade_keys)
                child_clade3=clade_from_tree(e.right2, clade_elements, clade_keys)
            #peut etre que dans le preprocess il vaudra mieux parcourir d'abord les bipartitions ?
            l_tmp=[child_clade1, child_clade2, child_clade3]
            l_tmp.sort()
            c1_tmp, c2_tmp,c3_tmp=l_tmp
            tripartition=c1_tmp, c2_tmp,c3_tmp
            if tripartition in tripartition_number:
                tripartition_number[tripartition]+=1
            else:
                tripartition_number[tripartition]=1

#renvoie le clade correspondant à l'union des deux clades données en entrées
def union_clade(c1,c2,clade_elements, clade_keys):
    l=list(clade_elements[c1])
    l+=list(clade_elements[c2])
    l.sort()
    c=clade_keys[tuple(l)]
    return c



#on part des feuilles
#il suffit de sort par la taille
def construct_clade_post_order(clade_frequencies, clade_elements):
    def order_clade(clade):
        return(len(clade_elements[clade]))
    l=[i for i in clade_elements.keys()]
    l.sort(key=order_clade, reverse=True)
    #ordre donné par la taille des clades est bien un ordre qui empeche les fils d'apparaitre avant leur parent
    #cependant on ne supprime pas les elements qui ne descendent pas (et donc ne mene pas) à la racine
    first_clade=0#on va chercher le plus grand des clades, en supposant qu'il y en a bien un qui contienne tout les autres
    for u in clade_elements:
        if len(clade_elements[u])>len(clade_elements[first_clade]):
            first_clade=u
    e=set()
    to_see=[first_clade]
    while len(to_see)>0:
        c=to_see.pop()
        e.add(c)
        for (cL,cR) in clade_frequencies[c]:
            if not cL in e:
                to_see.append(cL)
            if not cR in e:
                to_see.append(cR)
    #maintenant e contient tout les elements atteignable depuis la racine
    new_l=[]
    for c in l:
        if c in e:
            new_l.append(c)
    return new_l



#on fait attention a ce que les arbres ne sont pas enracinés
def compute_clade_frequencies(l_tree):
    bipartition_number=dict()
    tripartition_number=dict()
    clade_elements=dict()
    clade_keys=dict()
    #union_clade=dict()
    kcmpt=0
    for tree in l_tree:
        kcmpt+=1
        partition_counter_one_tree(tree, clade_elements, clade_keys,bipartition_number, tripartition_number)
        if kcmpt % 100 ==0:
            print(kcmpt, "/", len(l_tree))
    clade_frequencies=dict()
    n_clade=len(clade_elements)
    for c in range(n_clade):
        clade_frequencies[c]=dict()
    for c1,c2,c3 in tripartition_number:
        #pour faire mieux (plus rapide) il faudrait garder en mémoire les relations de parenté entre les clades
        n_tripartition=tripartition_number[(c1,c2,c3)]
        #a1 a2 sachant a, et ab représente a barre
        for tmp1,tmp2,tmp3 in [(c2,c3,c1),(c1,c3,c2),(c1,c2,c3)]:
            #a,a1,a2,ab = union_clade[(min(tmp1,tmp2), max(tmp1,tmp2))],tmp1,tmp2,tmp3
            #print("a")

            a,a1,a2,ab = union_clade(tmp1,tmp2,clade_elements, clade_keys),tmp1,tmp2,tmp3
            b_continue=True
            for u in [a,a1,a2]:
                if len(clade_elements[u])==0:
                    b_continue=False
            if b_continue:
                clade_frequencies[a][(min(a1,a2),max(a1,a2))]=n_tripartition/bipartition_number[(min(a,ab),max(a,ab))]
    #pour la fréquence de bipartition du premier clade il n'y a pas de sachant que, pas de tripartition
    s=sum([bipartition_number[u] for u in bipartition_number.keys()])
    for c1,c2 in bipartition_number:
        clade_frequencies[0][(min(c1,c2),max(c1,c2))]=bipartition_number[(min(c1,c2),max(c1,c2))]/s
    clade_post_order=construct_clade_post_order(clade_frequencies, clade_elements)
    return clade_post_order, clade_frequencies, clade_elements, clade_keys


def compute_clade_frequencies_multiple_families(l_tree_by_family):
    clades_data_list=[]
    for l_tree in l_tree_by_family:
        clade_data=compute_clade_frequencies(l_tree)
        #les noms des arbres sont "nom"+str(id_counter), donc pour avoir le nom, on enlève le 0 du nom du premier arbre
        clades_data_list.append(clade_data)
    return clades_data_list

def from_rooted_tree_to_clades(tree):
    clade_elements=dict()
    clade_keys=dict()
    clade_frequencies=dict()
    clade_to_tree=dict()
    c=0
    for e in tree.post_order_traversal():
        l_tmp=[u.name for u in e.leaves()]
        l_tmp.sort()
        leaves = tuple(l_tmp)
        clade_elements[c]=leaves
        clade_keys[leaves]=c
        clade_to_tree[e]=c
        if not e.isRoot():
            p=e.parent
            pl_tmp=[u.name for u in p.leaves()]
            pl_tmp.sort()
            pleaves=tuple(pl_tmp)
            pc=clade_keys[pleaves]
            cl_tmp=[u.name for u in p.left.leaves()]
            cr_tmp=[u.name for u in p.right.leaves()]
            cl_tmp.sort()
            cr_tmp.sort()
            clleaves=tuple(cl_tmp)
            crleaves=tuple(cr_tmp)
            if clleaves in clade_keys and crleaves in clade_keys:
                cl=clade_keys[clleaves]
                cr=clade_keys[crleaves]
                clade_frequencies[pc]=dict()
                clade_frequencies[pc][(min(cl,cr), max(cl,cr))]=1
        if e.isLeaf():
            clade_frequencies[c]=dict()
        c+=1
    clade_post_order=construct_clade_post_order(clade_frequencies, clade_elements)
    return (clade_post_order, clade_frequencies, clade_elements, clade_keys), clade_to_tree

#from rooted tree to clades data list, with no amalgamation or root uncertainty, for the intermediate level of 3 level reconciliation
def compute_clade_frequencies_rooted_tree(l_tree_by_family):
    clade_to_tree_list=[]
    clades_data_list=[]
    for tree in l_tree_by_family:
        clade_data, clade_to_tree=from_rooted_tree_to_clades(tree)
        clades_data_list.append(clade_data)
        clade_to_tree_list.append(clade_to_tree)
    return clades_data_list, clade_to_tree_list



def clade_to_name(clades_data_list, tree):
    d=dict()
    d1=tree.name_to_tree()
    clade_post_order, clade_frequencies, clade_elements, clade_keys = clades_data_list
    for c in clade_elements.keys():
        l_leaves=clade_elements[c]
        a_leaf=l_leaves[0]
        d[c]=(d1[a_leaf].n_leaves_ancestor(len(l_leaves))).name
    return d

def clade_to_name_by_fam(clades_data_list_by_fam, tree_list):
    d_list=[]
    for i_clades in range(len(tree_list)):
        d_list.append(clade_to_name(clades_data_list_by_fam[i_clades], tree_list[i_clades]))
    return d_list





