### Prepare the input, from repo containing trees in newick (or list of tree in newick) to tree class and/or clades amalgamation format

from os import walk
import numpy as np
from read_tree import read_tree, read_mult_tree
from read_clade_frequencies import compute_clade_frequencies_multiple_families, compute_clade_frequencies_rooted_tree
from rec_classes import Tree_list


### read the input datas  #####


## input trees
def construct_file_list(directory):
    liste_fichiers = []
    for (repertoire, sousRepertoires, fichiers) in walk(directory):
        liste_fichiers.extend(fichiers)
    return liste_fichiers

#selon si fonction considère liste d'arbre ou arbre (pour amalgamation)
def construct_tree_list(directory, amalgamation=False):
    tree_list=[]
    tree_file_list=construct_file_list(directory)
    for tree_file in tree_file_list:
        if amalgamation:
            #tree_list2=[read_tree(directory+tree_file)]
            tree_list2=read_mult_tree(directory+tree_file)
            tree_list.append(tree_list2)
            id_counter=0
            for tree in tree_list2:
                tree.tree_name=tree_file+str(id_counter)
                id_counter+=1
        else:
            tree=read_tree(directory+tree_file)
            tree.tree_name=tree_file
            tree_list.append(tree)
    return tree_list, tree_file_list

##  input matching dir or file
#return a  dict matching lower leaves to host leaves

#voir selon ce qu'on veut en faire apres
def construct_leaves_matching(matching_file):
    d=dict()

    #trying different names for the matching file (notably in the case of a dir where match needed between gene and their matching
    try:
        f = open(matching_file,"r")
    except IOError:
        s=matching_file
        new_match=s+".leafmap"
        try:
            f=open(new_match,"r")
        except IOError:
            if "." in s:
                new_match=s[:s.rindex(".")]+".leafmap"
            try:
                f=open(new_match,"r")
            except IOError:
                print("Matching file or files format non valid, if using dir, use gene file name + .leafmap or genefile name")

    s=f.read()
    s1=s.split(sep="\n")
    for  u in s1:
        if len(u)>0:
            s2=u.split(sep="\t")
            if not s2[0] in d:
                d[s2[0]]=[]
            if len(s2)==3:
                d[s2[0]].append([s2[1],s2[2]])
            else:
                if len(s2)==2:
                    d[s2[0]].append(s2[1])
                else:
                    print("Matching files format incorrect")
    return d

#file_list contain the name of the files, as explored for the trees, so that the gene list and matching lists correspond to one another
def construct_leaves_matching_dir(matching_directory, file_list):
    new_file_list=[]#some minor changes between leaf mapping and gene trees names
    for s in file_list:
        #print(s)
        new_file_list.append(s)
    leaf_matching_list=[]
    for i_file in range(len(new_file_list)):
        #pour correspondre aux noms donnés par sagephy
        d=construct_leaves_matching(matching_directory+new_file_list[i_file])
        leaf_matching_list.append(d)
    return leaf_matching_list


#à chaque nom de noeud renvoie le noeud correspondant
def construct_name_to_leaves(tree):
    d=dict()
    for u in tree.leaves():
        d[u.name]=u
    return d

def construct_name_to_clade(am_tree):
    d=dict()
    for u in am_tree.reverse_post_order:
        if u.is_leaf:
            d[u.clade_leaves[0]]=u
    return d

#à chaque nom d'arbre renvoie l'arbre corrspondant
def construct_name_to_tree(tree_list):
    d=dict()
    for t in tree_list:
        d[t.tree_name]=t
    return d

#always mult match output
def name_matching_to_tree_matching(name_to_node_host, name_to_node_lower,leaf_matching, tree=False):
    for lower_node_name in leaf_matching:
        if lower_node_name in name_to_node_lower : #sinon la feuille a été pruned de l'arbre
            lower_node=name_to_node_lower[lower_node_name]
            for host_info in leaf_matching[lower_node_name]:
                if len(host_info)==2:
                    host_node_name,host_tree_name=host_info
                    host_node=name_to_node_host[host_tree_name][host_node_name]
                else:
                    host_node_name=host_info
                    host_node=name_to_node_host[host_node_name]
                if tree:
                    lower_node.match=host_node
                else:
                    if lower_node.match is None:
                        lower_node.match=[]
                    lower_node.match.append(host_node)



def all_name_matching_to_tree_matching(symbiont_list, am_tree_list, leaf_matching_list):
    name_to_node_symbionts=dict()
    #ajouter un test pour savoir quel format est leaf matching list, arbre + feuille ou juste feuille ?
    for symbiont in symbiont_list:
        name_to_node_symbionts[symbiont.tree_name]=construct_name_to_leaves(symbiont)
    for clade_id in range(len(am_tree_list)):
        am_tree=am_tree_list[clade_id]
        lower_leaf_matching=leaf_matching_list[clade_id]
        name_to_node_lower=construct_name_to_clade(am_tree)
        name_matching_to_tree_matching(name_to_node_symbionts, name_to_node_lower, lower_leaf_matching)

#only one dict for all lower trees
def all_name_matching_to_tree_matching_onedict(symbiont_list, am_tree_list, lower_leaf_matching):
    name_to_node_symbionts=dict()
    for symbiont in symbiont_list:
        name_to_node_symbionts[symbiont.tree_name]=construct_name_to_leaves(symbiont)
    for clade_id in range(len(am_tree_list)):
        am_tree=am_tree_list[clade_id]
        name_to_node_lower=construct_name_to_clade(am_tree)
        name_matching_to_tree_matching(name_to_node_symbionts, name_to_node_lower, lower_leaf_matching)


#################################



#renomme tout les noeuds sauf les feuilles
def rename_tree(tree,root_name,n):
    if not tree.isLeaf():
        #if tree.name == None:
        tree.name=root_name+str(n)
        n+=1
        n=rename_tree(tree.left, root_name,n)
        n=rename_tree(tree.right, root_name,n)
    return n



def read_input_2levels(symbiont_directory, gene_directory=None, leaf_matching_directory=None, leaf_matching_file=None, inter=False, inter_list=None, inter_file_list=None):
    symbiont_list, symbiont_file_list=construct_tree_list(symbiont_directory, amalgamation=False)
    if gene_directory == None:
        gene_list, gene_file_list=inter_list, inter_file_list
    else:
        gene_list, gene_file_list=construct_tree_list(gene_directory, amalgamation=True)
    if inter:
        #non destructive construction of clades data list format for a rooted binary tree
        am_tree_list=compute_clade_frequencies_rooted_tree(inter_list)
    else:
        am_tree_list= compute_clade_frequencies_multiple_families(gene_list)
    if not leaf_matching_directory is None:
        leaf_matching_list=construct_leaves_matching_dir(leaf_matching_directory,gene_file_list)
        all_name_matching_to_tree_matching(symbiont_list,am_tree_list, leaf_matching_list)
    else:
        if not leaf_matching_file is None:
            leaf_matching=construct_leaves_matching(leaf_matching_file)
            all_name_matching_to_tree_matching_onedict(symbiont_list,am_tree_list, leaf_matching)
        else:
            symbiont_name_to_tree=dict()
            for symbiont in symbiont_list:
                for u in symbiont.leaves():
                    if not u.name in symbiont_name_to_tree:
                        symbiont_name_to_tree[u.name]=[]
                    symbiont_name_to_tree[u.name].append[u]
            for am_tree in am_tree_list:
                for u in am_tree.reverse_post_order:
                    if u.is_leaf():
                        u.match=symbiont_name_to_tree[u.clade_leaves[0]]
    for am_tree,gene_file in zip(am_tree_list,gene_file_list):
        am_tree.tree_name=gene_file
    return symbiont_list, am_tree_list

def rename_tree_list(tree_list,id_used,tree_name):
    i=0
    for tree in tree_list:
        id_used=rename_tree(tree,tree_name+str(i)+"_",id_used)
        i+=1
    return id_used

def rename_am_tree_list(am_tree_list,id_used,tree_name):
    i=0
    for am_tree in am_tree_list:
        for clade in am_tree.reverse_post_order:
            if not clade.is_leaf():
                clade.name=tree_name+str(i)+"_"+str(id_used)
                id_used+=1
        i+=1
    return id_used


def read_input(symbiont_directory, gene_directory, leaf_matching_directory=None, leaf_matching_file=None, host_directory=None, host_matching_file=None, host_matching_dir=None,inter_amalgamation=False):

    symbiont_list, am_tree_list=read_input_2levels(symbiont_directory, gene_directory=gene_directory, leaf_matching_directory=leaf_matching_directory, leaf_matching_file=leaf_matching_file, inter=False)

    id_used=0 #we give a unique id as a name to all non leaves nodes in all trees

    if host_directory is None:
        id_used=rename_tree_list(symbiont_list,id_used,"s")
        id_used=rename_am_tree_list(am_tree_list,id_used,"g")
        tree_list_symbiont_list=Tree_list(symbiont_list)

        return tree_list_symbiont_list, am_tree_list
    else:
        #we construct a new intermediate tree list, however,
        if inter_amalgamation:
            host_list, inter_am_tree_list=read_input_2levels(host_directory, inter_list=symbiont_list, leaf_matching_directory=host_matching_dir, leaf_matching_file=host_matching_file, inter=False,inter_file_list=symbiont_file_list)
        else:
            host_list, inter_am_tree_list=read_input_2levels(host_directory, inter_list=symbiont_list, leaf_matching_directory=host_matching_dir, leaf_matching_file=host_matching_file, inter=True, inter_file_list=symbiont_file_list)

        id_used=rename_tree_list(host_list,id_used,"h")
        id_used=rename_tree_list(symbiont_list,id_used,"st")
        id_used=rename_am_tree_list(inter_am_tree_list,id_used,"s")
        id_used=rename_am_tree_list(am_tree_list,id_used,"g")


        tree_list_host_list=Tree_list(host_list)

        return symbiont_list, am_tree_list, tree_list_host_list, inter_am_tree_list











### test
"""
data_dir="/home/hmenet/Documents/rewrite_ALE/test_erwinia/"
symbiont_dir="symbiont/"
gene_dir="genes/"
gene_match="lower_matching/"
file_links=[symbiont_dir, gene_dir, gene_match]
for i in range(len(file_links)):
    file_links[i]=data_dir+file_links[i]
symbiont_dir, gene_dir, gene_match=file_links


symbiont_list, clades_data_list, c_match_list=read_input(symbiont_dir, gene_dir, gene_match)
"""