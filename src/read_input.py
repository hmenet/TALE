### Prepare the input, from repo containing trees in newick (or list of tree in newick) to tree class and/or clades amalgamation format

from os import walk
import numpy as np
from read_tree import read_tree, read_mult_tree
from read_clade_frequencies import compute_clade_frequencies_multiple_families, compute_clade_frequencies_rooted_tree



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
#return a  dict matching inner leaves to host leaves

#voir selon ce qu'on veut en faire apres
def construct_leaves_matching(matching_file):
    d=dict()
    f=open(matching_file, "r")
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
        new_file_list.append(s[:-3]+"leafmap")
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

#always mult match output
def name_matching_to_tree_matching(name_to_node_host, name_to_node_inner,leaf_matching, tree=False):
    if not tree:
        d=dict()
    for inner_node_name in leaf_matching:
        if inner_node_name in name_to_node_inner : #sinon la feuille a été pruned de l'arbre
            inner_node=name_to_node_inner[inner_node_name]
            for host_info in leaf_matching[inner_node_name]:
                if len(host_info)==2:
                    host_node_name,host_tree_name=host_info
                    host_node=name_to_node_host[host_tree_name][host_node_name]
                    if tree:
                        inner_node.match=host_node
                    else:
                        if not inner_node in d:
                            d[inner_node]=[]
                        d[inner_node].append(host_node)
                else:
                    host_node_name=host_info
                    host_node=name_to_node_host[host_node_name]
                    inner_node.match=host_node
    if not tree:
        return d



def all_name_matching_to_tree_matching(symbiont_list, clades_data_list, leaf_matching_list):
    c_match_list=[]
    name_to_node_symbionts=dict()
    #ajouter un test pour savoir quel format est leaf matching list, arbre + feuille ou juste feuille ?
    for symbiont in symbiont_list:
        name_to_node_symbionts[symbiont.tree_name]=construct_name_to_leaves(symbiont)
    for clade_id in range(len(clades_data_list)):
        clade_keys=clades_data_list[clade_id][3]
        lower_leaf_matching=leaf_matching_list[clade_id]
        name_to_node_inner=construct_name_to_clade(clade_keys)
        d=name_matching_to_tree_matching(name_to_node_symbionts, name_to_node_inner, lower_leaf_matching)
        c_match_list.append(d)
    return c_match_list

#only one dict for all lower trees
def all_name_matching_to_tree_matching_onedict(symbiont_list, clades_data_list, lower_leaf_matching):
    c_match_list=[]
    name_to_node_symbionts=dict()
    for symbiont in symbiont_list:
        name_to_node_symbionts[symbiont.tree_name]=construct_name_to_leaves(symbiont)
    for clade_id in range(len(clades_data_list)):
        clade_keys=clades_data_list[clade_id][3]
        name_to_node_inner=construct_name_to_clade(clade_keys)
        d=name_matching_to_tree_matching(name_to_node_symbionts, name_to_node_inner, lower_leaf_matching)
        c_match_list.append(d)
    return c_match_list


#################################



#renomme tout les noeuds sauf les feuilles
def rename_tree(tree,root_name):
    def aux_rename_tree(tree, root_name, n):
        if not tree.isLeaf():
            tree.name=root_name+str(n)
            n+=1
            n=aux_rename_tree(tree.left, root_name,n)
            n=aux_rename_tree(tree.right, root_name,n)
        return n
    aux_rename_tree(tree,root_name,0)



def read_input_2levels(symbiont_directory, gene_directory=None, leaf_matching_directory=None, leaf_matching_file=None, inter=False, inter_list=None, inter_file_list=None):
    symbiont_list, symbiont_file_list=construct_tree_list(symbiont_directory, amalgamation=False)
    if gene_directory == None:
        gene_list, gene_file_list=inter_list, inter_file_list
    else:
        gene_list, gene_file_list=construct_tree_list(gene_directory, amalgamation=True)
    i=0
    for symbiont in symbiont_list:
        rename_tree(symbiont,symbiont_file_list[i])
        i+=1
    if inter:
        #non destructive construction of clades data list format for a rooted binary tree
        clades_data_list, clade_to_tree=compute_clade_frequencies_rooted_tree(inter_list)
    else:
        clades_data_list= compute_clade_frequencies_multiple_families(gene_list)
    if leaf_matching_directory != None:
        leaf_matching_list=construct_leaves_matching_dir(leaf_matching_directory,gene_file_list)
        c_match_list=all_name_matching_to_tree_matching(symbiont_list,clades_data_list, leaf_matching_list)
    else:
        if leaf_matching_file != None:
            leaf_matching=construct_leaves_matching(leaf_matching_file)
            c_match_list=all_name_matching_to_tree_matching_onedict(symbiont_list,clades_data_list, leaf_matching)
        else:
            symbiont_name_to_tree=dict()
            for symbiont in symbiont_list:
                for u in symbiont.leaves():
                    symbiont_name_to_tree[u.name]=u
            c_match_list=[]
            for i_clade in range(len(clades_data_list)):
                d=dict()
                for u in clades_data_list[i_clade][0]:
                    v=clades_data_list[i_clade][2][u]
                    if len(v)==1:
                        d[u]=symbiont_name_to_tree[v[0]]
                c_match_list.append(d)
    if inter:
        return symbiont_list, clades_data_list, c_match_list, gene_file_list, clade_to_tree
    else:
        return symbiont_list, clades_data_list, c_match_list, gene_file_list, symbiont_file_list


def read_input(symbiont_directory, gene_directory, leaf_matching_directory=None, leaf_matching_file=None, host_directory=None, host_matching_file=None, host_matching_dir=None):

    symbiont_list, clades_data_list, c_match_list, gene_file_list, symbiont_file_list=read_input_2levels(symbiont_directory, gene_directory=gene_directory, leaf_matching_directory=leaf_matching_directory, leaf_matching_file=leaf_matching_file, inter=False)
    if host_directory==None:
        return symbiont_list, clades_data_list, c_match_list, gene_file_list
    else:
        #we construct a new intermediate tree list, however,
        host_list, inter_clades_data_list, inter_c_match_list, inter_file_list, inter_clade_to_tree=read_input_2levels(host_directory, inter_list=symbiont_list, leaf_matching_directory=host_matching_dir, leaf_matching_file=host_matching_file, inter=True, inter_file_list=symbiont_file_list)
        return symbiont_list, clades_data_list, c_match_list, gene_file_list, host_list, inter_clades_data_list, inter_c_match_list, inter_file_list, inter_clade_to_tree







"""
def read_input_3Level(upper_dir, inter_dir, lower_dir, int_up_match_dir=None, low_int_match_dir=None):
    upper_list, upper_file_list=construct_tree_list(upper_dir, amalgamation=False)

    inter_list, inter_file_list=construct_tree_list(inter_dir)

    lower_list, lower_file_list=


    symbiont_list, symbiont_file_list=construct_tree_list(symbiont_directory, amalgamation=False)
    gene_list, gene_file_list=construct_tree_list(gene_directory, amalgamation=True)
    i=0
    for symbiont in symbiont_list:
        rename_tree(symbiont,symbiont_file_list[i])
        i+=1
    clades_data_list= compute_clade_frequencies_multiple_families(gene_list)
    if leaf_matching_directory != None:
        leaf_matching_list=construct_leaves_matching_dir(leaf_matching_directory,gene_file_list)
        c_match_list=all_name_matching_to_tree_matching(symbiont_list,clades_data_list, leaf_matching_list)
    else:
        symbiont_name_to_tree=dict()
        for symbiont in symbiont_list:
            for u in symbiont.leaves():
                symbiont_name_to_tree[u.name]=u
        c_match_list=[]
        for i_clade in range(len(clades_data_list)):
            d=dict()
            for u in clades_data_list[i_clade][0]:
                v=clades_data_list[i_clade][2][u]
                if len(v)==1:
                    d[u]=symbiont_name_to_tree[v[0]]
            c_match_list.append(d)
    return symbiont_list, clades_data_list, c_match_list, gene_file_list

"""





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