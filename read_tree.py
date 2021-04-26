#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 15:45:21 2020

@author: hmenet
"""

import arbre

################## standard newick input ####################



def read_tree_string(s, return_end_point=False, starting_point=0):
    tree=arbre.Tree()
    tree.id=1
    tree.root=tree
    node=tree
    k=starting_point
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
            if return_end_point:
                #print(s[k],k,len(s))
                #print("sk", s[k:])
                return tree,k+2
            else:
                return tree
        else:
            name=""
            while not c in ["(",")",",",":", ";"]:
                name+=c
                k+=1
                c=s[k]
            if c==":":
                while not c in ["(",")",",",";"]:
                    k+=1
                    c=s[k]
            if len(name)>0:
                node.name=name
            k-=1
        k+=1
    #if return_end_point:
    #    return (tree,k)
    #else:
    #    return tree

#return a tree
def read_tree(file):
    f=open(file, "r")
    s=f.read()
    return read_tree_string(s)

#liste d'arbre pour amalgamation, tous dans le même fichier
def read_mult_tree(file):
    f=open(file, "r")
    s=f.read()
    list_tree=[]
    k=0
    while k<len(s):
        #read_tree_string(s,return_end_point=True, starting_point=k)
        t,k=read_tree_string(s,return_end_point=True, starting_point=k)
        #print(k, len(s))
        list_tree.append(t)
    return list_tree

################## sagephy input ###########################

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









