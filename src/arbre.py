#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 10:22:42 2019

@author: hugo
"""

# Definition of the tree class used



class Tree:
    def __init__(self):
        self.left = None #left child
        self.right = None #right child
        self.right2 = None # to consider unrooted trees
        self.parent=None #point to the parent
        self.root=self #point to the root of the tree
        self.isUnrooted=False #say if the tree is unrooted, by default the tree is rooted
        
        
        self.isHostTree=False
        self.isParasiteTree=False
        self.isGeneTree=False
        self.isFreeLiving=False
        
        self.added_for_free_living=False
        
        #attributes used for simulation
        #not known for real data
        self.time = None #end time of the specie (of its branch)
        self.match=None #for a parasite, the host, for a gene, the specie (that can be parasite or host)
        self.birth_time=None #birth time of the specie
        self.parasites=[] # for a host, list of the parasites in this host (réciproque de match)
        
        self.event_list=None
        
        ###for visualisation and manipulation
        self.id=1 #an id to identificate the nodes, root is 1, if parent is n, children are 2n and 2n+1
        self.name=None # a name for the tree,to identificate trees with the root name and when the tree comes from real data its the name of the node in the newick file
        self.tree_name=None
        
        self.pruned_name_list=[]
        
        ###for prune (when one of the child of a node is a leaf and do not reach contemporary time 1)
        self.inTree=1 #says if the node is still in the tree or not when pruning
        
        #rates of branches used for simulation and for computation of likelihood
        self.birth_rate=None #rate of 
        self.death_rate=None # //
        self.transfer_rate=None #//
        self.speciation_rate=None #speciation de l'arbre
        self.parasite_birth_rate=None
        self.parasite_death_rate=None
        self.parasite_transfer_rate=None
        self.parasite_speciation_rate=None
        
    def isLeaf(self):
        return self.right == None
    
    def isRoot(self):
        return self.parent == None

    def isRooted(self):
        return self.right2 is None

    #return true if self is a strict ascendant of tree
    def isAscendant(self, tree):
        if tree.isRoot():
            return False
        else:
            if tree.parent == self:
                return True
            else:
                return self.isAscendant(tree.parent)
    
    #return the list of all nodes of the subtree starting at self
    def liste(self):
        if self.isLeaf():
            return [self]
        else:
            if self.right2 is None:
                return [self]+self.left.liste() + self.right.liste()
            else:
                return [self]+self.left.liste() + self.right.liste()+self.right2.liste()
    
    #return the lsit of all leaves of the subtree starting at self
    def leaves(self):
        if self.isLeaf():
            return [self]
        else:
            if self.right2 is None:
                return self.right.leaves() + self.left.leaves()
            else: 
                return self.right.leaves() + self.left.leaves() + self.right2.leaves()
        
    def leaves_unrooted(self):
        if self.isLeaf():
            return [self]
        else:
            if self.isRoot(): 
                return self.right.leaves() + self.left.leaves() + self.right2.leaves()
            else:
                return self.right.leaves() + self.left.leaves()
    
    def leaves_unrooted2(self):
        if self.isLeaf():
            return [self]
        else:
            if self.right2==None:
                return self.right.leaves_unrooted2() + self.left.leaves_unrooted2()
            else:
                return self.right.leaves_unrooted2() + self.left.leaves_unrooted2() + self.right2.leaves_unrooted2()
    
    ################### Function for simulation ######################

    #return the list of existing species at time t in the subtree starting at self
    def existing_species(self,t):
        if self.birth_time <= t and t <= self.time:
            return [self]
        elif self.time<t:
            if (not self.isLeaf()):
                return self.left.existing_species(t) + self.right.existing_species(t)
            else:
                return []
        else:
            return []
    
    #return the list of parasites in self at time t
    def parasiting_species(self,t):
        l=self.parasites
        current_parasites=[]
        for p in l:
            if p.birth_time < t and t <= p.time:
                current_parasites.append(p)
        return current_parasites

    #create a left son and a right son of birth time t, and set time of death of self at t
    def birth(self,t):
        self.right = Tree()
        self.left = Tree()
        self.right.root=self.root
        self.left.root=self.root
        self.time = t     
        self.right.birth_time=t
        self.left.birth_time=t
        self.right.id=2*self.id +1 
        self.left.id=2*self.id
        self.right.parent=self
        self.left.parent=self

    def undated_birth(self):
        self.right = Tree()
        self.left = Tree()
        self.right.root=self.root
        self.left.root=self.root
        self.right.id=2*self.id +1 
        self.left.id=2*self.id
        self.right.parent=self
        self.left.parent=self


    def birth_specie(self, t):
        self.birth(t)
        
    #birth + set the children match to be the same as the father, and add the children to the parasites of the match
    def birth_parasite(self, t):
        self.birth(t)
        self.right.match=self.match
        self.left.match=self.match
        self.match.parasites.append(self.left)
        self.match.parasites.append(self.right)
    
    #birth + set the children match to be the same as the father
    def birth_gene(self, t):
        self.birth(t)
        self.right.match=self.match
        self.left.match=self.match
    
    #set the death time of self at t
    def death(self, t):
        self.time=t

    #cospeciation of parasite and host
    def cospeciation_parasite(self,t):
        self.birth_parasite(self.match.time)
        self.time=self.match.time
        self.left.match.parasites.remove(self.left)
        self.left.match=self.match.left
        self.left.match.parasites.append(self.left)
        self.right.match.parasites.remove(self.right)
        self.right.match=self.match.right
        self.right.match.parasites.append(self.right)
    
    #cospeciation of gene and the species to which it is matched
    def cospeciation_gene(self,t):
        self.birth_gene(self.match.time)
        self.time=self.match.time
        self.left.match=self.match.left
        self.right.match=self.match.right
        
    def transfer_parasite(self, t, i): #time t, species i
        self.birth_parasite(t)
        self.right.match.parasites.remove(self.right)
        self.right.match=i
        self.right.match.parasites.append(self.right)
        
    def transfer_gene(self,t,i):
        self.birth_gene(t)
        self.right.match=i
        
    def reset_parasites(self):
        self.parasites=[]
        if not self.isLeaf():
            self.left.reset_parasites
            self.right.reset_parasites
    
    #delete all leaves that do not reach time t = 1
    #and fuse the nodes with only one child to their child
    def prune_with_match(self):
        self.prune_tree()
        self.aux_prune_with_match()
    
    def aux_prune_with_match(self):
        if self.match.inTree==-1:
            self.match=self.match.parent
            self.aux_prune_with_match()
        else:
            if not self.isLeaf():
                self.left.aux_prune_with_match()
                self.right.aux_prune_with_match()   
        
    def chemin(self, noeud_descendant):
        if self == noeud_descendant:
            return [noeud_descendant]
        else:
            c=self.chemin(noeud_descendant.parent)
            c.append(noeud_descendant)
            return c

    def copy(self):
        arbre=Tree()
        self.aux_copy(arbre)
        arbre.birth_time = self.birth_time
        return arbre
    
    def aux_copy(self, arbre):
        arbre.time=self.time
        arbre.name=self.name
        arbre.pruned_name_list=list(self.pruned_name_list)
        arbre.tree_name=self.tree_name
        if not self.isLeaf():
            arbre.birth(self.time)
            self.left.aux_copy(arbre.left)
            self.right.aux_copy(arbre.right)
            
    def copy_lower_match(self):
        arbre=Tree()
        self.aux_copy_lower_match(arbre)
        arbre.birth_time = self.birth_time
        return arbre
    
    def aux_copy_lower_match(self, arbre):
        arbre.time=self.time
        arbre.name=self.name
        arbre.pruned_name_list=list(self.pruned_name_list)
        arbre.tree_name=self.tree_name
        arbre.match=self.match
        if not self.isLeaf():
            arbre.birth(self.time)
            self.left.aux_copy_lower_match(arbre.left)
            self.right.aux_copy_lower_match(arbre.right)
    
    #to copy the host, correspondances between the original and the host have to be stored, to then be able to match the nodes of the copy to the nodes matched to the original
    def copy_host(self):
        d=dict()
        arbre=Tree()
        arbre.name=self.name+"copie"
        arbre.birth_time = self.birth_time
        d[self]=arbre
        self.aux_copy_host(arbre,d)
        return arbre,d
    
    def aux_copy_host(self, arbre, d):
        d[self]=arbre
        arbre.time=self.time
        if not self.isLeaf():
            arbre.birth(self.time)                
            self.left.aux_copy_host(arbre.left, d)
            self.right.aux_copy_host(arbre.right, d)
    
    #take as input a dictionary of the correspondances between copy and original
    #works for genes and symbiotes levels
    def copy_with_match(self, d):
        arbre=Tree()
        arbre.name=self.name+"copie"
        arbre.birth_time = self.birth_time
        d[self]=arbre
        self.aux_copy_with_match(arbre,d)
        return arbre
    
    def aux_copy_with_match(self, arbre,d):
        d[self]=arbre
        arbre.time=self.time
        arbre.match=d[self.match]
        if not self.isLeaf():
            arbre.birth(self.time)                
            self.left.aux_copy_with_match(arbre.left, d)
            self.right.aux_copy_with_match(arbre.right, d)
    
    #return a list of the nodes of the tree in a post order traversal
    #we will use it with pop, so in fact it's the reverse of a post order
    def post_order_traversal(self):
        if self.right2 is None:
            t=[self]
            if not self.isLeaf():
                t=t+self.left.post_order_traversal()+self.right.post_order_traversal()
        else:
            return self.post_order_traversal_unrooted()
        return t
    
    def post_order_traversal_unrooted(self):
        t=[self]
        if not self.isLeaf():
            if not self.isRoot():
                t=t+self.left.post_order_traversal()+self.right.post_order_traversal()
            else:
                t=t+self.left.post_order_traversal()+self.right.post_order_traversal()+self.right2.post_order_traversal()
        return t

    def sibling(self):
        if self.parent.left==self:
            return self.parent.right
        if self.parent.right==self:
            return self.parent.left
    
    #return the list of ancestors of a node
    def ancestor(self):
        if self.isRoot():
            return [self]
        else:
            l=self.parent.ancestor()
            l.append(self)
            return l
            
    def ancestor_dict(self):
        if self.isRoot():
            d=dict()
            d[self]=True
            return d
        else:
            d=self.parent.ancestor_dict()
            d[self]=True
            return d

        
    #modify the input     
    def prune_tree(self):
        #condition to not launch on a node already deleted
        if not(self.inTree==-1):
            if not self.isLeaf():
                if self.left.isLeaf() and self.left.time <1:
                    if self.isRoot():
                        if self.right.isLeaf():
                            self.graphviz_tree()
                            raise NameError("On prune et il ne reste que la racine et une feuille")
                        else:
                            self.time=self.right.time
                            self.pruned_name_list=self.pruned_name_list + self.right.pruned_name_list
                            #delete self.left et self.right
                            self.left.inTree=-1
                            self.right.inTree=-1
                            self.left=self.right.left
                            self.right=self.right.right
                            self.left.parent=self
                            self.right.parent=self
                            self.prune_tree()
                    else:
                        self.right.birth_time=self.birth_time
                        self.right.parent=self.parent
                        if self.parent.left == self:
                            self.parent.left=self.right
                        else:
                            self.parent.right=self.right
                        #delete self et self.left
                        self.inTree=-1
                        self.left.inTree=-1
                        self.right.pruned_name_list=self.right.pruned_name_list + self.pruned_name_list
                        self.parent.prune_tree()
                else:
                    if self.right.isLeaf() and self.right.time <1:
                        if self.isRoot():
                            if self.left.isLeaf():
                                self.graphviz_tree()
                                raise NameError("Prune but only the root and a leaf remaining")
                            else:
                                self.time=self.left.time
                                self.pruned_name_list=self.pruned_name_list + self.left.pruned_name_list
                                self.left.inTree=-1
                                self.right.inTree=-1
                                self.right=self.left.right
                                self.left=self.left.left
                                self.left.parent=self
                                self.right.parent=self
                                self.prune_tree()
                        else:
                            self.left.birth_time=self.birth_time
                            self.left.parent=self.parent
                            if self.parent.left == self:
                                self.parent.left=self.left
                            else:
                                self.parent.right=self.left
                            self.inTree=-1
                            self.right.inTree=-1
                            self.left.pruned_name_list=self.left.pruned_name_list + self.pruned_name_list
                            self.parent.prune_tree()
                    else:
                        self.left.prune_tree()
                        self.right.prune_tree()
    
    def name_to_tree(self):
        d=dict()
        for u in self.liste():
            d[u.name]=u
        return d
        
    def n_ancestor(self,n):
        if n==0:
            return self
        else:
            return self.parent.n_ancestor(n-1)
        
    #could be optimized
    def n_leaves_ancestor(self,n):
        if len(self.leaves())==n:
            return self
        else:
            return (self.parent).n_leaves_ancestor(n)



#node1 is ascendant (not strictly) to node 2, we want the distance between them in the tree
def distance_ascendant(node1, node2):
    if node1 == node2:
        return 0
    else:
        return 1 + distance_ascendant(node1,node2.parent)

#could be optimized
def distance(node1, node2):
    if node1.root==node2.root:
        #because isAscendant is strict
        if node1==node2:
            return 0
        if node1.isAscendant(node2):
            return distance_ascendant(node1, node2)
        else:
            if node2.isAscendant(node1):
                return distance_ascendant(node2, node1)
            else: #none of the nodes can be the root as one would be ascendant to the other in this case
                return 1 + distance(node1, node2.parent)
    else:
        #we do as if the two trees were joined at their roots
        return distance(node1,node1.root)+distance(node2,node2.root)+2
    
def tree_construct_real_match(tree):
    tmp=dict()
    match=dict()
    for i in tree.liste():
        match[i]=i.match
    for i in tree.liste():
        tmp[i]=[]
        if not i.isRoot():
            if match[i.parent].isAscendant(match[i]):
                if i.parent.left==i:
                    bro=i.parent.right
                else:
                    bro=i.parent.left
                c=match[i.parent].chemin(match[i])
                if match[i.parent].isAscendant(match[bro]):
                    c_bro=match[i.parent].chemin(match[bro])
                if not(match[i.parent].isAscendant(match[bro])) or c_bro[1]==c[1]: 
                    c=match[i.parent].chemin(match[i])
                    tmp[i]+=c[0:len(c)-1]
                else:
                    if distance_ascendant(match[i.parent], match[i])>1:
                        c=match[i.parent].chemin(match[i])
                        tmp[i]+=c[1:len(c)-1]
        tmp[i].append(match[i])
    return tmp




        

def aux_save_tree(f,tree):
    if not tree.isLeaf():
        f.write("(")
        aux_save_tree(f, tree.left)
        f.write(",")
        aux_save_tree(f,tree.right)
        if not tree.right2==None:
            f.write(",")
            aux_save_tree(f,tree.right2)
        f.write(")")
    if not tree.name==None:
        f.write(tree.name)

#save node names and tree topology in newick
def save_tree(tree, file):
    f=open(file, "w")
    aux_save_tree(f, tree)
    f.write(";")


def read_tree_string(s):
    tree=Tree()
    tree.id=1
    tree.root=tree
    node=tree
    k=0
    while k < len(s):
        c=s[k]
        if c=="(":
            node.left=Tree()
            node.left.parent=node
            node.left.root=node.root
            node.left.id=node.id*2
            node=node.left
        elif c==",":
            node=node.parent
            
            if node.right== None:
                node.right=Tree()
                node.right.parent=node
                node.right.root=node.root
                node.right.id=node.id*2+1
                node=node.right
            else:
                if node.right2==None:
                    node.right2=Tree()
                    node.right2.parent=node
                    node.right2.root=node.root
                    node.right2.id=node.id*2+1.5#additional node to model unrooted binary tree 
                    node=node.right2
                else:
                    print("Wasn't expecting multifurcating tree")
        elif c==")":
            node=node.parent
        elif c==";":
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
    return tree

#return a tree
def read_tree(file):
    f=open(file, "r")
    s=f.read()
    return read_tree_string(s)
    
def from_rooted_to_un(tree):
    if len(tree.leaves())<3 or not tree.right2==None:
        print("Au moins 3 feuilles dans from rooted to unrooted")
    else: 
        if not tree.left.isLeaf():
            tree.right2=tree.left.left
            tree.right2.parent=tree
            tree.left=tree.left.right
            tree.left.parent=tree
        else:
            tree.right2=tree.right.left
            tree.right2.parent=tree
            tree.right=tree.right.right
            tree.left.parent=tree
            
            
            
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
            
            
def read_tree_string_sagephy(s):
    nodes_info=dict()
    tree=Tree()
    tree.id=1
    tree.root=tree
    node=tree
    k=0
    while k < len(s):
        c=s[k]
        if c=="(":
            node.left=Tree()
            node.left.parent=node
            node.left.root=node.root
            node.left.id=node.id*2
            node=node.left
        elif c==",":
            node=node.parent
            
            if node.right== None:
                node.right=Tree()
                node.right.parent=node
                node.right.root=node.root
                node.right.id=node.id*2+1
                node=node.right
            else:
                if node.right2==None:
                    node.right2=Tree()
                    node.right2.parent=node
                    node.right2.root=node.root
                    node.right2.id=node.id*2+1.5#unrooted binary tree
                    node=node.right2
                else:
                    print("Wasn't expecting multifurcating tree")
        elif c==")":
            node=node.parent
        elif c==";":
            return tree, nodes_info
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
                    nodes_info[node]=l_word2
                    
                
            if len(name)>0:
                node.name=name
            k-=1
        k+=1

#return a tree
def read_tree_sagephy(file):
    f=open(file, "r")
    s=f.read()
    return read_tree_string_sagephy(s)





