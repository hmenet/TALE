

#simulation for trees mostly universal unicopy, and 3 level



import numpy as np
import random as rd

from arbre import Tree


def construct_species_tree(species_model):
    [birth_rate,death_rate]=species_model
    rate=birth_rate+death_rate
    root=Tree()
    root.birth_time=0
    leaves = [root]
    while len(leaves)>0:
        i=leaves.pop()
        t1=i.birth_time + np.random.exponential(scale=1/rate)
        if t1<1:
            x=rd.random()
            if x>(birth_rate/rate):
                i.death(t1)
            else:
                i.birth_specie(t1)
                leaves.append(i.left)
                leaves.append(i.right)
        else:
            i.time=1
    return root


def construct_parasite_inside_host_tree(parasites_model, host_tree, random_root=True):
    [transfer_rate, birth_rate, death_rate]=parasites_model
    rate=birth_rate+death_rate+transfer_rate
    root=Tree()
    #t0=rd.random()
    n=len(host_tree.liste())
    #par_exp=(1/n)**0.7
    #par_exp=np.exp(-n**0.7)
    #par_exp=1/np.log(n)
    if random_root:
        par_exp=0.5
        x=np.random.exponential(par_exp)
        while x>1:
            x=np.random.exponential(par_exp)
        t0=x
        #t0=rd.random()
        l=host_tree.existing_species(t0)
        if len(l)==0:
            raise NameError("L'arbre d'espece hote n'atteint pas le temps t=1")
        k=rd.randint(0, len(l)-1)
        i=l[k]
    else:
        t0=0
        i=host_tree
    #choose a node at random in the species tree (uniform ? inverse proportional to the depth ?)
    root.birth_time=t0
    root.match = i
    i.parasites.append(root)
    leaves=[root]
    while len(leaves)>0:
        i=leaves.pop()
        t1=i.birth_time + np.random.exponential(scale=1/rate)
        if t1 < i.match.time:
            x=rd.random()
            if x>((birth_rate+death_rate)/rate):
                #transfer
                #j=choose from all existing species at time t1 (uniform ? depend on distance ?)
                #existing_species fait un parcours complet de l'arbre jusqu'au temps t1
                #pour faire mieux on pourrait faire attention à dans quel ordre on parcourt les noeuds et tenir à jour une liste de toutes les especes encore vivantes
                l=host_tree.existing_species(t1)
                k=rd.randint(0, len(l)-1)
                j=l[k]
                i.transfer_parasite(t1,j)
                leaves.append(i.left)
                leaves.append(i.right)
                #print("transfer")
            elif x>birth_rate/rate:
                i.death(t1)
                #print("death")
            else:
                i.birth_parasite(t1)
                leaves.append(i.left)
                leaves.append(i.right)
                #print("birth")
        else:
            t1=i.match.time
            if i.match.isLeaf():
                i.death(t1)
                #print("species death/end of time")
            else:
                i.cospeciation_parasite(t1)
                leaves.append(i.left)
                leaves.append(i.right)
                #print("cospeciation")
    return root

def construct_gene(symbiont_list):

    root_choice()
    to_see=[gene]


    while len(to_see) > 0:

        gene=to_see.pop()
        gene.next_time()
        if gene.next_event() = T


