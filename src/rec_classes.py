
from arbre import Tree

import numpy as np

############

class Event_rates:
    def __init__(self, tr=0.01,lr=0.01,dr=0.01):
        self.tr=tr #transfer rate
        self.dr=dr #duplication rate
        self.lr=lr #loss rate
        self.sr_compute() #speciation rate, ensuring sum of rates equal to 1
        self.log_rates()

    def log_rates(self): #initialize log rates for computation in log space
        self.ltr=np.log(self.tr)
        self.ldr=np.log(self.dr)
        self.llr=np.log(self.lr)
        self.lsr=np.log(self.sr)

    def sr_compute(self):
        self.sr=1-self.tr - self.dr - self.lr #speciation rate, ensuring sum of rates equal to 1

    def reinit(self):#recompute speciation rate and log rates, to use after modification of one of the rate
        self.sr_compute()
        self.log_rates()

    def pp(self):
        return "S " + self.sr + "D " + self.dr + "T " + self.tr + "L " + self.lr

#############

class Tree_list:
    def __init__(self,t_list):
        self.tree_list=t_list
        self.post_order_traversal()
        self.first_element=t_list[0]

    def post_order_traversal(self):
        post_order=[]
        for tree in self.tree_list:
            for u in tree.post_order_traversal():
                post_order.append(u)
        self.post_order=post_order

##############

class Rec_problem:
    def __init__(self, symb_list,amal_genes):
        self.upper = symb_list
        self.lower = amal_genes
        self.output_path = "output/"
        self.n_sample = 100
        self.n_steps = 5
        self.n_rec_sample_rates = 100
        self.n_recphyloxml=10
        self.best_rec = False
        self.rates = Event_rates(tr=0.01,lr=0.01,dr=0.01)
        self.ncpu = 4
        self.multiprocess_fam=False
        self.multiprocess_sampling=False

        ### for three level rec
        self.third_level=False
        self.upper_rec=None
        self.heuristic=None
        self.mc_samples=10


#################

class Amalgamated_tree:
    def __init__(self):
        self.clade_leaves=None
        self.child_frequencies=None
        self.tree_name=None
        self.match=None
        self.corresponding_tree=None
        self.reverse_post_order=None

    def aux_init(self,clade_elements,clade_frequencies,clade_id,d,corresponding_clade_to_tree=None):
        self.clade_leaves=clade_elements[clade_id]
        if not corresponding_clade_to_tree is None:
            self.corresponding_tree=corresponding_clade_to_tree(clade_id)
        new_clade_frequencies=dict()
        for (cl,cr) in clade_frequencies[clade_id]:
            new_cs=[]
            for c in [cl,cr]:
                if not c in d:
                    new_c=Amalgamated_tree()
                    new_c.aux_init(clade_elements,clade_frequencies,c,d)
                    d[c]=new_c
                else:
                    new_c=d[c]
                new_cs.append(new_c)
            newcl,newcr=new_cs
            new_clade_frequencies[(newcl,newcr)]=clade_frequencies[clade_id][(cl,cr)]
        self.child_frequencies=new_clade_frequencies


    def initialize(self,clade_elements,clade_frequencies,corresponding_clade_to_tree=None):
        d=dict()
        d[0]=self
        self.aux_init(clade_elements,clade_frequencies,0,d,corresponding_clade_to_tree=corresponding_clade_to_tree)
        self.reverse_post_order_traversal()
        self.leaves_traversal()

    def is_leaf(self):
        return(len(clade_leaves))==1

    def aux_liste(self,l,d):
        for (cl,cr) in self.child_frequencies:
            for c in [cl,cr]:
                if not c in d:
                    d[c]=0
                    l.append(c)
                    c.aux_liste(l,d)

    def liste(self):
        d=dict()
        l=[self]
        self.aux_list(l,d)
        return l

    def reverse_post_order_traversal(self):
        def order_clade(clade):
            return(len(clade.clade_leaves))
        l=self.liste()
        l.sort(key=order_clade, reverse=True)
        #ordre donné par la taille des clades est bien un ordre qui empeche les fils d'apparaitre avant leur parent
        self.reverse_post_order=l
        return l

    def leaves_traversal(self):
        if self.reverse_post_order is None:
            self.reverse_post_order_traversal()
        l=[u for u in self.reverse_post_order if u.is_leaf()]
        self.leaves=l
        return l


##################$$


class Rec_event:
    event_name
    clade
    right_clade
    left_clade
    right_upper
    left_upper
                    event_name, f,v,g,w = event
                    r[v].append(f)
                    r[w].append(g)
                    clade_to_look.append(v)
                    clade_to_look.append(w)
                    event_to_append=(event_name,e,c,f,v,g,w)


class Rec_scenario
    event_list
    reconstructed_tree
    ancestral_matching

class Two_level_rec_sol

    likelihood,l_event_gene, l_scenarios=out_rec

    l_event_by_fam



class Three_level_rec_sol

    likelihood, l_event_gene, l_scenarios, l_scenarios_upper, log_likelihood_list
