
from arbre import Tree

import numpy as np

############

class Event_rates:
    def __init__(self, tr=0.01,lr=0.01,dr=0.01,ir=0):
        self.tr=tr #transfer rate
        self.dr=dr #duplication rate
        self.lr=lr #loss rate
        self.ir=ir #incomplete sorting rate
        self.sr_compute() #speciation rate, ensuring sum of rates equal to 1
        self.log_rates()

    def log_rates(self): #initialize log rates for computation in log space
        self.ltr=np.log(self.tr)
        self.ldr=np.log(self.dr)
        self.llr=np.log(self.lr)
        self.lsr=np.log(self.sr)
        if self.ir > 0:
            self.lir=np.log(self.ir)
        else:
            self.lir=None

    def sr_compute(self):
        self.sr=1-self.tr - self.dr - self.lr -self.ir #speciation rate, ensuring sum of rates equal to 1

    def reinit(self):#recompute speciation rate and log rates, to use after modification of one of the rate
        self.sr_compute()
        self.log_rates()

    def pp(self):
        if self.ir==0:
            return "S " + str(self.sr) + "\tD " + str(self.dr) + "\tT " + str(self.tr) + "\tL " + str(self.lr)
        else:
            return "S " + str(self.sr) + "\tD " + str(self.dr) + "\tT " + str(self.tr) + "\tL " + str(self.lr) +  "\tI " + str(self.ir)

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


    def __iter__(self):
        return self.tree_list.__iter__()

    def __next__(self):
        return self.tree_list.__next__()


    def append(self,new_tree):
        self.tree_list.append(new_tree)
        self.post_order_traversal()

    def name_to_tree(self):
        d=dict()
        for tree in self.tree_list:
            tmpd=tree.name_to_tree()
            for k in tmpd:
                d[k]=tmpd[k]
        return d

##############


class Rec_upper_tree_computation:
    def __init__(self):
        self.E=None
        self.Eavg_no_log=None
        self.P_transfer=None

class Rec_lower_tree_computation:
    def __init__(self):
        self.P=None
        self.P_TL=None
        self.log_l=None
        self.corr_size=None
        self.time=None

##############

class Rec_problem:
    def __init__(self, symb_list,amal_genes):
        self.upper = symb_list
        self.lower = amal_genes
        self.output_path = "output/"
        self.n_sample = 100
        self.n_steps = 5
        self.n_rec_sample_rates = 100
        self.n_output_scenario=10
        self.best_rec = False
        self.rates = Event_rates(tr=0.01,lr=0.01,dr=0.01)
        self.ncpu = 4
        self.multiprocess_fam=False
        self.multiprocess_sampling=False
        self.upper_tree_computation=Rec_upper_tree_computation()
        self.lower_tree_computation=Rec_lower_tree_computation()
        self.rate_inference=False #True when infering rates
        self.dd=False#distance dependent transfers

        ### for three level rec
        self.third_level=False
        self.upper_rec=None
        self.heuristic=None
        self.mc_samples=10




    def __iter__(self):
        self.iter_lower=0
        return self

    def __next__(self):
        if self.iter_lower<len(self.lower):


            d=vars(self)

            single_gene_rec=Rec_problem(self.upper,None)

            newd=vars(single_gene_rec)

            for key in d:
                if key == "lower":
                    newd["single_lower"]=self.lower[self.iter_lower]
                else:
                    if key != "single_lower":
                        newd[key]=d[key]

            self.iter_lower+=1
            return single_gene_rec
        else:
            raise StopIteration





#################

class Amalgamated_tree:
    def __init__(self):
        self.clade_leaves=None
        self.child_frequencies=None
        self.tree_name=None
        self.match=None
        self.corresponding_tree=None
        self.reverse_post_order=None
        self.leaves=None
        self.constant_match=None
        self.name=None
        self.root=self

    def aux_init(self,clade_elements,clade_frequencies,clade_id,d,corresponding_clade_to_tree=None):
        self.clade_leaves=clade_elements[clade_id]
        if not corresponding_clade_to_tree is None:
            self.corresponding_tree=corresponding_clade_to_tree[clade_id]
        new_clade_frequencies=dict()
        for (cl,cr) in clade_frequencies[clade_id]:
            new_cs=[]
            for c in [cl,cr]:
                if not c in d:
                    new_c=Amalgamated_tree()
                    new_c.root=self.root
                    new_c.aux_init(clade_elements,clade_frequencies,c,d,corresponding_clade_to_tree=corresponding_clade_to_tree)
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
        self.log_freq()

    def is_leaf(self):
        return(len(self.clade_leaves))==1

    def aux_liste(self,l,d):
        for (cl,cr) in self.child_frequencies:
            if self.child_frequencies[(cl,cr)]>0:
                for c in [cl,cr]:
                    if not c in d:
                        d[c]=0
                        l.append(c)
                        c.aux_liste(l,d)

    def liste(self):
        d=dict()
        l=[self]
        self.aux_liste(l,d)
        return l

    def reverse_post_order_traversal(self):
        def order_clade(clade):
            return(len(clade.clade_leaves))
        l=self.liste()
        l.sort(key=order_clade, reverse=True)
        #ordre donné par la taille des clades est bien un ordre qui empeche les fils d'apparaitre avant leur parent
        self.reverse_post_order=l
        return l

    def pruning(self):
        for c in self.reverse_post_order:
            if not c.is_leaf():
                m=max([c.child_frequencies[u] for u in c.child_frequencies])
                print(m,[c.child_frequencies[u] for u in c.child_frequencies])
                for u in c.child_frequencies:
                    if m > 0.1+c.child_frequencies[u]:
                        c.child_frequencies[u]=0
        self.reverse_post_order_traversal()



    def leaves_traversal(self):
        if self.reverse_post_order is None:
            self.reverse_post_order_traversal()
        l=[u for u in self.reverse_post_order if u.is_leaf()]
        self.leaves=l
        return l

    def leaves(self):
        if self.leaves is None:
            return self.leaves_traversal()
        else:
            return self.leaves

    #compute log of the frequencies
    def log_freq(self):
        for clade in self.reverse_post_order:
            log_d=dict()
            for k in clade.child_frequencies:
                log_d[k]=np.log(clade.child_frequencies[k])
            clade.log_child_frequencies=log_d

    def save_match(self):
        for clade in self.reverse_post_order:
            clade.constant_match=clade.match

    def name_to_tree(self):
        d=dict()
        for u in self.reverse_post_order:
            d[u.name]=u
        return d


##################$$


class Rec_event:
    def __init__(self):
        self.name=None
        self.lower=None
        self.lower_left=None
        self.lower_right=None
        self.upper=None
        self.upper_left_or_keeper_or_receiver=None
        self.upper_right_or_loser_or_donor=None
        self.upper_match=None
        self.upper_left_match=None
        self.upper_right_match=None

    def is_third_level(self):
        return self.upper_match is None

    def key(self):
        d=vars(self)
        l=[]
        for k in d:
            l.append(d[k])
        return tuple(l)

    def init_name(self,third_level=False):
        #if self.lower.root.tree_name!="gene.nwk" and self.upper.root.tree_name!="host.nwk" and self.upper.root.tree_name!=self.lower.root.tree_name:
        #    print(self.upper.root.tree_name,self.lower.root.tree_name)
        if self.upper.root.tree_name==self.lower.root.tree_name:
            self.upper = "FREE_LIVING"
        d=vars(self)
        if third_level:
            self.upper_match=self.upper.match
            if not self.upper_left_or_keeper_or_receiver is None:
                self.upper_left_match=self.upper_left_or_keeper_or_receiver.match
            if not self.upper_right_or_loser_or_donor is None:
                self.upper_right_match=self.upper_right_or_loser_or_donor.match
        for key in d:
            if (not (d[key] is None)) and not key=="name" and not d[key]=="FREE_LIVING":
                if type(d[key]) == list:
                    d[key]=tuple([u.name for u in d[key]])
                else:
                    d[key]=d[key].name



class Rec_scenario:
    def __init__(self):
        self.event_list=[]
        self.reconstructed_lower=Tree()
        self.reconstructed_lower.event_list=[]
        self.log_likelihood=0
        self.am_tree_to_reconstructed_tree=dict()

#l_scenario is a list of reconstructed_tree
class Rec_sol:
    def __init__(self,log_likelihood,event_list_by_fam,l_scenario):
        self.log_likelihood=log_likelihood
        self.event_list_by_fam=event_list_by_fam
        self.scenario_list=l_scenario
        self.upper_divided_sol=None
        self.log_likelihood_by_gene=None
        self.upper_scenario=None
        self.upper_log_likelihood=None
        self.log_likelihood_by_upper=None



