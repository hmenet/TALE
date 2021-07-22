
import random as rd
import arbre
import time
import matplotlib.pyplot as plt

from reconciliation import reconciliation
from rates_inference import gene_rates_ml
from out_recphyloxml import save_recphyloxml_from_l_event
from read_input import read_input
from arbre import save_tree, Tree


def origination_proba(l_event_by_family):
    orig_prob=dict()
    for l_event in l_event_by_family:
        for u in l_event:
            if u[2]==0:
                if not u[1] in orig_prob:
                    orig_prob[u[1]]=0
                orig_prob[u[1]]+=l_event[u]
    total=sum([orig_prob[u] for u in orig_prob])
    for u in orig_prob:
        orig_prob[u]/=total
    return orig_prob

# on calcul les fréquences
# il faut les fréquences d'arrêt dans une feuille ?
def events_frequencies(l_event_by_family, c_match_list, event_types=["S","T","TL", "SL", "E"]):
    n_by_branch=dict()
    i_clade=0
    for l_event in l_event_by_family:
        c_match=c_match_list[i_clade]
        for u in l_event:
            if not u[1] in n_by_branch:
                n_by_branch[u[1]]=dict()
            E=u[0]
            if E in event_types:
                if "T" in E or "L" in E:
                    E=(u[0], u[3])
                if not E in n_by_branch[u[1]]:
                    n_by_branch[u[1]][E]=0
                n_by_branch[u[1]][E]+=l_event[u]

            for e,v in [(u[3], u[4]),(u[5], u[6])]:
                if e.isLeaf() and v in c_match and e in c_match[v]:
                    E="E"#E like End
                    if not e in n_by_branch:
                        n_by_branch[e]=dict()
                    if not E in n_by_branch[e]:
                        n_by_branch[e][E]=0
                    n_by_branch[e][E]+=l_event[u]
        i_clade+=1
    freq_by_branch=dict()
    for b in n_by_branch:
        freq_by_branch[b]=dict()
        total=sum([n_by_branch[b][E] for E in n_by_branch[b]])
        for E in n_by_branch[b]:
            freq_by_branch[b][E]=n_by_branch[b][E]/total
    return freq_by_branch

def generate_tree(freq_by_branch, orig_prob, plus_S=0):
    #origination
    orig_branch=(rd.choices([u for u in orig_prob], weights=[orig_prob[u] for u in orig_prob]))[0]
    r=dict()
    tree=arbre.Tree()
    r[tree]=[orig_branch]
    to_see=[tree]
    while len(to_see)>0:
        u=to_see.pop()
        #e is the last host of u
        e=r[u][-1]
        weights_to_use=[]
        for E in freq_by_branch[e]:
            E2=E
            if type(E)==tuple:
                E2,h=E
            if E2=="SL" or E2=="S" or E2=="E":
                weights_to_use.append(freq_by_branch[e][E]+plus_S)
            else:
                weights_to_use.append(freq_by_branch[e][E])
        #E=rd.choices([E for E in freq_by_branch[e]], weights=[freq_by_branch[e][E] for E in freq_by_branch[e]])[0]
        E=rd.choices([E for E in freq_by_branch[e]], weights=weights_to_use)[0]

        if type(E)==tuple:
            E,h=E
        if E=="S":
            u.undated_birth()
            r[u.right]=[e.right]
            r[u.left]=[e.left]
            to_see.append(u.left)
            to_see.append(u.right)
        if E=="D":
            u.undated_birth()
            r[u.right]=[e]
            r[u.left]=[e]
            to_see.append(u.left)
            to_see.append(u.right)
        if E=="SL":
            r[u].append(h)
            to_see.append(u)
        if E=="T":
            u.undated_birth()
            r[u.left]=[e]
            r[u.right]=[h]
            to_see.append(u.left)
            to_see.append(u.right)
        if E=="TL":
            r[u].append(h)
            to_see.append(u)
        #if E=="E"#on s'arrête, rien à faire
    for u in tree.leaves():
        u.match=r[u][-1]
    return tree

#some sort of SPR simulation but with DTL events, may return a smaller tree, but unicopy
def generate_tree_SPR(freq_by_branch, orig_prob):
    #origination
    orig_branch=(rd.choices([u for u in orig_prob], weights=[orig_prob[u] for u in orig_prob]))[0]
    r=dict()
    tree=arbre.Tree()
    tree.time=0.1
    r[tree]=[orig_branch]
    to_see=[tree]
    upper_seen=dict()
    upper_seen[orig_branch]=True
    while len(to_see)>0:
        u=to_see.pop(rd.randrange(len(to_see)))
        #u=to_see.pop()
        #e is the last host of u
        e=r[u][-1]
        E=rd.choices([E for E in freq_by_branch[e]], weights=[freq_by_branch[e][E] for E in freq_by_branch[e]])[0]
        if type(E)==tuple:
            E,h=E
        if E=="S":
            if e.right in upper_seen:
                if not e.left in upper_seen:
                    E=="SL"
                    h=e.left
            if e.left in upper_seen:
                if not e.right in upper_seen:
                    E=="SL"
                    h=e.right
            if not(e.left in upper_seen or e.right in upper_seen):
                u.birth(0.5)
                u.left.time=0.5
                u.right.time=0.5
                r[u.right]=[e.right]
                r[u.left]=[e.left]
                to_see.append(u.left)
                to_see.append(u.right)
                upper_seen[e.right]=True
                upper_seen[e.left]=True
        if E=="T":
            if not h in upper_seen:
                u.birth(0.5)
                u.left.time=0.5
                u.right.time=0.5
                r[u.left]=[e]
                r[u.right]=[h]
                to_see.append(u.left)
                to_see.append(u.right)
                upper_seen[h]=True
        if E=="SL":
            if not h in upper_seen:
                r[u].append(h)
                to_see.append(u)
        if E=="TL":
            if not h in upper_seen:
                r[u].append(h)
                to_see.append(u)
                upper_seen[h]=True
        if E=="E":
            u.time=1
    for u in tree.leaves():
        u.match=r[u][-1]
    tree.prune_tree()
    return tree

#SPR guided by rates
def generate_tree_SPR2(host,freq_by_branch, orig_prob, plus_S=0):
    tree=host.copy()
    freq_by_branch_with_name=dict()
    for e in freq_by_branch:
        freq_by_branch_with_name[e.name]=dict()
        for E in freq_by_branch[e]:
            freq_by_branch_with_name[e.name][E]=freq_by_branch[e][E]
    freq_by_branch=freq_by_branch_with_name
    r=dict()
    name_to_tree_host=dict()
    for u in host.liste():
        name_to_tree_host[u.name]=u
        print(u.name)
    name_to_tree=dict()
    for u in tree.liste():
        print(u.name)
        name_to_tree[u.name]=u
    l_name=[e_name for e_name in freq_by_branch]
    rd.shuffle(l_name)
    for e_name in l_name:
        if e_name in [u.name for u in tree.liste()]:
            e=name_to_tree[e_name]
            weights_to_use=[]
            for E in freq_by_branch[e_name]:
                E2=E
                if type(E)==tuple:
                    E2,h=E
                if E2=="SL" or E2=="S" or E2=="E":
                    weights_to_use.append(freq_by_branch[e_name][E]+plus_S)
                else:
                    weights_to_use.append(freq_by_branch[e_name][E])
            E=rd.choices([E for E in freq_by_branch[e_name]], weights=weights_to_use)[0]
            if type(E)==tuple:
                E,h=E
            if E=="T" or E=="TL":
                print("once")
                print(e.name, e.isLeaf())
                if e.isLeaf():
                    e.left=arbre.Tree()
                    e.left.name=e.name
                    e.left.parent=e
                    e.right=h
                else:
                    tmpl=e.left
                    print(type(e.left))
                    tmpr=e.right
                    e.left=arbre.Tree()
                    e.left.left=tmpl
                    e.left.right=tmpl
                    tmpl.parent=e.left
                    tmpr.parent=e.left
                    e.right=h
                print(h.parent.name)
                print(h.name)
                print(h.parent.left.name, h.parent.right.name)
                if h.parent.left==h:
                    sibling=h.parent.right
                else:
                    sibling=h.parent.left
                if sibling.isLeaf():
                    sibling.parent.name=sibling
                    sibling.parent.left=None
                    sibling.parent.right=None
                else:
                    sibling.parent.left=sibling.left
                    sibling.parent.right=sibling.right
                    sibling.parent.left.parent=sibling.parent
                    sibling.parent.right.parent=sibling.parent
                h.parent=e

    for u in tree.leaves():
        u.match=name_to_tree_host[u.name]
    return tree


def test_number_leaf(c_match_list, intermediate_list):
    upper_to_inter=dict()
    for intermediate in intermediate_list:
        upper_seen=dict()
        for u in intermediate.leaves():
            if not u.match in upper_to_inter:
                upper_to_inter[u.match]=[]
            upper_to_inter[u.match].append(u)
    for c_match in c_match_list:
        for c in c_match:
            if type(c_match[c])==list:
                for u in c_match[c]:
                    if not u in upper_to_inter:
                        return False
            else:
                if not u in upper_to_inter:
                    return False
    return True


def matching_lower_to_intermediate(c_match_list, intermediate_list, host_list):
    upper_to_inter=dict()
    for intermediate in intermediate_list:
        for u in intermediate.leaves():
            if intermediate in host_list:
                if not u in upper_to_inter:
                    upper_to_inter[u]=[]
                upper_to_inter[u].append(u)
            else:
                if not u.match in upper_to_inter:
                    upper_to_inter[u.match]=[]
                upper_to_inter[u.match].append(u)
    new_c_match_list=[]
    for c_match in c_match_list:
        new_c_match=dict()
        for c in c_match:
            new_c_match[c]=[]
            if type(c_match[c])==list:
                for u in c_match[c]:
                    new_c_match[c]=new_c_match[c]+upper_to_inter[u]
            else:
                new_c_match[c]=upper_to_inter[u]
        new_c_match_list.append(new_c_match)
    return new_c_match_list

def prune_some_leaf(intermediate):
    seen_leaf=[]
    leaves=intermediate.leaves()
    rd.shuffle(leaves)
    for u in leaves:
        if u.name in seen_leaf:
            u.time=0.5
        else:
            u.time=1
            seen_leaf.append(u.name)
    intermediate.prune_tree()


####




def rec(symbiont_list, clades_data_list, c_match_list, n_sample=100, n_steps=5, n_rec_sample_rates=100, best_rec=False, n_recphyloxml=0, multiprocess=False, multiprocess_fam=False, out_file="output/rec"):
    parasite_post_order=[]
    for symbiont in symbiont_list:
        parasite_post_order+=symbiont.post_order_traversal()
    t1=time.perf_counter()
    rates=gene_rates_ml(parasite_post_order,clades_data_list,c_match_list, n_steps, init_rates_g=[0.01,0.01,0.01], n_rec_sample=n_rec_sample_rates, multi_process=multiprocess, multi_process_family=multiprocess_fam)
    print(rates)
    cmpt_time=time.perf_counter()-t1
    print("Rates estimated in ", cmpt_time, " s")
    print("Rates", rates)
    t1=time.perf_counter()
    if n_sample>0:
        likelihood,l_event_gene, l_scenarios=reconciliation(parasite_post_order, clades_data_list, c_match_list, rates, sample=True, n_sample=n_sample, best=best_rec, n_recphyloxml=n_recphyloxml, multi_process=multiprocess, multi_process_family=multiprocess_fam)
    else:
        likelihood=reconciliation(parasite_post_order, clades_data_list, c_match_list, rates, sample=False, n_sample=n_sample, best=best_rec, n_recphyloxml=n_recphyloxml, multi_process=multiprocess,multi_process_family=multiprocess_fam)
    cmpt_time=time.perf_counter()-t1
    print("Reconciliation ended in ", cmpt_time, " s")
    print("Log Likelihood: ", likelihood)
    print("Rates: ", rates)

    for i_recphyloxml in range(min(n_sample, n_recphyloxml)):
        scenario_by_family=l_scenarios[i_recphyloxml]
        if best_rec and i_recphyloxml==0:
            out_file_name=out_file+str(i_recphyloxml)+"_best"+".recphyloxml"
        else:
            out_file_name=out_file+str(i_recphyloxml)+".recphyloxml"
        save_recphyloxml_from_l_event(symbiont_list, scenario_by_family, out_file_name, c_match_list=c_match_list, clade_data_list=clades_data_list, clade=True)

    if n_sample>0:
        return likelihood, l_event_gene
    else:
        return likelihood



upper_dir="/home/hmenet/Documents/rewrite_ALE/example_data/ex_pylori/pop_pylori_1branch_2comp_same/"
lower_dir="/home/hmenet/Documents/rewrite_ALE/example_data/ex_pylori/genes_sub2/"
leaf_matching_file="/home/hmenet/Documents/rewrite_ALE/example_data/ex_pylori/matching_pop_genes_1branch2comp"


upper_dir="/home/hmenet/Documents/rewrite_ALE/example_data/ex_pylori/pop_pylori_1branch_2comp/"
lower_dir="/home/hmenet/Documents/rewrite_ALE/example_data/ex_pylori/genes_sub2/"
leaf_matching_file="/home/hmenet/Documents/rewrite_ALE/example_data/ex_pylori/matching_pop_genes_1branch2comp"



upper_dir="/home/hmenet/Documents/rewrite_ALE/example_data/ex_pylori/pop_pylori_1branch/"
lower_dir="/home/hmenet/Documents/rewrite_ALE/example_data/ex_pylori/genes_sub2/"
lower_dir="/home/hmenet/Documents/Stage_M2/These/script/h_pylori/genes_sub/"
leaf_matching_file="/home/hmenet/Documents/rewrite_ALE/example_data/ex_pylori/matching_pop_genes_1branch"

n_sample=100


symbiont_list, clades_data_list, c_match_list, gene_file_list=read_input(upper_dir, lower_dir, leaf_matching_file=leaf_matching_file)

likelihood, l_event_by_family=rec(symbiont_list, clades_data_list, c_match_list, n_sample=n_sample, n_steps=0, n_rec_sample_rates=n_sample, multiprocess=False, multiprocess_fam=True)
orig_prob=origination_proba(l_event_by_family)
freq_by_branch=events_frequencies(l_event_by_family, c_match_list, event_types=["S","T","TL", "SL", "E"])
for u in freq_by_branch:
    print(u.name, "\n")
    for v in freq_by_branch[u]:
        if len(v)>1:
            print(v[0], v[1].name, freq_by_branch[u][v])
        else:
            print(v, freq_by_branch[u][v])

    print("\n")
#for u in freq_by_branch:
#    print(u.name, u)

n_sample_comp=100000
n_inter=2
best_l=None
best_intermediate=None
list_l=[]
for i_sample in range(n_sample_comp):
    intermediate_list=[]
    total_leaves=0


    intermediate=generate_tree(freq_by_branch, orig_prob, plus_S=0.1)
    #intermediate=generate_tree_SPR2(symbiont_list[0], freq_by_branch, origination_proba, plus_S=1)
    intermediate_list.append(intermediate)
    intermediate0=intermediate
    total_leaves+=len(intermediate.leaves())
    already_seen=[]
    for u in intermediate.leaves():
        u.name=u.match.name
        if not u.name in already_seen:
            already_seen.append(u.name)
    if total_leaves<15 and len(already_seen)==8:
        prune_some_leaf(intermediate)
    total_leaves=len(intermediate.leaves())
    if len(already_seen)==8 and total_leaves==8:
        good=True
    else:
        good=False
    if good:


        intermediate0=symbiont_list[0]
        for i_sample2 in range(n_sample_comp//5):
            good=False
            total_leaves=0
            intermediate=generate_tree(freq_by_branch, orig_prob, plus_S=0.1)
            #intermediate=generate_tree_SPR2(symbiont_list[0], freq_by_branch, origination_proba, plus_S=1)
            total_leaves+=len(intermediate.leaves())
            already_seen=[]
            for u in intermediate.leaves():
                u.name=u.match.name
                if not u.name in already_seen:
                    already_seen.append(u.name)

            if total_leaves<15 and len(already_seen)==8:
                prune_some_leaf(intermediate)
                print("hey")
            total_leaves=len(intermediate.leaves())

            if len(already_seen)==8 and total_leaves==8:
                good=True
            if good:
                intermediate_list=[intermediate0, intermediate]
                #for symbiont in symbiont_list:
                #    intermediate_list.append(symbiont)
                new_c_match_list=matching_lower_to_intermediate(c_match_list, intermediate_list, symbiont_list)
                l=rec(intermediate_list, clades_data_list, new_c_match_list, n_sample=0, n_rec_sample_rates=n_sample, multiprocess=False, multiprocess_fam=True, n_steps=0)
                list_l.append(l)
                if best_l==None or l>best_l:
                    best_l=l
                    best_intermediate=intermediate_list
        print("\none done\n")
    if i_sample%1000 == 0:
        print(i_sample, "/", n_sample_comp)



""""
    #print(test_number_leaf(c_match_list, intermediate_list))
    #intermediate_list=list(symbiont_list)
    for symbiont in symbiont_list:
        intermediate_list.append(symbiont)

    #if test_number_leaf(c_match_list, intermediate_list) and len(intermediate_list[0].leaves()) in  [7,8,9,10,11,12] and len(intermediate_list[1].leaves()) in  [7,8,9,10,11,12]:
    if total_leaves==n_inter*8 and good:
        new_c_match_list=matching_lower_to_intermediate(c_match_list, intermediate_list, symbiont_list)
        print("hey", total_leaves)
        l=rec(intermediate_list, clades_data_list, new_c_match_list, n_sample=0, n_rec_sample_rates=n_sample, multiprocess=False, multiprocess_fam=True, n_steps=0)
        #l=rec(intermediate_list, clades_data_list, new_c_match_list, n_sample=1, n_rec_sample_rates=n_sample, multiprocess=False, multiprocess_fam=False, n_steps=0, n_recphyloxml=1, best_rec=True)
        list_l.append(l)
        if best_l==None or l>best_l:
            best_l=l
            best_intermediate=intermediate_list
    if i_sample%1000 == 0:
        print(i_sample, "/", n_sample_comp)
"""



print(len(list_l))
print(best_l)
plt.hist(list_l)
plt.show()
save_tree(best_intermediate[0], "/home/hmenet/Documents/rewrite_ALE/output/intermediate0")
save_tree(best_intermediate[1], "/home/hmenet/Documents/rewrite_ALE/output/intermediate1")

### script to create matching file for pylori

def create_matching_file_pylori():
    from read_tree import read_tree
    gene = read_tree("ex_pylori/genes_sub2/FAM000006.fasta_nuc.filtred_aligned.treefile")

    d_pop_geo_match=dict()
    for u in ["hpEastAsia","hpAsia2","hpSahul","hpNEAfrica","hpAfrica1","hpAfrica2","outgroup","hpEurope"]:
        d_pop_geo_match[u]=u
    for u in ["hspEAsia", "hpEAsia"]:
        d_pop_geo_match[u]="hpEastAsia"
    for u in ["hspMaori"]:
        d_pop_geo_match[u]="hpEastAsia"
    for u in ["hpAfrica","hpWAfrica","hspSAfrica","hspSAfrica1","hspWAfrica","hspWAfrica1"]:
        d_pop_geo_match[u]="hpAfrica1"


    f=open("ex_pylori/matching_pop_genes", "w")
    i=0
    for u in gene.leaves_unrooted():
        s=u.name
        pylori_pop=s[s.rfind("_")+1:]
        present_pop=d_pop_geo_match[pylori_pop]
        if present_pop=="hpEurope":
            f.write(s)
            f.write("\t")
            f.write("hpEurope_7")
            f.write("\t")
            f.write("pop")
            f.write("\n")
            f.write(s)
            f.write("\t")
            f.write("hpEurope_4")
            f.write("\t")
            f.write("pop")
        else:
            f.write(s)
            f.write("\t")
            f.write(present_pop)
            f.write("\t")
            f.write("pop")
        i+=1
        if i < len(gene.leaves_unrooted()):
            f.write("\n")
    f.close()

#create_matching_file_pylori()

def create_matching_file_pylori_1branch():
    from read_tree import read_tree
    gene = read_tree("ex_pylori/genes_sub2/FAM000006.fasta_nuc.filtred_aligned.treefile")

    d_pop_geo_match=dict()
    for u in ["hpEastAsia","hpAsia2","hpSahul","hpNEAfrica","hpAfrica1","hpAfrica2","outgroup","hpEurope"]:
        d_pop_geo_match[u]=u
    for u in ["hspEAsia", "hpEAsia"]:
        d_pop_geo_match[u]="hpEastAsia"
    for u in ["hspMaori"]:
        d_pop_geo_match[u]="hpEastAsia"
    for u in ["hpAfrica","hpWAfrica","hspSAfrica","hspSAfrica1","hspWAfrica","hspWAfrica1"]:
        d_pop_geo_match[u]="hpAfrica1"


    f=open("ex_pylori/matching_pop_genes_1branch", "w")
    i=0
    for u in gene.leaves_unrooted():
        s=u.name
        pylori_pop=s[s.rfind("_")+1:]
        present_pop=d_pop_geo_match[pylori_pop]
        if present_pop=="hpEurope":
            f.write(s)
            f.write("\t")
            f.write("hpEurope_4")
            f.write("\t")
            f.write("pop")
        else:
            f.write(s)
            f.write("\t")
            f.write(present_pop)
            f.write("\t")
            f.write("pop")
        i+=1
        if i < len(gene.leaves_unrooted()):
            f.write("\n")
    f.close()


#create_matching_file_pylori_1branch()


def create_matching_file_pylori_1branch2comp():
    from read_tree import read_tree
    gene = read_tree("ex_pylori/genes_sub2/FAM000006.fasta_nuc.filtred_aligned.treefile")

    d_pop_geo_match=dict()
    for u in ["hpEastAsia","hpAsia2","hpSahul","hpNEAfrica","hpAfrica1","hpAfrica2","outgroup","hpEurope"]:
        d_pop_geo_match[u]=u
    for u in ["hspEAsia", "hpEAsia"]:
        d_pop_geo_match[u]="hpEastAsia"
    for u in ["hspMaori"]:
        d_pop_geo_match[u]="hpEastAsia"
    for u in ["hpAfrica","hpWAfrica","hspSAfrica","hspSAfrica1","hspWAfrica","hspWAfrica1"]:
        d_pop_geo_match[u]="hpAfrica1"


    f=open("ex_pylori/matching_pop_genes_1branch2comp", "w")
    i=0
    for u in gene.leaves_unrooted():
        s=u.name
        pylori_pop=s[s.rfind("_")+1:]
        present_pop=d_pop_geo_match[pylori_pop]
        if present_pop=="hpEurope":
            f.write(s)
            f.write("\t")
            f.write("hpEurope_4")
            f.write("\t")
            f.write("pop")
        else:
            f.write(s)
            f.write("\t")
            f.write(present_pop)
            f.write("\t")
            f.write("pop")
        i+=1
        f.write("\n")
    i=0
    for u in gene.leaves_unrooted():
        s=u.name
        pylori_pop=s[s.rfind("_")+1:]
        present_pop=d_pop_geo_match[pylori_pop]
        if present_pop=="hpEurope":
            f.write(s)
            f.write("\t")
            f.write("hpEurope_4")
            f.write("\t")
            f.write("pop_Asie")
        else:
            f.write(s)
            f.write("\t")
            f.write(present_pop)
            f.write("\t")
            f.write("pop_Asie")
        i+=1
        if i < len(gene.leaves_unrooted()):
            f.write("\n")

    f.close()
#create_matching_file_pylori_1branch2comp()
