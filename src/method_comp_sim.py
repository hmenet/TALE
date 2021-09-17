#test method on simulation

from read_output import read_output_file_list
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_transfer_files(input_file_list):
    l_t_list=[]
    for input_file in input_file_list:
        f=open(input_file, "r")
        s=f.read()
        f.close()
        l_t=[]
        s1=s.split(sep="\n")
        for u in s1:
            s2=u.split(sep="\t")
            if len(s2)>1:
                l_t.append((s2[0],s2[1]))
        l_t_list.append(l_t)
    return l_t_list


def transfer_list_from_events(l_events_by_fam):
    l_t_list=[]
    for l_e in l_events_by_fam:
        l_t=dict()
        for u in l_e:
            if "T" in u[0]:
                if (u[1], u[2]) in l_t:
                    l_t[(u[1], u[2])]+=l_e[u]
                else:
                    l_t[(u[1], u[2])]=l_e[u]
        l_t_list.append(l_t)
    return l_t_list



def test_transfer(l_t_list,real_l_t_list):
    precision_list=[]
    recall_list=[]
    n_list=[]
    real_n_list=[]
    for i_gene in range(len(l_t_list)):
        l_t=l_t_list[i_gene]
        #real_l_t=l_t_list[i_gene]
        real_l_t=real_l_t_list[i_gene]
        if len(real_l_t)>0:
            tp=0#true positive
            for u in l_t:
                if u in real_l_t:
                    tp+=l_t[u]

            fp=sum([l_t[u] for u in l_t]) - tp #false positive
            fn=len(real_l_t)-tp #false negative
            print(sum([l_t[u] for u in l_t]), len(real_l_t))
            n_list.append(sum([l_t[u] for u in l_t]))
            real_n_list.append(len(real_l_t))
            if tp + fp == 0:
                precision=1
            else:
                precision=tp/(tp+fp)
            #si tp = 0, dur de mettre à 0 si ne trouve pas 0, alors que ce n'est pas un entier qu'on calcul
            #quoi mettre ?
            if tp+fn==0:
                recall=1
            else:
                recall=tp/(tp+fn)
            precision_list.append(precision)
            recall_list.append(recall)
    return precision_list, recall_list, n_list, real_n_list





heuristic_list=["2l","dec","MC"]
#heuristic_list=["2l","dec"]

path_dir="/home/hmenet/Documents/Stage_M2/These/script/simulation/output2/"
path_dir="/home/hmenet/Documents/Stage_M2/These/script/simulation/output/"
sim_number="_high1"
sim_number="1"

n_gene=5

path_dir_sim="/home/hmenet/Documents/Stage_M2/These/script/simulation/test_270721_sagephy/sim_290621_"


list_info=[]
#for sim_type in ["low","med","high"]:
for sim_type in ["highmed","highlow","low","med","high"]:
    for i_sim in range(1,6):
        #sim_type="med"
        #sim_type="high"
        #sim_type="low"
        sim_number=sim_type+str(i_sim)

        transfer_file_list=[path_dir_sim+sim_number+"/transfer_list/gene"+str(i)+".pruned.tree" for i in range(1,n_gene+1)]

        real_l_t_list=read_transfer_files(transfer_file_list)

        for heuristic in heuristic_list:
            print(heuristic)
            l_event_aggregate,l_events_by_fam=read_output_file_list([path_dir+sim_type+"/freq_"+sim_number+heuristic+"_"+"gene"+str(i)+".pruned.tree.nwk" for i in range(1,n_gene+1)])

            """
            gene_likelihood=[]
            for i in range(0,n_gene):
                log_file=path_dir+sim_type+"/freq_"+sim_number+heuristic+"_"+str(i)+"_likelihood"
                f=open(log_file,mode="r")
                s=f.read()
                log_l=float(s)
                gene_likelihood.append(log_l)
            """
            l_t_list=transfer_list_from_events(l_events_by_fam)
            precision, recall, n_list, real_n_list=test_transfer(l_t_list,real_l_t_list)
            print("precision",precision)
            print("recall",recall)
            for i_gene in range(len(precision)):
                list_info.append([precision[i_gene],recall[i_gene],n_list[i_gene],real_n_list[i_gene],heuristic,sim_type,sim_number])#,gene_likelihood[i_gene-1]])

list_columns=["precision","recall","inferred_n_transfer","sim_n_transfer","heuristic","sim_type","sim_number"]#,"log_likelihood"]
df=pd.DataFrame(list_info,columns=list_columns)


print(df.groupby("heuristic").mean())
#print(df.groupby("heuristic").precision.mean())

dfh=df.groupby("heuristic")

plt.figure()
for name, group in dfh:
    plt.scatter(group["recall"],group["precision"],label=name)

#df.groupby("heuristic").plot.scatter(x="precision",y="recall")

plt.legend()
plt.show()

df.boxplot(column="precision", by="heuristic")
plt.show()

df.boxplot(column="recall", by="heuristic")
plt.show()

#df.boxplot(column="log_likelihood", by="heuristic")
#plt.show()


print(df.groupby("sim_type").mean())

plt.figure()

for name, group in dfh:
    print(name)
    print(group.groupby("sim_type").mean())
    #plt.scatter(group["recall"],group["precision"],label=name)
    print("\n")

df.groupby("sim_type").boxplot(column="precision", by="heuristic")
plt.title("precision")
plt.show()

df.groupby("sim_type").boxplot(column="recall", by="heuristic")
plt.title("recall")
plt.show()

#df.groupby("sim_type").boxplot(column="log_likelihood", by="heuristic")
#plt.title("log likelihood")
#plt.show()


#plt.figure()
#for name, group in dfh:
#    plt.scatter(group["precision"],group["log_likelihood"],label=name)

#df.groupby("heuristic").plot.scatter(x="precision",y="recall")

#plt.legend()
#plt.show()

def test_pairwise(h1,h2):
    l2l=list(df[df.heuristic==h1].precision)
    ldec=list(df[df.heuristic==h2].precision)
    l3=[i - j for (i,j) in zip(l2l,ldec)]
    plt.figure()
    plt.hist(l3,bins=20)
    plt.show()
    print(len([i for i in l3 if i>0]))
    print(len([i for i in l3 if i<0]))
    print(len([i for i in l3 if i==0]))

    l2l=list(df[df.heuristic==h1].recall)
    ldec=list(df[df.heuristic==h2].recall)
    l3=[i - j for (i,j) in zip(l2l,ldec)]
    plt.figure()
    plt.hist(l3,bins=20)
    plt.show()
    print(len([i for i in l3 if i>0]))
    print(len([i for i in l3 if i<0]))
    print(len([i for i in l3 if i==0]))

test_pairwise("2l","dec")
test_pairwise("2l","MC")
test_pairwise("dec","MC")



#df.to_csv("/home/hmenet/Documents/Stage_M2/These/write/papier/figures/simulation/simulation_with_log.csv")

"""
plt.figure()
plt.xticks(fontsize=16)
plt.boxplot(y_precision_list)
plt.ylabel("Precision", fontsize=16)
plt.title("True positive detected HGTs / All detected HGTs", fontsize=18)
plt.show()

plt.figure()
plt.xticks(fontsize=16)
plt.boxplot(x_recall_list)
plt.ylabel("Recall", fontsize=16)
plt.title("True positive detected HGTs / All true HGTs", fontsize=18)
plt.show()
"""
