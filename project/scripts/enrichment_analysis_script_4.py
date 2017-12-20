## Script for analysing the gene enrichment of a set of differentially expressed genes.
## Input: 

#Last edited: 20/12/2017

#input paths
path_category = "../data/openMS/results/category_loc.csv"
path_sign = "../data/openMS/statistics/sign_genes.csv"
path_iden = "../data/openMS/statistics/diffacto_q_values.csv"

#output path
path_process = "../data/openMS/enrichment/enriched_processes_3.csv"

import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

#import file with significantly expressed genes:
sign_df = pd.read_csv(path_sign, header=0, index_col=0, usecols=["protein","f_value","p_value","q_value"])
iden_df = pd.read_csv(path_iden, header=0, index_col=0, usecols=["protein","f_value","p_value","q_value"])
cg_df = pd.read_csv(path_category, header=0, index_col=0, usecols=["GeneID","Process"])

#get counts for genes in the protome for all process types
process_counts_all = cg_df["Process"].value_counts()
process_all = process_counts_all.index.tolist()
counts_all = process_counts_all.values.tolist()
counts_dict = dict(zip(process_all, counts_all))
#get counts for sign proteins
process_list = []
##adding process to all genes
for index, row in sign_df.iterrows():	#row: "f_value","p_value","q_value"	#index: proteinID	
	process_list.append(cg_df.loc[index]["Process"])
sign_df = sign_df.assign(process = process_list)
##counting number for each process
process_counts = sign_df["process"].value_counts()
process_sign = process_counts.index.tolist()
counts_sign = process_counts.values.tolist()
#add to dictionary with all_counts
for i in range(len(process_sign)):
	counts_dict[process_sign[i]] = [counts_dict[process_sign[i]],counts_sign[i]]

#Statistics: Hypergeometric test
##What is the estimated counts for each process type based upon their frequencies in the genome?
##Is the seen counts significantly different from the expected?	
M= len(cg_df)		#total genes
n= len(sign_df)		#total significant genes
out = open(path_process,"w")
out.write("Process"+"\t"+"Genome_count"+"\t"+"significant_count"+"\t"+"p-value"+"\n")
count_sign_plot =[]
count_genome_plot =[]
process_plot =[]
for process, counts in counts_dict.items():
	x= counts[1]	#process counts in significant
	N= counts[0]	#process counts in whole genome
	pval = stats.hypergeom.sf(x-1, M, n, N)
	out.write(process+"\t"+str(N)+"\t"+str(x)+"\t"+str(pval)+"\n")
	#print(pval,process)
	count_sign_plot.append(x)
	count_genome_plot.append(N)
	process_plot.append(process)
out.close()

#Plotting
num = len(process_plot)
ind = np.arange(num)  # the x locations for the groups
width = 0.4       # the width of the bars
w=13
h=6
fig, ax = plt.subplots(figsize=(w,h)) 
rects1 = ax.barh(ind, count_genome_plot, width, color='teal')
rects2 = ax.barh(ind + width, count_sign_plot, width, color='orange')

# add some text for labels, title and axes ticks
ax.set_yticklabels(process_plot)
ax.set_xlabel('Number of genes annotated')
ax.set_title('Enrichment of Categories')
ax.set_yticks(ind + width/2)


ax.legend((rects1[0], rects2[0]), ('Proteome', 'Differentially expressed'))	
#rcParams.update({'figure.autolayout':True})
plt.subplots_adjust(left=0.4)
#plt.tight_layout()
plt.savefig("../data/openMS/enrichment/enrichment_process_plot_3.png")
plt.close(fig)


#Normalized plot:
count_genome_plot_norm = [x / M for x in count_genome_plot]
count_sign_plot_norm = [x / n for x in count_sign_plot]
num = len(process_plot)
ind = np.arange(num)  # the x locations for the groups
width = 0.4       # the width of the bars
w=13
h=6
fig, ax = plt.subplots(figsize=(w,h)) 
rects1 = ax.barh(ind, count_genome_plot_norm, width, color='teal')
rects2 = ax.barh(ind + width, count_sign_plot_norm, width, color='orange')

# add some text for labels, title and axes ticks
ax.set_yticklabels(process_plot)
ax.set_xlabel('Normalized number of genes annotated')
ax.set_title('Normalized Enrichment of Categories')
ax.set_yticks(ind + width / 2)


ax.legend((rects1[0], rects2[0]), ('Proteome', 'Differentially expressed'))	
#rcParams.update({'figure.autolayout':True})
plt.subplots_adjust(left=0.4)
#plt.tight_layout()
plt.savefig("../data/openMS/enrichment/enrichment_process_plot_norm_3.png")
plt.close(fig)







"""
#Plotting
num = len(process_list)



ind = np.arange(num)  # the x locations for the groups
width = 0.35       # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind, count_genome_plot, width, color='r')
rects2 = ax.bar(ind + width, count_sign_plot, width, color='y')

# add some text for labels, title and axes ticks
ax.set_ylabel('Number of genes annotated')
ax.set_title('Categories Enrichment')
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(process_plot)

ax.legend((rects1[0], rects2[0]), ('Proteome', 'Differentially expressed'))	







x= process_counts["Energy metabolism"]
M= len(cg_df)
n= len(sign_df)
N= process_counts_all["Energy metabolism"]
pval = stats.hypergeom.sf(x-1, M, n, N)
print(pval)

print(process_sign)
print(counts_sign)
print(processes_sign)
print(counts_sign)
#print(sign_df)
#print(cg_df)
for i in range(len(process_counts_all)):
	print(process_counts_all[i])
print()
print(process_counts_all)
#prb = stats.hypergeom.pmf(x, M, n, N)
print(prb)"""
