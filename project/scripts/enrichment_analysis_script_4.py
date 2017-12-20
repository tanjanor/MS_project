## Script for analysing the gene enrichment of a set of differentially expressed genes.
## Input: 

#Last edited: 20/12/2017

#input paths
path_category = "../data/openMS/results/category_loc.csv"
path_sign = "../data/openMS/statistics/sign_genes.csv"
path_iden = "../data/openMS/statistics/diffacto_q_values.csv"

#output path
path_sign_out = "../data/openMS/enrichment/sign_process.csv"
path_iden_out = "../data/openMS/enrichment/identified_process.csv"
path_process = "../data/openMS/enrichment/enriched_processes_3.csv"

import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

#import file with significantly expressed genes:
sign_df = pd.read_csv(path_sign, header=0, index_col=0, usecols=["protein","f_value","p_value","q_value"])
cg_df = pd.read_csv(path_category, header=0, index_col=0, usecols=["GeneID","Process"])
iden_df = pd.read_csv(path_iden, header=0, index_col=0, usecols=["protein","f_value","p_value","q_value"])

###proteome###
#get counts for proteins in the protome for all process types
process_counts_all = cg_df["Process"].value_counts()
process_all = process_counts_all.index.tolist()
counts_all = process_counts_all.values.tolist()
#counts_df = pd.DataFrame(zip(process_all, counts_all))
#counts_df = pd.DataFrame({'process':process_all,'counts_proteome':counts_all})
counts_dict = dict(zip(process_all, counts_all))

###ANOVA: significant###
#get counts for sign proteins
process_list_sign = []
#adding processes for each gene to dataframe
for index in sign_df.index.get_values():		#index: proteinID	
	process_list_sign.append(cg_df.loc[index]["Process"])
sign_df = sign_df.assign(process = process_list_sign)
#output df with processes
sign_df.to_csv(path_sign_out)
#counting number for each process
process_counts = sign_df["process"].value_counts()
process_sign = process_counts.index.tolist()
counts_sign = process_counts.values.tolist()
#add to dictionary with all_counts
for i in range(len(process_sign)):
	counts_dict[process_sign[i]] = [counts_dict[process_sign[i]],counts_sign[i]]

###Identified###
#get counts for iden proteins
process_list_iden = []
##adding process to all genes
for index in iden_df.index.get_values():		#index: proteinID	
	process_list_iden.append(cg_df.loc[index]["Process"])
iden_df = iden_df.assign(process = process_list_iden)
#output df with processes
iden_df.to_csv(path_iden_out)
##counting number for each process
process_counts = iden_df["process"].value_counts()
process_iden = process_counts.index.tolist()
counts_iden = process_counts.values.tolist()
#add to dictionary with all_counts
for i in range(len(process_iden)):
	counts_dict[process_iden[i]].append(counts_iden[i])	
	
	
###Statistics: Hypergeometric test###
M= len(cg_df)		#total proteins in proteome
N= len(sign_df)		#total significant genes
M_2 = len(iden_df)	#total identified genes
count_sign_plot =[]
count_iden_plot =[]
count_genome_plot =[]
process_plot =[]
for process, counts in counts_dict.items():
	n= counts[0]	#process counts in whole genome
	x= counts[1]	#process counts in significant
	n_2=counts[2]	#process counts in identified
	pval = stats.hypergeom.sf(x-1, M, n, N)
	pval_2 = stats.hypergeom.sf(x-1, M_2, n_2, N)
	counts_dict[process] = counts_dict[process] + [pval, pval_2] 
	
	count_genome_plot.append(n)
	count_sign_plot.append(x)
	count_iden_plot.append(n_2)
	process_plot.append(process)

#output process counts dataframe
df_counts = pd.DataFrame(counts_dict)
df_counts = df_counts.transpose()
df_counts.index.names = ["process"]
df_counts.columns =["counts_genome","counts_significant","counts_identified","p-value(sign vs genome)","p-value(sign vs identified)"]
df_counts.to_csv(path_process)




###Plots###
num = len(process_plot)
ind = np.arange(num)  # the x locations for the groups
width = 0.3       # the width of the bars
fig, ax = plt.subplots(figsize=(13,6)) 
rects1 = ax.barh(ind, count_genome_plot, width, color='teal')
rects2 = ax.barh(ind + width, count_iden_plot, width, color='green')
rects3 = ax.barh(ind + width*2, count_sign_plot, width, color='orange')

# add some text for labels, title and axes ticks
ax.set_yticklabels(process_plot)
ax.set_xlabel('Number of genes annotated')
ax.set_title('Enrichment of Categories')
ax.set_yticks(ind + width/2)

#save plot
ax.legend((rects1[0], rects2[0], rects3[0]), ('Proteome', 'Identified', 'ANOVA q<0.05'))	
plt.subplots_adjust(left=0.4)
plt.savefig("../data/openMS/enrichment/enrichment_process_plot_3.png")
plt.close(fig)


#Normalized plot
count_genome_plot_norm = [x / M for x in count_genome_plot]
count_iden_plot_norm = [x / M_2 for x in count_iden_plot]
count_sign_plot_norm = [x / N for x in count_sign_plot]
num = len(process_plot)
ind = np.arange(num)  # the x locations for the groups
width = 0.3       # the width of the bars
fig, ax = plt.subplots(figsize=(13,6)) 
rects1 = ax.barh(ind, count_genome_plot_norm, width, color='teal')
rects2 = ax.barh(ind + width, count_iden_plot_norm, width, color='green')
rects3 = ax.barh(ind + width*2, count_sign_plot_norm, width, color='orange')

# add some text for labels, title and axes ticks
ax.set_yticklabels(process_plot)
ax.set_xlabel('Normalized number of genes annotated')
ax.set_title('Normalized Enrichment of Categories')
ax.set_yticks(ind + width / 2)

#save plot
ax.legend((rects1[0], rects2[0], rects3[0]), ('Proteome', 'Identified', 'ANOVA q<0.05'))
plt.subplots_adjust(left=0.4)
plt.savefig("../data/openMS/enrichment/enrichment_process_plot_norm_3.png")
plt.close(fig)


