## Script for analysing the gene enrichment of a set of differentially expressed genes.
## Input: file with differentially expressed genes, file with identified genes in experiment
##		  file with process category-terms for all genes in species. 

#Last edited: 09/01/2018

#input paths
path_category = "../data/openMS/results/goterm.csv"
path_sign = "../data/openMS/statistics/sign_genes.csv"
path_iden = "../data/openMS/statistics/diffacto_q_values.csv"

#output path
path_sign_out = "../data/openMS/enrichment/new/sign_process.csv"
path_iden_out = "../data/openMS/enrichment/new/identified_process.csv"
path_process = "../data/openMS/enrichment/new/enriched_processes_new.csv"



import pandas as pd
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


#import file with significantly expressed genes:
sign_df = pd.read_csv(path_sign, header=0, index_col=0, usecols=["protein","f_value","p_value","q_value"])
cg_df = pd.read_csv(path_category, header=0, index_col=0, usecols=["gene_id","go_description","go_id"])
iden_df = pd.read_csv(path_iden, header=0, index_col=0, usecols=["protein","f_value","p_value","q_value"])

print(cg_df)
"""
###calc qvalue###
def get_q_value(pvals):
	q_sign = 0.05
	num = len(pvals)
	q_vals = []
	sign = []
	##starting from the lowest ranked p-value:
	prev_q = pvals[-1]
	q_vals.append(prev_q)
	if prev_q < q_sign:
		sign.append("Yes")
	else:
		sign.append("No")
	##minimum for those with higher rank
	for rank in range (num-1,0,-1):
		q_val = min(pvals[rank-1]*num/rank , prev_q)
		q_vals.append(q_val)
		prev_q = q_val
		##Find significant results
		if q_val < q_sign:
			sign.append("Yes")
		else:
			sign.append("No")
	return q_vals[::-1], sign[::-1]
"""
"""
##ignore
#Make a dictonary of GO-terms:
go_d = dict()
goid_d = dict()

for protein,go in zip(cg_df.index.tolist(),cg_df["go_description"]):
	if protein in go_d:
		go_d[protein] = go_d[protein] + go
	else:
		go_d[protein] = go
print(go_d)
"""


#new part for GO
GO_counts = cg_df["go_id"].value_counts()
print(GO_counts)
"""
###ANOVA: significant###
#get counts for sign proteins
process_list_sign = []
#adding processes for each gene to dataframe
for index in sign_df.index.get_values():		#index: proteinID	
	process_list_sign.append(cg_df.loc[index]["Process"])	#find process type
sign_df = sign_df.assign(process = process_list_sign)
#output df with processes
sign_df.to_csv(path_sign_out)
#counting number for each process
process_counts = sign_df["process"].value_counts()
process_sign = process_counts.index.tolist()
counts_sign = process_counts.values.tolist()
#add to dictionary with process counts
counts_df = pd.DataFrame(counts_sign,process_sign)
counts_df.index.names = ["process"]
counts_df.columns =["counts_significant"]
counts_df = counts_df.assign(counts_identified=([0]*17))
"""


"""
###Identified###
#get counts for iden proteins
process_list_iden = []
##adding process to all genes
for index in iden_df.index.get_values():		#index: proteinID	
	process_list_iden.append(cg_df.loc[index]["Process"])	#find process type
iden_df = iden_df.assign(process = process_list_iden)
#output df with processes
iden_df.to_csv(path_iden_out)
##counting number for each process
process_counts = iden_df["process"].value_counts()
for process,count in process_counts.items():
	counts_df.ix[process][1]=count
	
###Statistics: Hypergeometric test###
N= len(sign_df)		#total significant genes
M = len(iden_df)	#total identified genes
x=counts_df["counts_significant"]	#process counts in significant
n=counts_df["counts_identified"]	#process counts in identified
pvalsH = stats.hypergeom.sf(x-1, M, n, N, loc=0)	#perform hypergeometric test

pvalsL = []

#cdf wont accept my input as array/list/df. avoiding the struggle:
for xL,nL in zip(x,n):
	pvalL = stats.hypergeom.cdf(xL, M, nL, N, loc=0)
	pvalsL.append(pvalL)

counts_df = counts_df.assign(p_value_enrichment=pvalsH,p_value_depletion=pvalsL)


#Sort according to p-value_depletion and add q-values
counts_df_pval = counts_df.sort_values(by=["p_value_depletion"],ascending=True)
qvalsL, signL = get_q_value(counts_df_pval["p_value_depletion"])
counts_df_pval = counts_df_pval.assign(q_value_depletion=qvalsL)
counts_df_pval = counts_df_pval.assign(depletion_significant=signL)

#sort according to p-value enrichment and add q-values
counts_df_pval = counts_df_pval.sort_values(by=["p_value_enrichment"],ascending=True)
qvalsH, signH = get_q_value(counts_df_pval["p_value_enrichment"])
counts_df_pval = counts_df_pval.assign(q_value_enrichment=qvalsH)
counts_df_pval = counts_df_pval.assign(enrichment_significant=signH)

#output process counts dataframe	
counts_df_pval.to_csv(path_process)





###Plots########################################################################################################
counts_df_sign = counts_df.sort_values(by=["counts_significant"],ascending=True) #Sort according to #significant
num= len(counts_df)
fig, ax = plt.subplots(figsize=(13,7))
ind = np.arange(num)  # the x locations for the groups
width = 0.3       # the width of the bars
rects1 = ax.barh(ind, counts_df_sign["counts_identified"], width, color='green')
rects2 = ax.barh(ind + width, counts_df_sign["counts_significant"], width, color='orange')

# add some text for labels, title and axes ticks
ax.set_yticklabels(counts_df_sign.index)
ax.set_xlabel('Number of genes annotated')
ax.set_title('Enrichment of Categories')
ax.set_yticks(ind + width/2)

#save plot
ax.legend((rects1[0], rects2[0]), ('Identified', 'ANOVA q<0.05'))	
plt.subplots_adjust(left=0.4)
plt.savefig("../data/openMS/enrichment/new/enrichment_process_plot_new.png")
plt.close(fig)

###normalized plot###
norm=[x/y for x, y in zip( counts_df_pval["counts_significant"], counts_df_pval["counts_identified"])]
counts_df_pval = counts_df_pval.assign(significant_normalized=norm)
counts_df_norm = counts_df_pval.sort_values(by=["significant_normalized"],ascending=True) #sort according to norm-value
fig, ax = plt.subplots(figsize=(10,6))
rects = ax.barh(ind, counts_df_norm["significant_normalized"], width, color='#335B8E')
for i, bar in enumerate(rects):
	if counts_df_norm['enrichment_significant'][i] == 'Yes':
			bar.set_color("#B7DBDB")
# Shrink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

# Put a legend below current axis
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
#         fancybox=True, shadow=True, ncol=2)

# add some text for labels, title and axes ticks
ax.set_yticklabels(counts_df_norm.index)
ax.set_xlabel('Proportion of genes annotated')
ax.set_title('Normalised Enrichment of Categories')
ax.set_yticks(ind + width/2)

#save plot
colour_blue = matplotlib.patches.Patch(color='#335B8E', label='Nonsignificant')
colour_yel = matplotlib.patches.Patch(color="#B7DBDB", label='Significant enrichment')
#plt.legend(handles=[colour_blue,colour_yel], loc=8)
plt.legend(handles=[colour_blue,colour_yel], loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)
plt.subplots_adjust(left=0.5, bottom=0.2)
plt.savefig("../data/openMS/enrichment/new/enrichment_process_plot_new_norm.png")
plt.close(fig)
"""
"""
###normalized plot###V2
norm=[x/M for x in (counts_df_pval["counts_identified"])]
counts_df_pval = counts_df_pval.assign(significant_normalized=norm)
counts_df_norm = counts_df_pval.sort_values(by=["significant_normalized"],ascending=True) #sort according to norm-value
fig, ax = plt.subplots(figsize=(13,6))
rects = ax.barh(ind, counts_df_norm["significant_normalized"], width, color='blue')
for i, bar in enumerate(rects):
	if counts_df_norm['enrichment_significant'][i] == 'Yes':
			bar.set_color("#edd012")

# add some text for labels, title and axes ticks
ax.set_yticklabels(counts_df_norm.index)
ax.set_xlabel('Proportion of genes annotated')
ax.set_title('Normalised Enrichment of Categories')
ax.set_yticks(ind + width/2)

#save plot
colour_blue = matplotlib.patches.Patch(color='blue', label='Nonsignificant')
colour_yel = matplotlib.patches.Patch(color="#edd012", label='Significant enrichment')
plt.legend(handles=[colour_blue,colour_yel])
plt.subplots_adjust(left=0.4)
plt.savefig("../data/openMS/enrichment/new/enrichment_process_plot_new_norm_V2.png")
plt.close(fig)
"""