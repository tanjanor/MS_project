##Script to do an ANOVA on the output from openMS data
##to find differentially expressed genes

#Last edited: 19/12/2017


#Input paths
path_res = "../data/openMS/results/diffacto_weightedsum.csv"


#output paths
outpath_q = "../data/openMS/statistics/diffacto_q_values.csv"
outpath_sign = "../data/openMS/statistics/sign_genes.csv"

#q-value significance level
q_sign = 0.05

import pandas as pd
from scipy import stats

#import openMS file:
res_df = pd.read_csv(path_res, header=0, index_col=1)

#Remove first row in file (numbering not needed)
res_df.drop(res_df.columns[[0]], axis=1, inplace=True)

#statistics: oneway ANOVA
##tests the null hypothesis that two or more groups have the same population mean.
##scipy.stats.f_oneway(sample 1,sample 2....)
##returns F-value and p-value 
f_vals =[]
p_vals =[]
for row in range(len(res_df.index)):
	f, p = stats.f_oneway(res_df.iloc[row][0:4],res_df.iloc[row][4:8],res_df.iloc[row][8:12], res_df.iloc[row][12:16],res_df.iloc[row][16:])
	f_vals.append(f)
	p_vals.append(p)
	
res_df = res_df.assign(f_value = f_vals)
res_df = res_df.assign(p_value = p_vals)


#To correct for multiple hypothesis testing - use q-values
##QVALUE = min(p-value * total_number_of_peptides/rank_of_current_peptide, qvalue_for_peptide_one_row_below)
res_df.sort_values(by = ["p_value"], axis = 0, ascending = True, inplace = True)
q_vals = []
sign = []
##starting from the lowest ranked p-value:
prev_q = res_df.iloc[-1][-1]
q_vals.append(prev_q)
if prev_q < q_sign:
	sign.append("Yes")
else:
	sign.append("No")
##minimum for those with higher rank
for rank in range (len(res_df.index)-1,0,-1):
	q_val = min(res_df.iloc[rank-1][-1]*len(res_df.index)/rank , prev_q)
	q_vals.append(q_val)
	prev_q = q_val
	##Find significant results
	if q_val < q_sign:
		sign.append("Yes")
	else:
		sign.append("No")

#add lists to dataframe (lists need to reversed)		
res_df = res_df.assign(q_value = q_vals[::-1])
res_df = res_df.assign(significant = sign[::-1])

#output dataframe
res_df.to_csv(outpath_q)
sign_df = res_df.loc[res_df['significant'] == "Yes"]
sign_df.to_csv(outpath_sign)

