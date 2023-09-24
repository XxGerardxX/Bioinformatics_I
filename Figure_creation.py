import csv
import pandas as pd
import scipy.stats as sp
import matplotlib.pyplot as plt

# Methylation correlations between CpG sites were assessed by the absolute
# value of Pearson’s correlation coefficient and MED

# x to 200 bp are the correlations calculated
'''
Make a figure similar to Fig. 1A in Zhang et al. to see how the correlation between the methylation profiles of
neighbouring CpGs varies with the distance between them. No need to make a distinction between ‘CGI’,
‘CGIshore&shelf’ and ‘nonCGI’ but the background estimate indicated by the solid line is useful.
'''

# creating pandas dataframe
brca = pd.read_csv("brca.csv")
assay = pd.read_csv("assay.csv")

# # data for the same chromosome is only relevant
# chromosome_column = brca['Chromosome']
#
# # unique chromosomes used to calculate pearson correlation
# chromosome_list = sorted(list(set(chromosome_column)))
#
#
# # values connected to the specific chromosome
# matrix_list = []
#
# for value in chromosome_list:
#     chromosome_filter_matrix = brca[chromosome_column == str(value)]
#     chromosome_filter_matrix = chromosome_filter_matrix.sort_values(by=["Genomic_Coordinate"])
#     matrix_list.append(chromosome_filter_matrix)
#
#
#
# # names list for searching in assay matrix
# ordered_names_list = list(matrix_list[0]["Unnamed: 0"])
#
#
# ### not all values of brca (ordered names list) are in assay>??!!
# filtered_df = assay[assay["Unnamed: 0"].isin(ordered_names_list)]
# filtered_unnamed_list = list(filtered_df["Unnamed: 0"])
#
# filtered_matrix_list = brca[brca["Unnamed: 0"].isin(filtered_unnamed_list)]["Genomic_Coordinate"]


# filtered_df.insert(1,"Genomic_Coordinate",filtered_matrix_list)
#########################################################################################################################

# merging dataframes on important values
merged_df = assay.merge(brca[["Unnamed: 0", "Genomic_Coordinate", "Chromosome"]], on="Unnamed: 0")

sorted_df = merged_df.sort_values(["Chromosome"])

# filtering dataframe 1 on chromosome value and genomic coordinate sorting
# ch_1_df = sorted_df[sorted_df["Chromosome"] == "1"]
ch_1_df = sorted_df



ch_1_df_sorted_val = ch_1_df.sort_values(["Genomic_Coordinate"])
print(ch_1_df_sorted_val)


# calculating when or when not neighbors
distance = []
correlation_v = []



for i in range(len(ch_1_df_sorted_val["Genomic_Coordinate"]) - 1):
    first_coord = ch_1_df_sorted_val.iloc[i]
    second_coord = ch_1_df_sorted_val.iloc[i + 1]


    if 200 <= second_coord["Genomic_Coordinate"] - first_coord["Genomic_Coordinate"] <= 6000:
        dis = second_coord["Genomic_Coordinate"] - first_coord["Genomic_Coordinate"]
        distance.append(dis)
        first_calc_val = first_coord.iloc[1:50]
        # print(first_calc_val)
        second_calc_val = second_coord.iloc[1:50]

        #calc pearson
        cor_v, pval = sp.pearsonr(first_calc_val,second_calc_val)
        correlation_v.append(abs(cor_v))


distance.sort()
# correlation_v.sort(reverse=True)
correlation_v.sort()
print(distance)
print(correlation_v)

plt.plot(distance,correlation_v)
plt.show()




