import csv

import numpy as np
import pandas as pd
import scipy.stats as sp
import scipy.optimize as spo
import matplotlib.pyplot as plt

# creating pandas dataframe
brca = pd.read_csv("brca.csv")
assay = pd.read_csv("assay.csv")


# merging dataframes on important values
merged_df = assay.merge(brca[["Unnamed: 0", "Genomic_Coordinate", "Chromosome"]], on="Unnamed: 0")

sorted_df = merged_df.sort_values(["Chromosome"])

# filtering dataframe 1 on chromosome value and genomic coordinate sorting
# ch_1_df = sorted_df[sorted_df["Chromosome"] == "1"]   ###########change this for specific CHROMOSOME!!!!#######
ch_1_df = sorted_df


ch_1_df_sorted_val = ch_1_df.sort_values(["Genomic_Coordinate"])
print(ch_1_df_sorted_val)


# calculating when or when not neighbors
distance = []
correlation_v = []



for i in range(len(ch_1_df_sorted_val["Genomic_Coordinate"]) - 1):
    first_coord = ch_1_df_sorted_val.iloc[i]
    second_coord = ch_1_df_sorted_val.iloc[i + 1]


    if  second_coord["Genomic_Coordinate"] - first_coord["Genomic_Coordinate"] <= 8000:
        dis = second_coord["Genomic_Coordinate"] - first_coord["Genomic_Coordinate"]
        distance.append(dis)
        first_calc_val = first_coord.iloc[1:50]
        # print(first_calc_val)
        second_calc_val = second_coord.iloc[1:50]

        #calc pearson
        cor_v, pval = sp.pearsonr(first_calc_val,second_calc_val)
        correlation_v.append(abs(cor_v))


# distance.sort()
# correlation_v.sort()

distances = np.array(distance)
correlation_v = np.array(correlation_v)
# print(distance)




distance = np.array(distance, dtype=float)

# slope, intercept, r, p, se = sp.linregress(distance,correlation_v)
slope, intercept, r_value, p_value, std_err = sp.linregress(distance, correlation_v)

# Calculate the regression line
regression_line = slope * distance + intercept

''' Background plot'''



# plt.scatter(distance, correlation_v, label='Data')
plt.plot(distance, regression_line, color='red', label='Linear Regression Line')
plt.ylim(0.1,0.6)
plt.xlabel('Distance')
plt.ylabel('Correlation Value')
plt.title("Methylation neighbours distance vs correlation value")
plt.legend()
plt.savefig("Regression_line.png")
plt.show()






