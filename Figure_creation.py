import csv

import numpy as np
import pandas as pd
import scipy.stats as sp
import scipy.optimize as spo
import matplotlib.pyplot as plt
import random
import statsmodels.api as sm


# creating pandas dataframe
brca = pd.read_csv("brca.csv")
assay = pd.read_csv("assay.csv")

# merging dataframes on important values
merged_df = assay.merge(brca[["Unnamed: 0", "Genomic_Coordinate", "Chromosome"]], on="Unnamed: 0")
sorted_df = merged_df.sort_values(["Chromosome"])

# filtering dataframe 1 on chromosome value and genomic coordinate sorting
ch_1_df = sorted_df[sorted_df["Chromosome"] == "1"]  ###########change this for specific CHROMOSOME!!!!#######
# ch_1_df = sorted_df

ch_1_df_sorted_val = ch_1_df.sort_values(["Genomic_Coordinate"])


########################################################################################################################
'''Calculating correlation values'''
distance = []
correlation_v = []

for i in range(len(ch_1_df_sorted_val["Genomic_Coordinate"]) - 1):
    first_coord = ch_1_df_sorted_val.iloc[i]
    second_coord = ch_1_df_sorted_val.iloc[i + 1]

    if second_coord["Genomic_Coordinate"] - first_coord["Genomic_Coordinate"] <= 8000:
        dis = second_coord["Genomic_Coordinate"] - first_coord["Genomic_Coordinate"]
        distance.append(dis)
        first_calc_val = first_coord.iloc[1:50]
        # print(first_calc_val)
        second_calc_val = second_coord.iloc[1:50]

        # calc pearson
        cor_v, pval = sp.pearsonr(first_calc_val, second_calc_val)
        correlation_v.append(abs(cor_v))


distances = np.array(distance)
correlation_v = np.array(correlation_v)

distance = np.array(distance, dtype=float)

# slope, intercept, r, p, se = sp.linregress(distance,correlation_v)
slope, intercept, r_value, p_value, std_err = sp.linregress(distance, correlation_v)

# Calculate the regression line
regression_line = slope * distance + intercept
####################################################################################################################
''' Background plot'''


# 50k samples random
random_dataframe_1 = sorted_df.sample(n=50000)
random_dataframe_2 = sorted_df.sample(n=50000)
values = [i for i in range(50000)]

# get column names
random_DF_column_names = random_dataframe_1.columns
random_DF_column_names = random_DF_column_names[1:51]

# sort the two dataframes on distances
random_dataframe_1 = random_dataframe_1.sort_values(["Genomic_Coordinate"])
random_dataframe_2 = random_dataframe_2.sort_values(["Genomic_Coordinate"])

cor_val_list = []


# calculating p values

for i in range(50000):
    first_df = random_dataframe_1.iloc[i]
    second_df = random_dataframe_2.iloc[i]

    first_df_calc_val = first_df.iloc[1:50]
    second_df_calc_val = second_df.iloc[1:50]
    cor_v_2, pval_2 = sp.pearsonr(first_df_calc_val, second_df_calc_val)
    cor_val_list.append(cor_v_2)


cor_val = np.mean(cor_val_list)
print(cor_val)

# location of coordinates lists
LoC_DF_1 = random_dataframe_1[["Genomic_Coordinate"]]
LoC_DF_2 = random_dataframe_2["Genomic_Coordinate"]

LoC_DF_1 = list(LoC_DF_1)
LoC_DF_2 = list(LoC_DF_2)




'''LOWESS'''
# Calculate the LOWESS regression line
lowess = sm.nonparametric.lowess(correlation_v, distance, frac=0.3)

# Extract the LOWESS smoothed values
lowess_x = lowess[:, 0]
lowess_y = lowess[:, 1]

# Sort the LOWESS values by distance
lowess_x, lowess_y = zip(*sorted(zip(lowess_x, lowess_y)))

# Plot the LOWESS regression line
plt.plot(lowess_x, lowess_y, color='green', label='LOWESS Regression Line')





# plotting the data
plt.plot(distance, regression_line, color='red', label='Linear Regression Line')
plt.axhline(y=cor_val, color='blue', linestyle='--', label=f'Correlation Value: {cor_val}')
plt.ylim(0.1, 0.6)
plt.xlabel('Distance')
plt.ylabel('Correlation Value')
plt.title("Methylation neighbours distance vs correlation value")
plt.legend()
plt.savefig("Regression_line.png")
plt.show()
