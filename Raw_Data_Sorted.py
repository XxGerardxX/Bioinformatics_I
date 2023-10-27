import csv
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)

import scipy.stats as sp
import scipy.optimize as spo
import matplotlib.pyplot as plt


'''MAKE FUNCTIONS FOR CODE SO ITS REUSABLE INSTEAD OF HOT GARBAGE ::))))'''
# importing csv files
brca = pd.read_csv("brca.csv")
assay = pd.read_csv("assay.csv")
boolean_val = pd.read_csv("Boolean_assay.csv")
import_nine_values = pd.read_csv("nine_features.csv")
import_48_values = pd.read_csv("48_peaks_allignment.csv")




'''creating specifically for a chromosome the merged dataset for the features (52)'''




'''merging boolean and non boolean dataframes and sorting them on chromosome 1'''
#merging boolean dataframes



merged_boolean_df = boolean_val.merge(brca[["Unnamed: 0", "Genomic_Coordinate", "Chromosome"]], on="Unnamed: 0")
sorted_boolean_df = merged_boolean_df.sort_values(["Chromosome"]) # sorts values based on what chromosome
sorted_boolean_df = sorted_boolean_df[sorted_boolean_df["Chromosome"] == "1"] # throws filter on df that only chromosome 1 is taken for
sorted_on_distance_boolean = sorted_boolean_df.sort_values("Genomic_Coordinate")
sorted_on_distance_boolean.to_csv("test_boolean.csv", encoding='utf-8', index=False)


#merging non_boolean_dataframes
merged_df = assay.merge(brca[["Unnamed: 0", "Genomic_Coordinate", "Chromosome"]], on="Unnamed: 0")
sorted_df = merged_df.sort_values(["Chromosome"])
sorted_df = sorted_df[sorted_df["Chromosome"] == "1"]
sorted_on_distance = sorted_df.sort_values("Genomic_Coordinate")

sorted_on_distance.to_csv("test.csv", encoding='utf-8', index=False)

##########################################################################
'''Lists for new columns for the new dataframe'''
# Lists for columns
b_values_patient = []
upstream_distances = []
downstream_distances = []
upstream_m_stats = []
downstream_m_stats = []


# setting up a patient between 1 and 50
patient_number = 0

'''Getting for the specific patient the beta values of the non-boolean dataframe'''
df_column_names = sorted_on_distance.columns
df_column_names = df_column_names[1:51]

one_patient_df = sorted_on_distance[df_column_names[patient_number]]
for i in one_patient_df:
    b_values_patient.append(i)


'''Getting the distances for a specific patient between the chromosomes'''

for coordinate in range(len(sorted_on_distance_boolean["Genomic_Coordinate"])-1):

    # setting up distance calculations for up and downstream distances
    main_coordinate = sorted_on_distance_boolean.iloc[coordinate]
    upstream_coordinate = sorted_on_distance_boolean.iloc[coordinate - 1]
    downstream_coordinate = sorted_on_distance_boolean.iloc[coordinate + 1]


    #calculating the upstream and downstream value
    upstream_value = abs(upstream_coordinate["Genomic_Coordinate"] - main_coordinate["Genomic_Coordinate"])
    downstream_value = abs(downstream_coordinate["Genomic_Coordinate"] - main_coordinate["Genomic_Coordinate"])

    upstream_distances.append(upstream_value)
    downstream_distances.append(downstream_value)



'''Calculating the binary methylation values for a specific patient'''

for methylation_val in range(len(sorted_on_distance_boolean[df_column_names[patient_number]])-1):   # change df columnnames to change the patient
    main_m_index = sorted_on_distance_boolean.iloc[methylation_val]
    upstream_m_index = sorted_on_distance_boolean.iloc[methylation_val - 1]
    downstream_m_index = sorted_on_distance_boolean.iloc[methylation_val + 1]

    upstream_m_value = upstream_m_index[df_column_names[patient_number]]
    downstream_m_value = downstream_m_index[df_column_names[patient_number]]

    upstream_m_stats.append(upstream_m_value)
    downstream_m_stats.append(downstream_m_value)


# print(len(b_values_patient))
# print(len(upstream_m_stats))
# print(len(downstream_m_stats))
# print(len(upstream_distances))

#appending for 1 extra value

upstream_distances.append(None)
downstream_distances.append(None)
upstream_m_stats.append(None)
downstream_m_stats.append(None)

'''making the new dataframe'''

data = {"B_Val": b_values_patient,"Upstream_distance": upstream_distances, "Downstream_distance": downstream_distances,
        "Upstream_methylation": upstream_m_stats,"Downstream_methylation": downstream_m_stats  }

final_df = pd.DataFrame(data)
final_df.to_csv("DF_With_Distances.csv", encoding='utf-8', index=False)

