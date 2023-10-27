import pandas as pd

'''Reading in the csv.files'''

# creating pandas dataframe
brca = pd.read_csv("brca.csv")
boolean_val= pd.read_csv("Boolean_assay.csv")
assay_values = pd.read_csv("assay.csv")


'''Sorting and calculating the assay csv'''

merged_df = boolean_val.merge(brca[["Unnamed: 0", "Genomic_Coordinate", "Chromosome"]], on="Unnamed: 0")
sorted_df = merged_df.sort_values(["Chromosome"]) # sorts values based on what chromosome
sorted_df = sorted_df[sorted_df["Chromosome"] == "1"] # throws filter on df that only chromosome 1 is taken for
# calculations
distance_df = sorted_df.sort_values("Genomic_Coordinate")


df_column_names = distance_df.columns
df_column_names = df_column_names[1:51]

# lists for the columns
b_values_patient = []
upstream_distance_l = []
downstream_distance_l = []
upstream_m_stat_l = []
downstream_m_stat_l = []



'''Calculating up and downstream distances'''
for coordinate in range(len(distance_df["Genomic_Coordinate"])-1):

    # setting up distance calculations for up and downstream distances
    main_coordinate = distance_df.iloc[coordinate]
    upstream_coordinate = distance_df.iloc[coordinate - 1]
    downstream_coordinate = distance_df.iloc[coordinate + 1]


    #calculating the upstream and downstream value
    upstream_value = abs(upstream_coordinate["Genomic_Coordinate"] - main_coordinate["Genomic_Coordinate"])
    downstream_value = abs(downstream_coordinate["Genomic_Coordinate"] - main_coordinate["Genomic_Coordinate"])

    upstream_distance_l.append(upstream_value)
    downstream_distance_l.append(downstream_value)


    #calculating upstream and downstream m status
    #Taking patient 1 to 50


for methylation_val in range(len(distance_df[df_column_names[0]])-1):   # change df columnnames to change the patient
    main_m_index = distance_df.iloc[methylation_val]
    upstream_m_index = distance_df.iloc[methylation_val - 1]
    downstream_m_index = distance_df.iloc[methylation_val + 1]

    upstream_m_value = upstream_m_index[df_column_names[0]]
    downstream_m_value = downstream_m_index[df_column_names[0]]

    upstream_m_stat_l.append(upstream_m_value)
    downstream_m_stat_l.append(downstream_m_value)







#appending NaN values
upstream_distance_l.append(None)
downstream_distance_l.append(None)

upstream_m_stat_l.append(None)
downstream_m_stat_l.append(None)

upstream_distance_l[0] = None

#appending the columns to the dataframe
distance_df["Upstream_distance"] = upstream_distance_l
distance_df["Downstream_distance"] = downstream_distance_l
distance_df["Upstream_methylation"] = upstream_m_stat_l
distance_df["Downstream_methylation"] = downstream_m_stat_l




data = {"B_Val": b_values_patient,"Upstream_distance": upstream_distance_l, "Downstream_distance": downstream_distance_l,
        "Upstream_methylation": upstream_m_stat_l,"Downstream_methylation": downstream_m_stat_l  }

final_df = pd.DataFrame(data)
#creating the dataframe csv

final_df.to_csv("DF_With_Distances_2.csv", encoding='utf-8', index=False)
