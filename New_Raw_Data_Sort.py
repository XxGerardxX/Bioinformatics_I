import pandas as pd

pd.set_option('display.max_columns', None)

'''reading in all the CSV files'''
LIHC = pd.read_csv("brca.csv")
assay = pd.read_csv("assay.csv")
boolean_assay = pd.read_csv("Boolean_assay.csv")
nine_values = pd.read_csv("final_nine.csv")
rest_values = pd.read_csv("48_peaks_allignment.csv")


# TODO add this function as a subfunction to the merge_df main function
# TODO make the sort function for multiple chromosomes possible
# TODO make filter for multiple patients
def dataframe_filter(dataframe, chromosome_filter = "no", patient_filter="no", b_values_binary = "no" , creating_csv = "no", sorting = "no"):
    """Input: Dataframe
    Optional: chromosome filtering (default no) patient filtering (default = no)
    Output: dataframe filtered on chromosome or specific patient column"""

    dataframe = dataframe
    if chromosome_filter == "yes":
        try:
            all_chromosomes = dataframe["Chromosome"]
            all_chromosomes = set(all_chromosomes.tolist())
            which_chromosome = input(f'from the following list of chromosomes: {all_chromosomes} which chromosome do you want? ')
            dataframe = dataframe[dataframe["Chromosome"] == which_chromosome]
        except Exception as e:
            print(f"The following error occurred {e}")

    if patient_filter =="yes":
        dataframe_columns = dataframe.columns.tolist()
        which_patient_index = int(input(f'Which patient do you want to choose from {dataframe_columns}?\n Give a number '))
        name_column_index = int(input(f"which column contains the index names for ever row? (The CsG's) " ))
        try:
            chosen_column = dataframe.iloc[:,which_patient_index]
            name_column = dataframe.iloc[:,name_column_index]
            dataframe = pd.DataFrame({"CsG": name_column, "Patient_1": chosen_column})
        except Exception as e:
            print(f"The following error occurred {e}")

    if b_values_binary == "yes":
        dataframe['Patient_1'] = dataframe['Patient_1'].apply(lambda x: 0 if x < 0.5 else 1)

    if sorting == "yes":
        try:
            distance_column_name = "Genomic_Coordinate"
            dataframe = dataframe.sort_values(by=distance_column_name)
        except Exception as e:
            print(f"The following error occurred {e}")

    # creating a CSV output

    if creating_csv == "yes":
        name_csv = input("What is the name of this csv? ")

        try:
            dataframe.to_csv(name_csv, encoding='utf-8', index=False)
        except:
            print("The csv creation failed")
        else:
            print("The csv is created!")



    return dataframe




# TODO: Make merge.df cleaner by implementing new method for changing Unnamed:0 (eg put Unnamed:0 in a list and every
#  item in that list gets changed to x where x is another name)
def merge_df(two_df_list, values_added_to_merge_df_list, merging_column_index, creating_csv="no"):
    """Merges Dataframes and creates the necessary csv
    Input: two dataframes, the values of df two that need to be merged on to df one, the merging column index parameter
    Optional: sorting the data (default = no), creating csv (default = no)
    Output: dataframe two merged with dataframe one"""
    merging_column_index = int(merging_column_index)
    creating_csv = str.casefold(creating_csv)

    if all(isinstance(item, str) for item in values_added_to_merge_df_list) == False:
        print("not all values in the second list are strings")

    if "Unnamed: 0" in values_added_to_merge_df_list:
        unnamed_index = values_added_to_merge_df_list.index("Unnamed: 0")
        values_added_to_merge_df_list[unnamed_index] = "CsG"



    '''Important for changing the unnamed column to a named one'''
    if "Unnamed: 0" in two_df_list[0].columns.tolist():
        try:
            two_df_list[0].rename(columns={"Unnamed: 0": "CsG"}, inplace=True)
        except:
            print("name change didn't work for the first dataframe")
        else:
            print("Unnamed:0 has been changed to CsG in dataframe 1")


    if "Unnamed: 0" in two_df_list[1].columns.tolist():
        try:
            two_df_list[1].rename(columns={"Unnamed: 0": "CsG"}, inplace=True)
        except:
            print("name change didnt work for the second dataframe")
        else:
            print("Unnamed:0 has been changed to CsG in dataframe 2")

    column_names_first_df = two_df_list[0].columns  # called it twice because when the value from the first column changes the variable value changes as well

    # merging the dataframes
    merged_df = None

    try:
        merge_on_what_val = column_names_first_df[merging_column_index]
        merged_df = two_df_list[0].merge(two_df_list[1][values_added_to_merge_df_list], on=merge_on_what_val)
    except Exception as e:
        print(f"Dataframes can't be merged \n this is because: {e}")
    else:
        print("The merging worked!")

    # creating a CSV output

    if creating_csv == "yes":
        name_csv = input("What is the name of this csv? ")

        try:
            merged_df.to_csv(name_csv, encoding='utf-8', index=False)
        except:
            print("The csv creation failed")
        else:
            print("The csv is created!")

    return merged_df

def neighbouring_distances(dataframe_with_distances, distance_column_name, name_column = "CsG"):
    """Calculates the distances between the neighbouring chromosomes
    Input: dataframe with distances, the column name where the distances are in
    Output: 2 lists, up and downstream methylation"""
    dataframe_with_distances = dataframe_with_distances.sort_values(by = distance_column_name)
    CsG_Names = dataframe_with_distances.loc[:,name_column]
    upstream_distances = []
    downstream_distances = []




    for distance in range(len(dataframe_with_distances[distance_column_name]) - 1):
        # setting up distance calculations for up and downstream distances
        main_coordinate = dataframe_with_distances.iloc[distance]
        upstream_coordinate = dataframe_with_distances.iloc[distance - 1]
        downstream_coordinate = dataframe_with_distances.iloc[distance + 1]

        # calculating the upstream and downstream value
        upstream_value = abs(upstream_coordinate[distance_column_name] - main_coordinate[distance_column_name])
        downstream_value = abs(downstream_coordinate[distance_column_name] - main_coordinate[distance_column_name])

        upstream_distances.append(upstream_value)
        downstream_distances.append(downstream_value)

    upstream_distances.append(None)
    downstream_distances.append(None)

    calculated_distances_df = pd.DataFrame({"CsG":CsG_Names,"Upstream_Dis": upstream_distances, "Downstream_Dis": downstream_distances})

    return calculated_distances_df


def neighbouring_methyl_val(dataframe_with_methylation_vals, methylation_column_name, name_column="CsG"):
    CsG_Names = dataframe_with_methylation_vals[name_column]
    upstream_m_values = []
    downstream_m_values = []

    for m_val in range(len(dataframe_with_methylation_vals[methylation_column_name]) - 1):
        # setting up distance calculations for up and downstream distances
        main_m_val = dataframe_with_methylation_vals.iloc[m_val]
        upstream_m_val = dataframe_with_methylation_vals.iloc[m_val - 1]
        downstream_m_val = dataframe_with_methylation_vals.iloc[m_val + 1]

        upstream_m_values.append(upstream_m_val[methylation_column_name])
        downstream_m_values.append(downstream_m_val[methylation_column_name])

    upstream_m_values.append(None)
    downstream_m_values.append(None)

    calculated_M_Vals_df = pd.DataFrame({
        name_column: CsG_Names,
        "Upstream_M_Vals": upstream_m_values,
        "Downstream_M_Vals": downstream_m_values
    })

    return calculated_M_Vals_df



if __name__ == "__main__":

    #0 which patient do you want?
    ## filter the assay data (csv was created)
    filtered_assay_data = dataframe_filter(assay, chromosome_filter="no",patient_filter="yes",b_values_binary="yes") # assay with binary value and patient 1 is created

    # 1) which patient and which chromosome?
    ## merge LIHC and the filtered assay (csv was created to test)
    merge_assay_columns = filtered_assay_data.columns.tolist()
    merged_LIHC_filtered_assay = merge_df([LIHC,filtered_assay_data],merge_assay_columns,0)

    # 2) append the rest of the 9+48 values (csv was created to test) for both merge and join
    merged_9_columns = nine_values.columns.tolist()
    merged_9_to_df = merge_df([merged_LIHC_filtered_assay,nine_values],merged_9_columns,0)

    merged_9_48 = merged_9_to_df.join(rest_values)


    # 3) filter on chromosome 1 (the csv was created to test)
    final_df_chr_1 = dataframe_filter(merged_9_48, chromosome_filter="yes") #creates unnamed column because of first column 48 features for this drop the unnamed column
    # print(final_df_chr_1.columns)
    final_df_chr_1 = final_df_chr_1.drop("Unnamed: 0",axis=1)

    ## the data sorting on distance
    final_df_sorted_distance = neighbouring_distances(final_df_chr_1, distance_column_name="Genomic_Coordinate") # creating df with distances


    #4) calculate distances and methylation values
    final_distance_names = final_df_sorted_distance.columns.tolist()
    merged_distances = merge_df([final_df_chr_1,final_df_sorted_distance],final_distance_names, 0)
    merged_distances = dataframe_filter(merged_distances, sorting="yes")

    # calculating the m values

    adding_to_final_methylation_status = neighbouring_methyl_val(merged_distances,"Patient_1")
    adding_to_final_methylation_status.to_csv("aaa.csv")
    final_methylation_names = adding_to_final_methylation_status.columns.tolist()
    merged_methylation = merge_df([merged_distances,adding_to_final_methylation_status],final_methylation_names,0,creating_csv="yes")









