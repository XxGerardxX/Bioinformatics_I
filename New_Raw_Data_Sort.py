import pandas as pd

pd.set_option('display.max_columns', None)

'''reading in all the CSV files'''
brca = pd.read_csv("brca.csv")
assay = pd.read_csv("assay.csv")
boolean_assay = pd.read_csv("Boolean_assay.csv")
nine_values = pd.read_csv("final_nine.csv")
rest_values = pd.read_csv("48_peaks_allignment.csv")

# TODO: Make merge.df cleaner by implementing new method for changing Unnamed:0 (eg put Unnamed:0 in a list and every
#  item in that list gets changed to x where x is another name)

def merge_df(two_df_list, values_added_to_merge_df_list, merging_column_index, filtering="no", sorting="no", creating_csv="no"):
    """Merges Dataframes and creates the necessary csv
    Input: two dataframes, the values of df two that need to be merged on to df one, the merging column index parameter
    Optional: sorting the data (default = no), creating csv (default = no)
    Output: dataframe two merged with dataframe one"""
    merging_column_index = int(merging_column_index)
    sorting = str.casefold(sorting)
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

    # filtering the dataframe on a specific chromosome and checking if sorting is necessary and if dataframe value isn't none
    ask_filter = str.casefold(filtering)

    if ask_filter == "yes":
        filtering_chr_value = str(input("On which chromosome needs the data to be filtered? "))

        try:
            merged_df = merged_df[merged_df["Chromosome"] == filtering_chr_value]
            if sorting == "yes":
                merged_df = merged_df.sort_values("Genomic_Coordinate")
        except:
            print("the filtering failed (check hardcoded names)")

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




def neighbouring_distances(dataframe_with_distances, distance_column_name):
    """Calculates the distances between the neighbouring chromosomes
    Input: dataframe with distances, the column name where the distances are in
    Output: 2 lists, up and downstream methylation"""
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
    return upstream_distances, downstream_distances


''' The function neighbouring methylation value needs cleanup'''

def neighbouring_methylation_value(dataframe_with_methylation_values, non_bool_dataframe_methyl_values):
    """Give a dataframe with beta values for multiple patients.
     The first columns are the cg names and thus are not taken into account
     Input: Dataframe with methylation values from different patients,
     Output: Upstream and downstream methylation values in 2 lists and the b values"""

    b_values_patient = []
    upstream_methylation_values = []
    downstream_methylation_values = []

    patient_number = 0
    try:
        patient_number = int(input("which patient do you want? "))
    except:
        print("This is not a valid input!")

    '''Getting for the specific patient the beta values of the non-boolean dataframe'''
    df_column_names = non_bool_dataframe_methyl_values.columns
    df_column_names = df_column_names[0:51]


    one_patient_df = non_bool_dataframe_methyl_values[df_column_names[patient_number]]
    for i in one_patient_df:
        b_values_patient.append(i)

    '''Calculating the binary methylation values for a specific patient'''

    for methylation_val in range(len(dataframe_with_methylation_values[df_column_names[patient_number]]) - 1):
        main_m_index = dataframe_with_methylation_values.iloc[methylation_val]
        upstream_m_index = dataframe_with_methylation_values.iloc[methylation_val - 1]
        downstream_m_index = dataframe_with_methylation_values.iloc[methylation_val + 1]

        upstream_m_value = upstream_m_index[df_column_names[patient_number]]
        downstream_m_value = downstream_m_index[df_column_names[patient_number]]

        upstream_methylation_values.append(upstream_m_value)
        downstream_methylation_values.append(downstream_m_value)

    upstream_methylation_values.append(None)
    downstream_methylation_values.append(None)




    return upstream_methylation_values, downstream_methylation_values, b_values_patient


# TODO add this function as a subfunction to the merge_df main function
# TODO make the sort function for multiple chromosomes possible
# TODO make a patient filter


def dataframe_filter(dataframe, chromosome_filter = "no", patient_filter="no" ):
    dataframe = dataframe
    if chromosome_filter == "yes":
        all_chromosomes = dataframe["Chromosome"]
        all_chromosomes = set(all_chromosomes.tolist())
        which_chromosome = input(f'from the following list of chromosomes: {all_chromosomes} which chromosome do you want?')
        dataframe = dataframe[dataframe["Chromosome"] == which_chromosome]

    return dataframe





if __name__ == "__main__":

    '''merging non boolean dataframes'''
    need_to_be_merged = [assay, brca]
    merged_on_values = ["CsG", "Genomic_Coordinate", "Chromosome"]
    merged_on = 0
    merged_df = merge_df(need_to_be_merged, merged_on_values, merged_on, filtering="no", sorting="yes", creating_csv="no")



    '''merging boolean dataframes'''
    need_to_be_merged_bool = [boolean_assay, brca]
    merged_on_values_bool = ["CsG", "Genomic_Coordinate", "Chromosome"]
    merged_on_bool = 0
    merged_bool_df = merge_df(need_to_be_merged_bool, merged_on_values_bool, merged_on_bool, filtering="no",
                              sorting="yes", creating_csv="no")


    '''creating the big dataframe the 9 features, the patient data and the distances in it'''
    # ntbm_non_bool_nine = [merged_bool_df, nine_values]
    # merging_on_values_nine = nine_values.columns.tolist()


    # merged_nine_and_bool = merge_df(ntbm_non_bool_nine, merging_on_values_nine, merged_on, filtering="no", sorting="yes", creating_csv = "yes")

    nine_with_patient = pd.read_csv("Coordinates_Patient_Nine_Features.csv")
    '''Creating the up and downstream methylation data'''

    # using merged dataframes
    # upstream_distances, downstream_distances = neighbouring_distances(dataframe_with_distances=merged_bool_df,
    #                                                                   distance_column_name="Genomic_Coordinate")
    # upstream_methylation, downstream_methylation, beta_values_patient = neighbouring_methylation_value(merged_bool_df,
    #                                                                                                    merged_df)
    #
    # # creating the dataframe
    # data = {"B_Val": beta_values_patient, "Upstream_distance": upstream_distances,
    #         "Downstream_distance": downstream_distances,
    #         "Upstream_methylation": upstream_methylation, "Downstream_methylation": downstream_methylation}
    #
    # # creating the csv out of the dataframe. CSV is meant for the first basic random forest model
    # final_df = pd.DataFrame(data)
    # final_df.to_csv("DF_With_Distances_2.csv", encoding='utf-8', index=False)

    '''Adding the 48 columns to the rest of the data'''
    # columns_nine = nine_with_patient.loc[:,"CsG"]
    # rest_values.iloc[:,0] = columns_nine
    # named_first_column = pd.DataFrame(rest_values)
    # named_first_column.to_csv("Changed_with_48.csv", encoding='utf-8', index=False)
    #
    # changed_48 = pd.read_csv("Changed_with_48.csv")



    #merging 48 to the rest of the data
    # final_merge_datasets = [nine_with_patient,changed_48]
    # values_added_final = changed_48.columns.tolist()
    # merged_on_final = 0

    # final_merge = merge_df(final_merge_datasets, values_added_final , merged_on_final, filtering="no", sorting="yes", creating_csv = "yes")

    ''' calculating methylation values'''
    # neighbouring_methylation_value(boolean_assay, assay)
    dataframe_filter(brca, chromosome_filter="yes")



