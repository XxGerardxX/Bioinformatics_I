import csv
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)

# TODO: add if name == Main

'''reading in all the CSV files'''

brca = pd.read_csv("brca.csv")
assay = pd.read_csv("assay.csv")
boolean_val = pd.read_csv("Boolean_assay.csv")
import_nine_values = pd.read_csv("nine_features.csv")
import_48_values = pd.read_csv("48_peaks_allignment.csv")


def merge_df(two_df_list,values_added_to_merge_df_list,merge_on_what_val):
    merge_on_what_val = str(merge_on_what_val)
    if all(isinstance(item,str) for item in values_added_to_merge_df_list) == False:
        print("not all values in the second list are strings")

    try:
        merged_df = two_df_list[0].merge(two_df_list[1][values_added_to_merge_df_list], on=merge_on_what_val)
    except:
        print("Lists can't be merged")


    ask_filter = input("Do you want to fiter the data? ")
    ask_filter = str.casefold(ask_filter)

    if ask_filter == "yes":
        filtering_chr_value = str(input("On which chromosome needs the data to be filtered? "))
        sorting_ask = input("Does the chromosome data need to be sorted on distance? ")

        try:
            merged_df = merged_df[merged_df["Chromosome"] == filtering_chr_value]
            if sorting_ask == "yes":
                merged_df = merged_df.sort_values("Genomic_Coordinate")
        except:
            print("the filtering failed")



    return merged_df




'''merging boolean and non boolean dataframes'''
need_to_be_merged_bool = [boolean_val, brca]
merged_on_values_bool = ["Unnamed: 0", "Genomic_Coordinate", "Chromosome"]
merged_on_bool = "Unnamed: 0"
merged_df = merge_df(need_to_be_merged_bool, merged_on_values_bool,merged_on_bool)

merged_df.to_csv("test_boolean_2.csv", encoding='utf-8', index=False)








