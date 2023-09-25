import csv
import numpy as np
import pandas as pd
import scipy.stats as sp
import scipy.optimize as spo
import matplotlib.pyplot as plt


assay = pd.read_csv("assay.csv")

# get column names
column_names = assay.columns
# print(column_names)





# creating new dataframe without first column
column_names = column_names[1:]
altered_df = assay[column_names]



# Define a function to apply to each cell
def threshold_function(cell_value):
    if cell_value < 0.5:
        return 0
    else:
        return 1

# Apply the function to each cell in the DataFrame
new_df = altered_df.applymap(threshold_function)
new_df.to_csv("Boolean_assay.csv", encoding='utf-8', index=False)

# Display the new DataFrame
print(new_df)




