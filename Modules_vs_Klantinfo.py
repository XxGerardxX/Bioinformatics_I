'''Importing nessecary modules'''
import pandas as pd

# doel -> als ze module hebben moeten ze ervoor betalen
try:
    modules = pd.read_excel("modulen_humanwave_b.v._2023-11-01_2023-11-30.xlsx", header=3)
    prijzen = pd.read_excel("prijzen_humanwave_b.v._2023-11-01_2023-11-30.xlsx", header=3)
except AssertionError as e:
    print(e)
else:
    modules = pd.DataFrame(modules)
    prijzen = pd.DataFrame(prijzen)

modules_columns = modules.columns
prijzen_columns = prijzen.columns
# extracting name and comparing them

'''WORKS! only nessecary when file column headers change of names, else unnecessary double for loop'''
# test_organisation = prijzen.loc[0, "Organisatie"]
# test_finder = modules[modules["Organisatie"] == test_organisation]

# checksum if column has data or not
# if pd.isna(modules.loc[1, "Salarisverwerking"]):
#     print("yes")


##############################


# creating dictionary for items # key = what needs to be paid for; value = active module
modulen_prijzen_namen = {
    "HW1100": ["Bedrijfsmiddelen", "Verlof", "Opleiding", "Declaraties", "Werknemer portaal", "Documenten e-mailen",
               "Nieuwsberichten", "Zoeken"], "HW1110": "Salaris", "HW1300": "Kostenverdeling", "HW1400": "Bitcare API",
    "HW1410": "Easy secure", "HW1500": "Ondertekenen documenten", "HW3100": "Salarisverwerking"}

# new items to be added
to_add = {"HW8000": "Standaard drempel", "HW8010": "Drempel full service"}

# merging the two dataframes on corresponding same columns
merged_df = modules.merge(prijzen, on=['Organisatie', 'Account', 'Van', 'Tot en met', "Status"])
merged_df.to_csv('test.csv')

count = 0
# iterating over keys (items the company must pay for)
for key, value in modulen_prijzen_namen.items():
    '''Loops over each key-value pair. For each pair it checks the values for the corresponding rows in the dataframe'''
    # HW1100 check
    for index, row in merged_df.iterrows():
        # gets the n-th row with the corresponding values of the necessary columns
        checked_value = merged_df.loc[index, value]
        checked_key = merged_df.loc[index, key]
        # x is type str nan is type float

        #checks if value is longer then 1. this is because the first value set is a lost (works!)
        try:
            if len(checked_value) > 1:
                print("longer then 1!")
                count += 1
        except:
            continue





print(count)