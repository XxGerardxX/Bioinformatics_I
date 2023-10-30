import pandas as pd

# imporing the nessecart csv
other_body = pd.read_csv("Other.csv")
SSI = pd.read_csv("island_data.csv")
SNPs = pd.read_csv("SNPs.csv")


def right_csv():
    merged_df_1 = pd.merge(other_body, SSI)
    merged_df_2 = pd.merge(merged_df_1, SNPs)

    # print(merged_df_2.columns)

    dataframe = merged_df_2[['Unnamed: 0','UCSC_RefGene_Group','Islands_Name', 'Relation_to_Island','Probe_SNPs', 'Probe_SNPs_10']]
    # print(dataframe)
    dataframe.to_csv("nine_features.csv", encoding='utf-8', index=False)


# importing newly created csv
nine_f_data = pd.read_csv("nine_features.csv")


#creating new dataframe
nine_features_df = pd.DataFrame(columns=['CsG','Island', 'Shore', 'Shelf', 'OpenSea','Probe','Probe 10','Promotor','Gene Body', 'Intergenic Region' ])

#first column
nine_features_df["CsG"] = nine_f_data['Unnamed: 0']

# Shore Shelf Island and OpenSea
unique_Island_items = ['S_Shore','N_Shore','S_Shelf','N_Shelf','Island','OpenSea']

nine_features_df['Shore'] = nine_f_data['Relation_to_Island'].apply(lambda x: 1 if x in ['S_Shore', 'N_Shore'] else 0)
nine_features_df['Shelf'] = nine_f_data['Relation_to_Island'].apply(lambda x: 1 if x in ['S_Shelf','N_Shelf'] else 0)
nine_features_df['Island'] = nine_f_data['Relation_to_Island'].apply(lambda x: 1 if x == 'Island' else 0)
nine_features_df['OpenSea'] = nine_f_data['Relation_to_Island'].apply(lambda x: 1 if x == 'OpenSea' else 0)

# SNPs
nine_features_df['Probe'] = nine_f_data['Probe_SNPs'].apply(lambda x: 1 if type(x) == str else 0)
nine_features_df['Probe 10'] = nine_f_data['Probe_SNPs_10'].apply(lambda x: 1 if type(x) == str else 0)

# promotor, gene body and intergenic region
nine_features_df['Promotor'] = nine_f_data['UCSC_RefGene_Group'].apply(lambda x: 1 if 'TSS1500' in str(x) or 'TSS200' in str(x) else 0)
nine_features_df['Gene Body'] = nine_f_data['UCSC_RefGene_Group'].apply(lambda x: 1 if "5'UTR" in str(x) or "3'UTR" in str(x) or "Body" in str(x) or "Exon" in str(x) else 0)

# for intergenic region
not_intergenic = ['TSS1500','TSS200',"5'UTR" ,"3'UTR","Body","Exon"]
nine_features_df['Intergenic Region'] = nine_f_data['UCSC_RefGene_Group'].apply(lambda x: 0 if any(item in str(x) for item in not_intergenic) else 1)




nine_features_df.to_csv("final_nine.csv", encoding='utf-8', index=False)
print(nine_features_df)











