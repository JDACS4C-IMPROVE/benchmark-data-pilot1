from pathlib import Path
import os
import pandas as pd
import numpy as np



fdir = Path(__file__).resolve().parent
# pd.set_option('display.max_columns', None)

# drug_data_dir = '../../Drug_Data_AID/'
drug_data_dir = fdir/'../../drug_data_aid/'

# benchmark_data_dir = '../CSA_Data/'
benchmark_data_dir = fdir/'../csa_data/raw_data'
os.makedirs(benchmark_data_dir, exist_ok=True)

x_data_dir = benchmark_data_dir/'x_data'
os.makedirs(x_data_dir, exist_ok=True)



# Load drug fingerprint data
# df = pd.read_parquet(drug_data_dir + 'ecfp4_nbits512.parquet', engine='pyarrow')
df = pd.read_parquet(drug_data_dir/'ecfp4_nbits512.parquet', engine='pyarrow')
print("ECFP4", df.shape)
df.index = df.DrugID
df = df.iloc[:, 1:]
print(df.drop_duplicates().shape[0]) # 1565

# Load drug descriptor data
# dd = pd.read_parquet(drug_data_dir + 'mordred.parquet', engine='pyarrow')
dd = pd.read_parquet(drug_data_dir/'mordred.parquet', engine='pyarrow')
print("Mordred", dd.shape)
dd.index = dd.DrugID
dd = dd.iloc[:, 1:]
print(dd.drop_duplicates().shape[0]) # 1566

# Load drug information
# drug_info = pd.read_csv(drug_data_dir + 'drug_meta.tsv', sep='\t', engine='c', na_values=['na', '-', ''],
#                         header=0, index_col=None, low_memory=False)
drug_info = pd.read_csv(drug_data_dir/'drug_meta.tsv', sep='\t', engine='c', na_values=['na', '-', ''],
                        header=0, index_col=None, low_memory=False)
print("Drug info", drug_info.shape)
print(len(np.unique(drug_info.SMILES))) # 1781
print(len(np.unique(drug_info.canSMILES))) # 1664


# # This code checks whether the same canSMILES always have the same fingerprint and descriptor profiles. Result is Yes.
# drug_info.index = drug_info.DrugID
# index = 0
# for i in np.unique(drug_info.canSMILES):
#     index = index + 1
#     print(index)
#     idi = np.where(drug_info.canSMILES == i)[0]
#     iddf = np.where(np.isin(df.index, drug_info.index[idi]))[0]
#     if len(iddf) > 1:
#         if df.iloc[iddf, :].drop_duplicates().shape[0] > 1:
#             print('More than one fingerprint profiles associated with' + str(i))
#     iddd = np.where(np.isin(dd.index, drug_info.index[idi]))[0]
#     if len(iddd) > 1:
#         if dd.iloc[iddd, :].drop_duplicates().shape[0] > 1:
#             print('More than one descriptor profiles associated with' + str(i))


# # This code identifies a descriptor profile that associates with more than one canSMILES. Result shows that
# # such a descriptor profile correspond to one drug, although the canSMILES can be different.
# index = 0
# dd_str = dd.iloc[:, 0].map(str)
# for i in range(1, dd.shape[1]):
#     dd_str = dd_str + '|' + dd.iloc[:, i].map(str)
# for i in np.unique(dd_str):
#     index = index + 1
#     print(index)
#     idi = np.where(dd_str == i)[0]
#     idinfo = np.where(np.isin(drug_info.DrugID, dd.index[idi]))[0]
#     if len(np.unique(drug_info.iloc[idinfo, :].canSMILES)) > 1:
#         print(drug_info.iloc[idinfo, :])


# # This code identifies the fingerprint profile that associates with more than on descriptor profile.
# # Result shows that the two descriptor profiles correspond to the same drug.
# index = 0
# df_str = df.iloc[:, 0].map(str)
# for i in range(1, df.shape[1]):
#     df_str = df_str + '|' + df.iloc[:, i].map(str)
# for i in np.unique(df_str):
#     index = index + 1
#     print(index)
#     idi = np.where(df_str == i)[0]
#     iddd = np.where(np.isin(dd.index, df.index[idi]))[0]
#     if len(iddd) > 1:
#         if dd.iloc[iddd, :].drop_duplicates().shape[0] > 1:
#             print(dd.iloc[iddd, :].index)

# So we decide to use fingerprint profiles to define unique drugs.


# # This code checks whether the same SMILES always have the same fingerprint and descriptor profiles.
# uni_smiles = np.unique(drug_info.SMILES)
# # for u in uni_smiles:
# #     uid = np.where(drug_info.SMILES == u)[0]
# #     if len(uid) > 1:
# #         dd_id = np.where(np.isin(dd.index, drug_info.iloc[uid, 0]))[0]
# #         print(dd.iloc[dd_id, :].drop_duplicates().shape)
# #         df_id = np.where(np.isin(df.index, drug_info.iloc[uid, 0]))[0]
# #         print(df.iloc[df_id, :].drop_duplicates().shape)



# Based on fingerprints, generate unique drug IDs
df_str = df.iloc[:, 0].map(str)
for i in range(1, df.shape[1]):
    df_str = df_str + '|' + df.iloc[:, i].map(str)
uni_df_str = np.unique(df_str)
drug_uni_id = pd.Series(['' for i in range(df.shape[0])])
for i in range(len(uni_df_str)):
    uid = np.where(df_str == uni_df_str[i])[0]
    drug_uni_id.iloc[uid] = ['Drug_' + str(i) for j in range(len(uid))]
drug_uni_id.index = df.index
drug_info['improve_chem_id'] = drug_uni_id.loc[drug_info.DrugID].values
# drug_info.to_csv(benchmark_data_dir + 'drug_info.txt', header=True, index=None, sep='\t', line_terminator='\r\n')
drug_info.to_csv(x_data_dir/'drug_info.tsv', header=True, index=None, sep='\t', line_terminator='\r\n')
uni_improve_chem_id, uni_index = np.unique(drug_info['improve_chem_id'], return_index=True)
drug_info = drug_info.iloc[uni_index, ]

df = df.loc[drug_info.DrugID, :]
print(np.sum(df.index != drug_info.DrugID))
df.index = drug_info.improve_chem_id

dd = dd.loc[drug_info.DrugID, :]
print(np.sum(dd.index != drug_info.DrugID))
dd.index = drug_info.improve_chem_id

# drug_info.loc[:, ['improve_chem_id', 'canSMILES']].to_csv(benchmark_data_dir + 'drug_SMILES.txt', header=True, index=False,
#                                                        sep='\t', line_terminator='\r\n')
# df.to_csv(benchmark_data_dir + 'drug_fingerprint.txt', header=True, index=True, sep='\t', line_terminator='\r\n')
# dd.to_csv(benchmark_data_dir + 'drug_descriptor.txt', header=True, index=True, sep='\t', line_terminator='\r\n')

drug_info.loc[:, ['improve_chem_id', 'canSMILES']].to_csv(x_data_dir/'drug_SMILES.tsv', header=True, index=False,
                                                       sep='\t', line_terminator='\r\n')
df.to_csv(x_data_dir/'drug_ecfp4_nbits512.tsv', header=True, index=True, sep='\t', line_terminator='\r\n')
dd.to_csv(x_data_dir/'drug_mordred.tsv', header=True, index=True, sep='\t', line_terminator='\r\n')


print("\nFinished generating drug data.")
