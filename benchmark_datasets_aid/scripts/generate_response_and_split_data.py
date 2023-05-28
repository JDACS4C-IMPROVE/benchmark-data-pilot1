import pandas as pd
import numpy as np
from collections import Counter

from Functions import generate_cross_validation_partition



pd.set_option('display.max_columns', None)

omics_data_dir = '../../Data_Curation_final_AID/Curated_CCLE_Multiomics_files/'

response_data_dir = '../../Response_Data_AID/'

drug_data_dir = '../../Drug_Data_AID/'

auxiliary_data_dir = '../Auxiliary_Data/'

benchmark_data_dir = '../CSA_Data/'



# Load response data
res = pd.read_csv(response_data_dir + '/combined_single_response_rescaled_agg', sep='\t', engine='c',
                  na_values=['na', '-', ''], header=0, index_col=None, low_memory=False)
res.columns = [str(i).lower() for i in res.columns]
id = np.where(np.invert(pd.isna(res.auc)))[0]
res = res.iloc[id, :]

# Counter(res.source)
# Counter({'CCLE': 11670,
#          'CTRP': 395263,
#          'GDSC': 225480,
#          'NCI60': 3780148,
#          'SCL': 28275,
#          'SCLC': 36789,
#          'gCSI': 6455,
#          'GDSC1': 322450,
#          'GDSC2': 196770})

# Need to ask Alex, whether this is the correct data for cross-study analysis.
id = np.where(np.isin(res.source, ['CCLE', 'CTRP', 'GDSC1', 'GDSC2', 'gCSI']))[0]
res = res.iloc[id, :]
id = np.where(res.source == 'CTRP')[0]
res.iloc[id, 0] = 'CTRPv2'
id = np.where(res.source == 'GDSC1')[0]
res.iloc[id, 0] = 'GDSCv1'
id = np.where(res.source == 'GDSC2')[0]
res.iloc[id, 0] = 'GDSCv2'

# Load the list of cell lines with gene expressions, mutations, copy numbers, and DNA methylations
ccl_mapping = pd.read_csv(auxiliary_data_dir + 'ccl_ge_mu_cnv_me.txt', sep='\t', engine='c',
                               na_values=['na', '-', ''], header=None, index_col=None, low_memory=False)

# Load drug information
drug_info = pd.read_csv(benchmark_data_dir + 'drug_info.txt', sep='\t', engine='c', na_values=['na', '-', ''],
                        header=0, index_col=None, low_memory=False)

id = np.intersect1d(np.where(np.isin(res.cell, ccl_mapping.iloc[:, 0]))[0],
                    np.where(np.isin(res.drug, drug_info.DrugID))[0])
res = res.iloc[id, :]
# Need to map ANL IDs to DepMap ID and unique drug ID
ccl_mapping.index = ccl_mapping.iloc[:, 0]
drug_info.index = drug_info.DrugID
res.cell = ccl_mapping.loc[res.cell, :].iloc[:, 1].values
res.drug = drug_info.loc[res.drug, :].improve_chem_id.values
res.columns = [str(res.columns[0]), 'improve_sample_id', 'improve_chem_id'] + [str(i) for i in res.columns[3:]]
res.to_csv(benchmark_data_dir + 'response.txt', header=True, index=False, sep='\t', line_terminator='\r\n')

study = np.unique(res.source)
for s in study:
    ids = np.where(res.source == s)[0]
    pd.DataFrame(ids).to_csv(benchmark_data_dir + s + '_all.txt', header=False, index=False, sep='\t',
                             line_terminator='\r\n')
    p = generate_cross_validation_partition(list(range(len(ids))), n_folds=10, n_repeats=1, portions=[8, 1, 1],
                                            random_seed=1)
    for i in range(len(p)):
        pd.DataFrame(ids[p[i][0]]).to_csv(benchmark_data_dir + s + '_split_' + str(i) + '_train.txt', header=False,
                                          index=False, sep='\t', line_terminator='\r\n')
        pd.DataFrame(ids[p[i][1]]).to_csv(benchmark_data_dir + s + '_split_' + str(i) + '_val.txt', header=False,
                                          index=False, sep='\t', line_terminator='\r\n')
        pd.DataFrame(ids[p[i][2]]).to_csv(benchmark_data_dir + s + '_split_' + str(i) + '_test.txt', header=False,
                                          index=False, sep='\t', line_terminator='\r\n')
print(Counter(res.source))
print('Number of cell lines is ' + str(len(np.unique(res.improve_sample_id))))
print('Number of drugs is ' + str(len(np.unique(res.improve_chem_id))))


# Counter({'CTRPv2': 286665, 'GDSCv1': 171940, 'GDSCv2': 114644, 'CCLE': 9519, 'gCSI': 4941})
# Number of cell lines is 785
# Number of drugs is 749
