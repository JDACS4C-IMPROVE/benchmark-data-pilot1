from pathlib import Path
import pandas as pd
import numpy as np



fdir = Path(__file__).reslove().parent
# pd.set_option('display.max_columns', None)

# omics_data_dir = '../../Data_Curation_final_AID/Curated_CCLE_Multiomics_files/'
omics_data_dir = fdir/'../../Data_Curation_final_AID/Curated_CCLE_Multiomics_files/'
# auxiliary_data_dir = '../Auxiliary_Data/'
auxiliary_data_dir = fdir/'../auxiliary_data/'

# Load gene expression data
ge = pd.read_csv(omics_data_dir + '/CCLE_AID_expression_full.csv', sep=',', engine='c', na_values=['na', '-', ''],
                 header=None, index_col=None, low_memory=False)
ge_ccl = ge.iloc[3:, :2]
ge = None

# Load mutation count data
mu_count = pd.read_csv(omics_data_dir + '/Mutation_AID_count.csv', sep=',', engine='c', na_values=['na', '-', ''],
                       header=None, index_col=None, low_memory=False)
mu_ccl = mu_count.iloc[3:, :2]
mu_count = None

# Load copy number data
cnv = pd.read_csv(omics_data_dir + '/CCLE_AID_gene_cn_discretized.csv', sep=',', engine='c', na_values=['na', '-', ''], header=None,
                  index_col=None, low_memory=False)
cnv_ccl = cnv.iloc[3:, :2]
cnv = None

mapping = pd.concat((ge_ccl, mu_ccl, cnv_ccl), axis=0)
mapping = mapping.drop_duplicates()
ccl_ge_mu_cnv = np.intersect1d(np.intersect1d(ge_ccl.iloc[:, 1], mu_ccl.iloc[:, 1]), cnv_ccl.iloc[:, 1])
mapping = mapping.iloc[np.where(np.isin(mapping.iloc[:, 1], ccl_ge_mu_cnv))[0], ]
pd.DataFrame(mapping).to_csv(auxiliary_data_dir + 'ccl_ge_mu_cnv.txt', header=False, index=False, sep='\t',
                                   line_terminator='\r\n')

# Load DNA methylation data
me = pd.read_csv(omics_data_dir + '/CCLE_AID_RRBS_TSS_1kb_20180614.csv', sep=',', engine='c', na_values=['na', '-', ''],
                 header=None, index_col=0, low_memory=False)
me_ccl = me.index[4:]
me = None
mapping = mapping.iloc[np.where(np.isin(mapping.iloc[:, 0], me_ccl))[0], ]
pd.DataFrame(mapping).to_csv(auxiliary_data_dir + 'ccl_ge_mu_cnv_me.txt', header=False, index=False, sep='\t',
                                   line_terminator='\r\n')

