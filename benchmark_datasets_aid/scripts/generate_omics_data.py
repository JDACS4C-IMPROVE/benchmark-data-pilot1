import pandas as pd
import numpy as np
from collections import Counter
from Functions import replace_ccl_name



pd.set_option('display.max_columns', None)

omics_data_dir = '../../Data_Curation_final_AID/Curated_CCLE_Multiomics_files/'

benchmark_data_dir = '../CSA_Data/'

auxiliary_data_dir = '../Auxiliary_Data/'



# Load the list of cell lines with gene expressions, mutations, copy numbers, and DNA methylations
ccl_info = pd.read_csv(auxiliary_data_dir + 'ccl_ge_mu_cnv_me.txt', sep='\t', engine='c',
                               na_values=['na', '-', ''], header=None, index_col=None, low_memory=False)
ccl_info.index = ccl_info.iloc[:, 0]



# Load gene expression data
ge = pd.read_csv(omics_data_dir + 'CCLE_AID_expression_full.csv', sep=',', engine='c',
                        na_values=['na', '-', ''], header=None, index_col=None, low_memory=False)
ge = ge.iloc[:, 1:]
ge.iloc[0, 0] = np.nan
id = np.array(range(3, ge.shape[0]))
uni_improve_sample_id, uni_index = np.unique(ge.iloc[id, 0], return_index=True)
ge = ge.iloc[np.concatenate(([0, 1, 2], np.sort(id[uni_index]))), :]
ge.to_csv(benchmark_data_dir + 'cancer_gene_expression.txt', header=False, index=False, sep='\t', line_terminator='\r\n')
ge = None



# Load copy number data
cn = pd.read_csv(omics_data_dir + 'CCLE_AID_gene_cn.csv', sep=',', engine='c',
                        na_values=['na', '-', ''], header=None, index_col=None, low_memory=False)
cn = cn.iloc[:, 1:]
cn.iloc[0, 0] = np.nan
id = np.array(range(3, cn.shape[0]))
uni_improve_sample_id, uni_index = np.unique(cn.iloc[id, 0], return_index=True)
cn = cn.iloc[np.concatenate(([0, 1, 2], np.sort(id[uni_index]))), :]
cn.to_csv(benchmark_data_dir + 'cancer_copy_number.txt', header=False, index=False, sep='\t', line_terminator='\r\n')
cn = None



# Load discretized copy number data
discretized_cn = pd.read_csv(omics_data_dir + 'CCLE_AID_gene_cn_discretized.csv', sep=',', engine='c',
                        na_values=['na', '-', ''], header=None, index_col=None, low_memory=False)
discretized_cn = discretized_cn.iloc[:, 1:]
discretized_cn.iloc[0, 0] = np.nan
id = np.array(range(3, discretized_cn.shape[0]))
uni_improve_sample_id, uni_index = np.unique(discretized_cn.iloc[id, 0], return_index=True)
discretized_cn = discretized_cn.iloc[np.concatenate(([0, 1, 2], np.sort(id[uni_index]))), :]
discretized_cn.to_csv(benchmark_data_dir + 'cancer_discretized_copy_number.txt', header=False, index=False, sep='\t',
          line_terminator='\r\n')
discretized_cn = None



# Load DNA methylation data
me = pd.read_csv(omics_data_dir + 'CCLE_AID_RRBS_TSS_1kb_20180614.csv', sep=',', engine='c',
                        na_values=['na', '-', ''], header=None, index_col=None, low_memory=False)
me.iloc[0, 0] = np.nan
id = np.setdiff1d(range(me.shape[1]), 1)
me = me.iloc[:, id]
id_data = np.array(range(4, me.shape[0]))
id = np.sort(np.where(np.isin(me.iloc[id_data, 0], ccl_info.index))[0])
id_data = id_data[id]
uni_improve_sample_id, uni_index = np.unique(ccl_info.loc[me.iloc[id_data, 0], :].iloc[:, 1], return_index=True)
id_data = id_data[uni_index]
me = me.iloc[np.concatenate(([0, 1, 2, 3], np.sort(id_data))), :]
id_data = np.array(range(4, me.shape[0]))
me.iloc[id_data, 0] = ccl_info.loc[me.iloc[id_data, 0], :].iloc[:, 1].values
me.to_csv(benchmark_data_dir + 'cancer_DNA_methylation.txt', header=False, index=False, sep='\t',
          line_terminator='\r\n')
me = None



# Load miRNA expression data
miRNA = pd.read_csv(omics_data_dir + 'CCLE_AID_miRNA_20180525.csv', sep=',', engine='c',
                        na_values=['na', '-', ''], header=0, index_col=0, low_memory=False)
miRNA.index = miRNA.index.set_names('')
miRNA = miRNA.iloc[:, 1:]
id = np.where(np.isin(miRNA.index, ccl_info.index))[0]
miRNA = miRNA.iloc[id, :]
uni_improve_sample_id, uni_index = np.unique(ccl_info.loc[miRNA.index, :].iloc[:, 1], return_index=True)
miRNA = miRNA.iloc[uni_index, :]
miRNA.index = ccl_info.loc[miRNA.index, :].iloc[:, 1]
miRNA.index = miRNA.index.set_names('')
miRNA.to_csv(benchmark_data_dir + 'cancer_miRNA_expression.txt', header=True, index=True, sep='\t',
          line_terminator='\r\n')
miRNA = None



# Load protein expression data
rppa = pd.read_csv(omics_data_dir + 'CCLE_AID_RPPA_20180123.csv', sep=',', engine='c',
                        na_values=['na', '-', ''], header=0, index_col=0, low_memory=False)
rppa.index = rppa.index.set_names('')
rppa = rppa.iloc[:, 1:]
id = np.where(np.isin(rppa.index, ccl_info.index))[0]
rppa = rppa.iloc[id, :]
uni_improve_sample_id, uni_index = np.unique(ccl_info.loc[rppa.index, :].iloc[:, 1], return_index=True)
rppa = rppa.iloc[uni_index, :]
rppa.index = ccl_info.loc[rppa.index, :].iloc[:, 1]
rppa.index = rppa.index.set_names('')
rppa.to_csv(benchmark_data_dir + 'cancer_RPPA.txt', header=True, index=True, sep='\t',
          line_terminator='\r\n')
rppa = None



# Load mutation count data
mu_count = pd.read_csv(omics_data_dir + 'Mutation_AID_count.csv', sep=',', engine='c',
                        na_values=['na', '-', ''], header=None, index_col=None, low_memory=False)
mu_count = mu_count.iloc[:, 1:]
mu_count.iloc[0, 0] = np.nan
id = np.array(range(3, mu_count.shape[0]))
uni_improve_sample_id, uni_index = np.unique(mu_count.iloc[id, 0], return_index=True)
mu_count = mu_count.iloc[np.concatenate(([0, 1, 2], np.sort(id[uni_index]))), :]
mu_count.to_csv(benchmark_data_dir + 'cancer_mutation_count.txt', header=False, index=False, sep='\t',
          line_terminator='\r\n')
mu_count = None



# Load mutation data
mu = pd.read_parquet(omics_data_dir + 'Mutation_AID_binary.parquet', engine='pyarrow')
id_data = np.array(range(14, mu.shape[1]))
uni_improve_sample_id, uni_index = np.unique(mu.iloc[0, id_data], return_index=True)
id_data = id_data[uni_index]
mu = mu.iloc[:, np.concatenate((list(range(14)), np.sort(id_data)))]
id_data = np.array(range(14, mu.shape[1]))
mu.columns = np.concatenate((mu.columns[:14].values, mu.iloc[0, id_data].values))
mu = mu.iloc[np.setdiff1d(range(mu.shape[0]), 0), :]
mu.to_parquet(benchmark_data_dir + 'cancer_mutation.parquet', engine='pyarrow')
mu = None



# Load_long_format mutation data
mu_long = pd.read_csv(omics_data_dir + 'Mutation_AID_long_format.csv', sep=',', engine='c',
                        na_values=['na', '-', ''], header=0, index_col=None, low_memory=False)
mu_long_com = mu_long.Argonne_ID + '|' + mu_long.DepMap_ID
uni_mu_long_com = np.unique(mu_long_com)
uni_improve_sample_id, uni_index = np.unique([str(i).split('|')[1] for i in uni_mu_long_com], return_index=True)
uni_mu_long_com = uni_mu_long_com[uni_index]
id = np.sort(np.where(np.isin(mu_long_com, uni_mu_long_com))[0])
mu_long = mu_long.iloc[id, ]
mu_long = mu_long.iloc[:, 1:]
mu_long.to_csv(benchmark_data_dir + 'cancer_mutation_long_format.txt', header=True, index=False, sep='\t',
               line_terminator='\r\n')