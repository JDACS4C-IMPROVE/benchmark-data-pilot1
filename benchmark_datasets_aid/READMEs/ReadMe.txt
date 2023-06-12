Order of running scripts:
1. get_common_ccl.py
2. generate_omics_data.py and generate_drug_data.py 
3. generate_response_and_split_data.py

Multiomics data of cancer cell lines include:
cancer_copy_number.txt: continuous copy numbers.
cancer_discretized_copy_number.txt: -2, -1, 0, 1, 2 indicate deep deletion, heterozygous deletion, neutral, copy number gain, copy number amplification, respectively.
cancer_DNA_methylation.txt: average DNA methylation values in transcription start sites.
cancer_gene_expression.txt: gene expression values.
cancer_miRNA_expression.txt: miRNA expression values.
cancer_mutation_count.txt: gene-level mutation counts.
cancer_mutation_long_format.txt: mutation data in long format
cancer_mutation.parquet: mutation data in matrix format. Columns are cell lines. Rows are mutation.
cancer_RPPA.txt: protein expressions measured using RPPA.

Multimodal data of drugs include:
drug_descriptor.txt: numeric descriptors of drugs.
drug_fingerprint.txt: binary fingerprints of drugs.
drug_SMILES.txt: SMILES strings of drugs.

response.txt: response data

Take CCLE as an example to explain data partitioning index files. Indices are row indices of response data frame:
CCLE_all.txt: indices of all experiments in CCLE
CCLE_split_2_test.txt: indices of experiments in the testing set of split 2 of CCLE
CCLE_split_2_train.txt: indices of experiments in the training set of split 2 of CCLE 
CCLE_split_2_val.txt: indices of experiments in the validation set of split 2 of CCLE

drug_info.txt provides information of drugs. But it is not used for model evaluation.

