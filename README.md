# SMR-Antipsychotics-VTE
This project is the script required to implement TwoSample Mendelian randomisation of drug target gene expression and venous thromboembolism.
Methods of obtaining the blood expression data:
# Input files:
# 1. eQTL data in SMR query format (provided in the smr_query_files folder)
# 2. GWAS data gcta_format
# 3. List of genes for which analysis needs to be done (ENSG_IDs_for_drug_target_genes.txt)
# EXTRACT EQTL SUMMARY DATA FOR GENES OF INTEREST USING SMR  
## NOTE: eqtl summary files are provided in the smr_query_files folder 
# SMR command
system("smr --beqtl-summary cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense --query 1.0e-5 --gene gene_id --out myquery
")
## where myeqtl_data is the summary-level data from a eQTL study in smr format which can be downloaded from https://cnsgenomics.com/software/smr/#DataResource;
## --gene is the gene_id as provided in the eqtl data.
Brain tissue eQTL data were downloaded directly from the GTEx webpage.
In the SMR section for single gene expression and venous thromboembolism associations, only one is listed as an example (each_gene_example_ADRA2A.txt) due to the large number of individual gene expression files.
