# DNM1

This repository contains a script which reflects the analysis and data visualization performed to assess exon 10-related differences in disease-causing variants of *DNM1*-related disorders. 

Required files/inputs are specified in the scripts. 

For transcript quantification analysis, the raw quantities input file is obtained by: 
1) performing transcript quantification for each sample using Kallisto
2) subsetting the output for *DNM1* transcripts as per Ensembl
3) merging/row-binding all samples into a single table 
Additional sample information tied to sample IDs can be included in a separate table, such as found in Table S1. 

Transcript quantification is supplemented by comparing with adult control samples obtained from the GTEx Portal v8, https://gtexportal.org/home/datasets. 
The GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct file is subsetted for gene_id = 'ENSG00000106976.20' for DNM1 expression. 
Sample-level information is taken from the GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt file to filter by tissue type. 

Splicing and intronic variants are analyzed with output from SpliceAI, https://spliceailookup.broadinstitute.org. 

Missense variation in the population is analyzed using gnomAD data. 
Transcript-specific constraint metrics are downloaded from Supplementary Dataset 11 of Karczewski et al., doi:10.1038/s41586-020-2308-7. 
*DNM1* missense variants are downloaded as an export from the gnomAD browser, https://gnomad.broadinstitute.org/gene/ENSG00000106976. 
