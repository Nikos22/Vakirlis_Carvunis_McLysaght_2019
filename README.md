This repository contains all data necessary to reproduce all figures and statistics of Vakirlis, Carvunis and McLysaht 2019 .

The scripts can be found in the scripts/ folder and the figures to which they relate should be evident from the names of the files.
The figures, which are the output of the scripts, can be found in the figures/ folder. 

Figure5B.csv is an extended version of Figure 5B of the article in text table format.

All source data that are referenced in the article can be found here:

Figure_3-source_data_1.csv : Data on undetectable homologies for different E-value cut-offs used to generate the top panel of Figure 3A.
Column names: "total" : total number of genes in conserved micro-synteny, "not_found" : number of genes without significant sequence similarity, "div" : time since divergence from focal species

Figure_3-source_data_2.csv : Data on false homologies for different E-value cut-offs used to generate the bottom panel of Figure 3A.
Column names: "found" : number of genes with significant sequence similarity, The rest are the same as in the file above. 

Figure_5_source_data_1.csv : dN and dS data used to generate Figure 5 and the accompanying stats. See Methods section for how these data were generated.
Column names: "micro-synteny" : whether the gene satisfies our conserved micro-synteny criteria with the relevant species (see Methods).

Figure_6-source_data_1.xls : An excel file with one dataset per sheet, containing the similarity and micro-synteny conservation information for every focal - target species comparison. The tables can also be found separately in text format in the synt_simil_tables/ folder.

Figure_7_source_data_1.csv : Data on common Pfam matches and gene/protein properties used to generate Figure 5.
Column names: "Gene_focal" : Name of the focal species gene, "Gene_ortho" : name of the target species gene, "same" : whether a common Pfam match was found or these genes are in the same OrthoDB group, "focal" : the value for the property in the focal gene, "ortho" : the value for the property in the target species gene, "var" : the name of the property

The table figure supplements are also provided here, details can be found in the article supplementary material.

Furthermore, we provide some additional/raw data used by the scripts:

all_gene_pairs/ : Data (for each dataset separately and combined in one file) for pairs of undetectable homologues used to calculate correlations after removal of pairs with high percentage of low complexity.
Column names should be self-explanatory.

Pfam_search_raw_data/ : PfamScan search output files for focal and target species proteins of interest.

synt_simil_tables/ : Individual data tables that make up Figure_6-source_data_1.xls . See readme inside the folder. 

divergence_times/ : Divergence times from the focal species for each target species.

dnds/ : Raw data used for the dN, dS analyses.

full_lists_genes/ : lists with gene names and IDs used in the analyses.

nr_results/ : Parsed results of similarity searches in NCBI's NR database, see Methods for details.

oo_dfs/ : Raw files containing all the focal-target gene pairs found in conserved micro-synteny.







  
