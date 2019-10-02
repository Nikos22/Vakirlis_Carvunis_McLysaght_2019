This repository contains all data to reproduce all figures and statistics of Vakirlis, Carvunis and McLysaht 2019 .


Figure_3A1_data.csv : Data on undetectable homologies for different E-value cut-offs used to generate the top panel of Figure 3A.
Column names: "total" : total number of genes in conserved micro-synteny, "not_found" : number of genes without significant sequence similarity, "div" : time since divergence from focal species

Figure_3A2_data.csv : Data on false homologies for different E-value cut-offs used to generate the bottom panel of Figure 3A.
Column names: "found" : number of genes with significant sequence similarity, The rest are the same as in the file above. 

Supp_figure_4_data.csv : dN, dS and dN/dS data used to generate Supplementary Figure 4 and the accompanying stats. See Methods section for how these data were generated.
Column names: "micro-synteny" : whether the gene satisfies our conserved micro-synteny criteria with the relevant species (see Methods).

Figure_5_data.csv : Data on common Pfam matches and gene/protein properties used to generate Figure 5.
Column names: "Gene_focal" : Name of the focal species gene, "Gene_ortho" : name of the target species gene, "Common Pfam match or OrthoDB group" : whether a common Pfam match was found
or these genes are in the same OrthoDB group, "value focal" : the value for the property in the focal gene, "value target" : the value for the property in the target species gene,
"property" : the name of the property

All_undetectable_homologue_pairs_data.csv : Data for pairs of undetectable homologues used to calculate correlations after removal of pairs with high percentage of low complexity.
Column names should be self-explanatory.

Supplementary_Table_1_for_paper.csv : Supp. Table 1 from the paper.

Supp_table_3.csv : Supp. Table 3 of the paper.

synt_simil_tables/ : data for Figure 4 and  Supp. Figure 6, see readme in the folder

Pfam_search_raw_data/ : PfamScan search output files for focal and target species proteins of interest



  
