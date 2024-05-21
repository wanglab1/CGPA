# CGPA single
This folder contains the CGPA single gene search at the pan-cancer level.



* __CGPA_single_app.R__: The main Shiny function for building the single gene search Shiny app.

* __Pan-cancer_summary.R__: Function for the "Pan-Cancer Summary" tab, including prognostic marker summary, gene expression profile, Kaplan-Meier plots with 
      optimal cutoff, and PPI network (STRING).
* __Multivariable_Analysis_uni.R__: Corresponds to the "Multivariable Analysis" tab for user-defined cutoffs Kaplan-Meier plots and univariable Cox regression.
* __Multivariable_Analysis_multi.R__: Corresponds to the "Multivariable Analysis" tab for covariate-adjusted Kaplan-Meier plots and multivariable Cox model.
* __Multivariable_Analysis_Top_prognostic_genes.R__: Corresponds to the "Top Prognostic Genes" under the "Multivariable Analysis" tab, for checking top prognostic genes and lncRNAs within each 
      cancer type. Users can also explore the subnetwork for top-ranked protein-coding genes.
* __Gene-hallmark_interation.R__: Corresponds to the "Gene-Hallmark Interaction" tab, for testing the interaction between the gene and hallmarks.
   -- code_Splicing.R: Function for the "ProgSplicing" tab.
* __lncRNA_top_prog_across_cancers.R__: Corresponds to the "Multivariable Analysis" tab for user-defined cutoffs Kaplan-Meier plots and univariable Cox regression.
* __lncRNA_TIDE.R__ and * __TlncRNA_TIDE_km.R__ : For the "TIDE - OS Interaction" tab under "LncRNA Exploration", using the Tumor Immune Dysfunction and Exclusion (TIDE) computational framework (Jiang et al., 2018).
* __lncRNA_exploration_coexpression.R__: For the "Coexpression Analysis" tab under "LncRNA Exploration", identifying protein-coding genes most correlated with the selected lncRNA.

