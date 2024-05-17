# CGPA single
This folder contains the CGPA single gene search at the pan-cancer level.



* __app.R__: The main Shiny function for building the single gene search Shiny app.

* __pan_cancer_dashboard.R__: Function for the "Pan-Cancer Summary" tab, including prognostic marker summary, gene expression profile, Kaplan-Meier plots with 
      optimal cutoff, and PPI network (STRING).
* __KM_surv.R__: Corresponds to the "Multivariable Analysis" tab for user-defined cutoffs Kaplan-Meier plots and univariable Cox regression.
* __KM_surv_adj.R__: Corresponds to the "Multivariable Analysis" tab for covariate-adjusted Kaplan-Meier plots and multivariable Cox model.
* __top_prog.R__: Corresponds to the "Top Prognostic Genes" under the "Multivariable Analysis" tab, for checking top prognostic genes and lncRNAs within each 
      cancer type. Users can also explore the subnetwork for top-ranked protein-coding genes.
* __GHI_model.R__: Corresponds to the "Gene-Hallmark Interaction" tab, for testing the interaction between the gene and hallmarks.
   -- code_Splicing.R: Function for the "ProgSplicing" tab.
* __top_prog_across_cancers.R__: Corresponds to the "Multivariable Analysis" tab for user-defined cutoffs Kaplan-Meier plots and univariable Cox regression.
* __TIDE.R__ and * __TIDE_km.R__ : For the "TIDE - OS Interaction" tab under "LncRNA Exploration", using the Tumor Immune Dysfunction and Exclusion (TIDE) 
      computational framework (Jiang et al., 2018).
* __corr_genes.R__: For the "Coexpression Analysis" tab under "LncRNA Exploration", identifying protein-coding genes most correlated with the selected lncRNA.

