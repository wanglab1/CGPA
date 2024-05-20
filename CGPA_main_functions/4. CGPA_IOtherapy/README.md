# CGPA IOtherapy
This folder contains the CGPA IOtherapy.

* __CGPA_IOtherapy_app.R__: The main Shiny function for building the CGPA IOtherapy Shiny app.

* __Single-gene_panCancer_summary.R__: Function for the "Pan-Cancer Summary" tab, including prognostic marker summary, gene expression profile, Kaplan-Meier plots with different cutoffs, and PPI network (STRING).
* __Single-gene_univariable_analysis.R__: Corresponds to the "Multivariable Analysis" tab for user-defined cutoffs Kaplan-Meier plots and univariable Cox regression.
* __Single-gene_multivariable_analysis.R__: Corresponds to the "Multivariable Analysis" tab for covariate-adjusted Kaplan-Meier plots and multivariable Cox model.
* __Single-gene_GeneHallmark_interaction.R__: Corresponds to the "Gene-Hallmark Interaction" tab, for testing the interaction between the gene and hallmarks.
