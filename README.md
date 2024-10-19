# RStudio_Projects

RStudio scripts for daily processing.

## How to use the Metabolomics Data Analysis Pipeline

### MetabolomicsDataAnalysisPipeline_1_FeatureStatistics.R

It contains the instructions and packages for generating the feature (All MS1) statistics for untargeted metabolomics studies.<br>
You first have to generate a folder where the quantification file from MZmine (v2.53 was consiered for this pipeline) will be inside:<br>
  * **MS1.csv**: Raw feature quantification matrix, includes all MS1.<br>
In this section, MetaboAnalyst web server is used for multivariate statistics (PCA, heatmap, etc.), and NormalyzerDE for univariate statistics.<br>

This pipeline generates:
 * A *NormalyzerDE project folder*.<br>
 * A normalized, imputed quantification table as **MetaboAnalyst.csv** file for multivariate statistics in MetaboAnalyst.<br>
 * **Volcano plot** depicting metabolites with altered abundance between groups.<br>

### MetabolomicsDataAnalysisPipeline_2_ChemicalClassStatistics.R

Contains the instructions and packages for generating graphs and plots related to chemical class annotations statistics from SIRIUS.<br>
First, create a folder with the following project folders:<br>
 * The **sirius project** folder, with the **canopus_formula_summary.tsv** file from SIRIUS.<br>
 * The **DiffExp** folder, with the **DiffExp_stats.tsv** file from the *NormalyzerDE project folder*.<br>

 This pipeline generates:
  * **Bar plots** illustrating *chemical classes with differential abundance* between groups.<br>
  * **Violin plots** representing the chemical diversity in each group through the *Shannon index*.<br>
  * **Hypothesis-testing results (_P_-values)** for chemical diversity differences.<br>
  * **Bar plots** showing the *total abundance of chemical classes in each group*, both absolute and *log2* transformed.<br>

  
