# RStudio_Projects
RStudio scripts for daily processing

#How to use the MetabolomicsDataAnalysisPipeline

"1_FeatureStatistics"
It contains the instructions and packages for generating the feature (All MS1) statistics for metabolomic studies.
You first have to generate a folder where all the MZmine files will be included:
  MS1.csv: Raw feature quantification matrix, includes all MS1.
  FBMN_quant.csv: Raw feature quantification matrix for MS1 with associated MS1 (for Feature-Based Molecular Networking)
  FBMN.mgf: Spectral file associated with FBMN_quant.csv (for Feature-Based Molecular Networking)
  SIRIUS.mgf: Spectral file with isotopic information (For SIRIUS)
