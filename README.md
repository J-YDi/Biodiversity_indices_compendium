# QUANTIFY PHYTOPLANKTON BIODIVERSITY
## Dias Jean-Yves $^1$, Vincent Dorothée $^2$, Goberville Eric $^1$
#### $^1$ Laboratoire de Biologie des Organismes et des Écosystèmes Aquatiques-BOREA, Muséum national d’Histoire naturelle (MNHN), SU, CNRS, IRD, UA, F-75005 Paris, France ; $^2$ PatriNat (OFB, MNHN), Brest, France ; Office Français de la Biodiversité
##### Correspond to the R scripts and data for the INDIBIO project report (2026) and Chapter 1 of Jean-Yves Dias' PhD thesis (2025-2028)

#### Github repository organization

##### script folder :
+ 01_Listing_functions_alpha_indices.R : Script to create dataset of values from functions that calculate alpha diversity indices
+ 02_Discriminate_alpha_functions.R : Script to discrimate functions that return a false value of alpha diverisity indices
+ 03_Alpha_index_analysis.R : Script to analyze alpha biodiversity indices
+ 04_Listing_functions_beta_indices.R : Script to create dataset of values from functions that calculate beta diversity indices
+ 05_Discriminate_beta_functions.R : Script to discrimate functions that return a false value of beta diverisity indices
+ 06_Beta_index_analysis.R : Script to analyze beta biodiversity indices
+ 07_Calculate_beta_indices_transformed_data.R : Script to create dataset of values from functions that calculate beta diversity indices with data transformation
+ 07_Compare_transformations_beta_indices.R : Script to represent the influence on transformations on beta biodiversity indices
+ 08_Bibliometry.R : Script to represent the number of publications involving biodiversity and/or marine and/or policy

##### data folder : 
+ contains all the original raw datasets and additionnal files supporting the data processing.

##### output folder : 
+ contains all the figures generated

 Listing_fonctions_R_Indices_div.xlsx : it contains an Excel-formatted list of all functions for computing alpha and beta diversity, including those that return erroneous values; it explains why when possible, indicates whether beta indices calculate dissimilarity or similarity, and also lists the packages.

###### Contact : jean-yves.dias@sorbonne-universite.fr
