# CELLDEX_geospatial
Data and code for geospatial modeling of CELLDEX data
Tiegs et al. *Predicting global organic-matter decomposition in flowing waters*.


Raw data files in repository:
  * `LeRoy.ExpandedDataset.Kvalues.csv`
  * `LeafConditionKey2.csv`
  * `TRY_traits.csv`
  * `Litter_traits_review.csv`
  * `CELLDEX_deploy_dates.csv`
  * `DRP_Conc.tif`
  * `NO3_Conc.tif`
  * `gbc1549-sup-0009-ts02.csv`
  * `gbc20320-sup-0002-supinfo.csv`
  * `var_names.csv`
  * `validation_variables.csv`


Derived data files in repository:
  * `litter_processed.csv`
  * `traits.csv`


Code files in repository:
  * `CELLDEX_geospatial_BRT.R` code for statistical analyses and figure generation
  * `litter_process.R` code for cleaning and compiling litter decay and traits

Metadata:

__LeRoy.ExpandedDataset.Kvalues.csv__

*Source*: LeRoy et al. (2020) https://doi.org/10.1111/1365-2745.13262

*Repository*: https://github.com/andrew-hipp/decomposition-phylogeny-2019

See repository for metadata

__`TRY_traits.csv`__

*Source*: Data request from TRY Plant Trait Database (https://try-db.org/) 

|Parameter     |Definition   |Units  |
| ------------- |-----------| -----|
|Genus|Plant genus||
|N_Leaf_Mn|Mean N content of fresh leaves|% dry mass|
|N_Leaf_Med|Median N content of fresh leaves|% dry mass|
|P_Leaf_Mn|Mean P content of fresh leaves|% dry mass|
|P_Leaf_Med|Median P content of fresh leaves|% dry mass|
|NtoP_Leaf_Mn|Mean N:P of fresh leaves|molar|
|NtoP_Leaf_Med|Median N:P of fresh leaves|molar|
|Thick_Mn|Mean fresh leaf thickness||    
|Thick_Med|Median fresh leaf thickness|| 
|C_Leaf_Mn|Mean C content of fresh leaves|% dry mass|
|C_Leaf_Med|Median C content of fresh leaves|% dry mass|
|CtoN_Leaf_Mn|Mean C:N of fresh leaves|molar|
|CtoN_Leaf_Med|Median C:N of fresh leaves|molar|
|Ca_Leaf_Mn|Mean Ca content of fresh leaves|% dry mass|
|Ca_Leaf_Med|Median Ca content of fresh leaves|% dry mass|
