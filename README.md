# CELLDEX geospatial
Data and code for geospatial modeling of CELLDEX data
Tiegs et al. *Predicting global organic-matter decomposition in flowing waters*.

## List of all files in repository

Raw data files in repository:\
  * `LeRoy.ExpandedDataset.Kvalues.csv`\
  * `LeafConditionKey2.csv`\
  * `TRY_traits.csv`\
  * `Litter_traits_review.csv`\
  * `CELLDEX_deploy_dates.csv`\
  * `DRP_Conc.tif`\
  * `NO3_Conc.tif`\
  * `gbc1549-sup-0009-ts02.csv`\
  * `gbc20320-sup-0002-supinfo.csv`\
  * `var_names.csv`\
  * `validation_variables.csv`\


Derived data files in repository:\
  * `litter_processed.csv`\
  * `traits.csv`


Code files in repository:\
  * `CELLDEX_geospatial_BRT.R` code for statistical analyses and figure generation\
  * `litter_process.R` code for cleaning and compiling litter decay and traits

___

## Metadata for files used in `litter_process.R`

__`LeRoy.ExpandedDataset.Kvalues.csv`__

*Source*: LeRoy et al. (2020) https://doi.org/10.1111/1365-2745.13262

*Repository*: https://github.com/andrew-hipp/decomposition-phylogeny-2019

See repository for metadata


__`LeafConditionKey2.csv`__

Additional information about leaf condition not included in `LeRoy.ExpandedDataset.Kvalues`. Provided by Jenn Follstad Shah. Linked by `Sorting.code`.


__`TRY_traits.csv`__

*Source*: Data request from TRY Plant Trait Database (https://try-db.org/) 

|Parameter     |Definition   |Units  |
| ------------- |-----------| -----|
|Genus|Plant genus||
|N_Leaf_Mn|Mean nitrogen content of fresh leaves|% dry mass|
|N_Leaf_Med|Median nitrogen content of fresh leaves|% dry mass|
|P_Leaf_Mn|Mean phosphorus content of fresh leaves|% dry mass|
|P_Leaf_Med|Median phosphorus content of fresh leaves|% dry mass|
|NtoP_Leaf_Mn|Mean N:P of fresh leaves|molar|
|NtoP_Leaf_Med|Median N:P of fresh leaves|molar|
|Thick_Mn|Mean fresh leaf thickness||    
|Thick_Med|Median fresh leaf thickness|| 
|C_Leaf_Mn|Mean carbon content of fresh leaves|% dry mass|
|C_Leaf_Med|Median carbon content of fresh leaves|% dry mass|
|CtoN_Leaf_Mn|Mean C:N of fresh leaves|molar|
|CtoN_Leaf_Med|Median C:N of fresh leaves|molar|
|Ca_Leaf_Mn|Mean calcium content of fresh leaves|% dry mass|
|Ca_Leaf_Med|Median calcium content of fresh leaves|% dry mass|


__`Litter_traits_review.csv`__

Database of litter chemistry from stream decomp experiments from literature. Provided by Jenn Follstad Shah.

|Parameter     |Definition   |Units  |
| ------------- |-----------| -----|
|Citation|Citation from which data were extracted||
|Sci_name|Taxonomic name of litter used||
|Genus|Taxonomic genus||
|Species|Taxonomic species||
|Family|Taxonomic family||
|Stream_names|Name of study stream||
|perC|Litter carbon content|% dry mass|
|perN|Litter nitrogen content|% dry mass|
|perP|Litter phosphorus content|% dry mass|
|CtoN|Litter C:N|Molar ratio|
|CtoP|Litter C:P|Molar ratio|
|NtoP|Litter N:P|Molar ratio|
|How_ratios_were_determined| If ratios were reported in paper or calculated||
|perLignin|Litter lignin content|% dry mass|
|perCellulose|Litter cellulose content|% dry mass|


__`litter_processed.csv`__

Derived dataset of average litter decomposition rates in streams with at least 3 measurements for each genera.

|Parameter     |Definition   |Units  |
| ------------- |-----------| -----|
|Sorting.code|Code to link back to original data||
|Mesh.size|Size of mesh in litter bag|See note 1| 
|Leaf.condition|Condition of leaves in litter bag|See note 2|
|Genus|Taxonomic genus||
|latitude|Reported latitude|Decimal degrees|
|longitude|Reported longitude|Decimal degrees|
|mean_kd|Mean decomposition rate|1/d|
|mean_mean_daily_temp|Mean reported average temperature|°C|
|Citation|Citation from which data were extracted||

1. Data were screened to include only coarse or fine mesh bags.\
2. Data were screened to include only senesced or green leaves. Air-dried leaves were classified as senesced.


  * `traits.csv`
