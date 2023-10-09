# CELLDEX geospatial
Data and code for geospatial modeling of CELLDEX data
Tiegs et al. *Predicting global organic-matter decomposition in flowing waters*.

## List of all files in repository

Raw data files in repository:
  * [`CELLDEX_deploy_dates.csv`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#celldex_deploy_datescsv)
  * `DRP_Conc.tif`
  * `NO3_Conc.tif`
  * `gbc1549-sup-0009-ts02.csv`
  * `gbc20320-sup-0002-supinfo.csv`
  * `var_names.csv`
  * `validation_variables.csv`
  * `LeRoy.ExpandedDataset.Kvalues.csv`
  * `LeafConditionKey2.csv`
  * `TRY_traits.csv`
  * `Litter_traits_review.csv`


Derived data files in repository:
  * `litter_processed.csv`
  * `traits.csv`


Code files in repository:
  * [`CELLDEX_geospatial_BRT.R`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#metadata-for-files-used-in-celldex_geospatial_brtr) 
  code for statistical analyses and figure generation
  * [`litter_process.R`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#metadata-for-files-used-in-litter_processr) 
  code for cleaning and compiling litter decay and traits

___

## Metadata for files used in `CELLDEX_geospatial_BRT.R`

### __`CELLDEX_deploy_dates.csv`__

|Parameter     |Definition   |Units  |
| ------------- |-----------| -----|
|partnerid|Numerical partner ID||
|stream|Numerical stream ID||
|deploy_date|Date of cotton deployment|month/day/year|


### __`DRP_Conc.tif`__ and __`NO3_Conc.tif`__
*Source*: McDowell et al. (2021) https://doi.org/10.1002/gdj3.111

*Repository*: https://doi.org/10.25400/lincolnuninz.11894697

Rasters of dissolved reactive phosphorus (DRP) and nitrate (NO3) yield from watersheds (kg (NO~3~-N or DRP-P)/ha/yr). Coarse-scale maps are provided in the above source, and fine-scale rasters were provided by authors.


### `gbc1549-sup-0009-ts02.csv`
*Source*: Mahowald et al. (2008)  https://doi.org/10.1029/2008GB003240

Phosphorus dry deposition (kg P/km^2^/yr). See supporting information at above link for metadata. 

### `gbc20320-sup-0002-supinfo.csv`
Brahney et al. (2015) https://doi.org/10.1002/2015GB005137

Phosphorus dry deposition (kg P/km^2^/yr). See supporting information at above link for metadata. 

  * `var_names.csv`
  * `validation_variables.csv`

___

## Metadata for files used in `litter_process.R`

### __`LeRoy.ExpandedDataset.Kvalues.csv`__

*Source*: LeRoy et al. (2020) https://doi.org/10.1111/1365-2745.13262

*Repository*: https://github.com/andrew-hipp/decomposition-phylogeny-2019

See repository for metadata


### __`LeafConditionKey2.csv`__

Additional information about leaf condition not included in `LeRoy.ExpandedDataset.Kvalues`. Provided by Jenn Follstad Shah. Linked by `Sorting.code`.


### __`TRY_traits.csv`__

*Source*: Data request from TRY Plant Trait Database (https://try-db.org/) 

|Parameter     |Definition   |Units  |
| ------------- |-----------| -----|
|Genus|Plant genus||
|N_Leaf_Mn|Mean nitrogen content of fresh leaves|% dry mass|
|N_Leaf_Med|Median nitrogen content of fresh leaves|% dry mass|
|P_Leaf_Mn|Mean phosphorus content of fresh leaves|% dry mass|
|P_Leaf_Med|Median phosphorus content of fresh leaves|% dry mass|
|NtoP_Leaf_Mn|Mean N:P of fresh leaves||
|NtoP_Leaf_Med|Median N:P of fresh leaves||
|Thick_Mn|Mean fresh leaf thickness|mm|    
|Thick_Med|Median fresh leaf thickness|mm| 
|C_Leaf_Mn|Mean carbon content of fresh leaves|% dry mass|
|C_Leaf_Med|Median carbon content of fresh leaves|% dry mass|
|CtoN_Leaf_Mn|Mean C:N of fresh leaves||
|CtoN_Leaf_Med|Median C:N of fresh leaves||
|Ca_Leaf_Mn|Mean calcium content of fresh leaves|% dry mass|
|Ca_Leaf_Med|Median calcium content of fresh leaves|% dry mass|


### __`Litter_traits_review.csv`__

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


### __`litter_processed.csv`__

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

1. Data were screened to include only coarse or fine mesh bags.
2. Data were screened to include only senesced or green leaves. Air-dried leaves were classified as senesced.


### __`traits.csv`__

Derived dataset averaging and merging leaf and litter traits by genus. Parameter definitions and units are the same as `TRY_traits.csv` and `Litter_traits_review.csv`. All senesced litter traits are reported as genus-level means. 

