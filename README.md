# CELLDEX geospatial
Data and code for geospatial modeling of CELLDEX data
Tiegs et al. *Predicting global organic-matter decomposition in flowing waters*.

## List of all files in repository

Raw data files in repository:
  * [`CELLDEX_deploy_dates.csv`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#celldex_deploy_datescsv)
  * [`DRP_Conc.tif`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#drp_conctif-and-no3_conctif)
  * [`NO3_Conc.tif`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#drp_conctif-and-no3_conctif)
  * [`gbc1549-sup-0009-ts02.csv`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#gbc1549-sup-0009-ts02csv)
  * [`gbc20320-sup-0002-supinfo.csv`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#gbc20320-sup-0002-supinfocsv)
  * [`var_names.csv`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#var_namescsv)
  * [`validation_variables.csv`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#validation_variablescsv)
  * [`LeRoy.ExpandedDataset.Kvalues.csv`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#leroyexpandeddatasetkvaluescsv)
  * [`LeafConditionKey2.csv`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#leafconditionkey2csv)
  * [`TRY_traits.csv`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#try_traitscsv)
  * [`Litter_traits_review.csv`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#litter_traits_reviewcsv)


Derived data files in repository:
  * [`litter_processed.csv`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#litter_processedcsv)
  * [`traits.csv`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#traitscsv)
  * `global_kd.tif` Raster of global predictions of ln cotton decomposition rates (1/d) 

Code files in repository:
  * [`CELLDEX_geospatial_BRT.R`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#metadata-for-files-used-in-celldex_geospatial_brtr) 
  code for statistical analyses and figure generation
  * [`litter_process.R`](https://github.com/dmcostello/CELLDEX_geospatial/tree/main#metadata-for-files-used-in-litter_processr) 
  code for cleaning and compiling litter decay and traits

___

## Metadata for files used in `CELLDEX_geospatial_BRT.R`

### `CELLDEX_deploy_dates.csv`

|Parameter     |Definition   |Units  |
| ------------- |-----------| -----|
|partnerid|Numerical partner ID||
|stream|Numerical stream ID||
|deploy_date|Date of cotton deployment|month/day/year|


### `DRP_Conc.tif` and `NO3_Conc.tif`
Source: McDowell et al. (2021) https://doi.org/10.1002/gdj3.111

Repository: https://doi.org/10.25400/lincolnuninz.11894697

Rasters of dissolved reactive phosphorus (DRP) and nitrate (NO<sub>3</sub>) yield from watersheds (kg (NO<sub>3</sub>-N or DRP-P)/ha/yr). Coarse-scale maps are provided in the above source, and fine-scale rasters were provided by authors.


### `gbc1549-sup-0009-ts02.csv`
Source: Mahowald et al. (2008)  https://doi.org/10.1029/2008GB003240

Phosphorus dry deposition (kg P/km<sup>2</sup>/yr). See supporting information at above link for metadata. 

### `gbc20320-sup-0002-supinfo.csv`
Brahney et al. (2015) https://doi.org/10.1002/2015GB005137

Phosphorus dry deposition (kg P/km<sup>2</sup>/yr). See supporting information at above link for metadata. 

### `var_names.csv`
Names of predictor variables that were included in the boosted regression tree model predicting cotton strip decomposition.

|Parameter     |Definition   |
| ------------- |-----------|
|Variables|Name of variable|
|Transform|Transformation needed (see note)|
|Description|Longer description of variable|
|Source|Source of data|

Transformation codes: log = log(x); log1 = log(x+1); xhund = x\*100; xten = x\*10


### `validation_variables.csv`
Names of predictor variables that were included in the boosted regression tree validation model predicting leaf litter decomposition.

|Parameter     |Definition   |
| ------------- |-----------|
|Variables|Name of variable|
|Source|Source of data|


### Other data files stored in online repositories

__Other CELLDEX cotton data__

Source: Tiegs et al. 2019 https://doi.org/10.1126/sciadv.aav0486

Repository: https://github.com/dmcostello/CELLDEX2018 

Files: `CELLDEX_SITE_DATA.csv` `CELLDEX_TEMPERATURE.csv` `str_k.csv`


__Climate and weather__

Source: https://www.worldclim.org/

Repository: `worldclim` package in R

Files: `wc2.1_10m_tmin_12.tif`


__N deposition__

Source: Ackerman et al. (2019) https://doi.org/10.1029/2018GB005990

Repository: https://hdl.handle.net/11299/197613

Files: `inorganic_N_deposition.csv`


__HydroBASINS__

Source: Linke et al. (2019) https://doi.org/10.1038/s41597-019-0300-6

Repository: https://www.hydrosheds.org/products/hydroatlas

Files: `BasinATLAS_v10_lev12.shp`




___

## Metadata for files used in `litter_process.R`

### `LeRoy.ExpandedDataset.Kvalues.csv`

Source: LeRoy et al. (2020) https://doi.org/10.1111/1365-2745.13262

Repository: https://github.com/andrew-hipp/decomposition-phylogeny-2019

See repository for metadata


### `LeafConditionKey2.csv`

Additional information about leaf condition not included in `LeRoy.ExpandedDataset.Kvalues`. Provided by Jenn Follstad Shah. Linked by `Sorting.code`.


### `TRY_traits.csv`

Source: Data request from TRY Plant Trait Database (https://try-db.org/) 

|Parameter     |Definition   |Units  |
| ------------- |-----------| -----|
|Genus|Plant genus||
|N_Leaf_Mn|Mean nitrogen content of fresh leaves|% dry mass|
|N_Leaf_Med|Median nitrogen content of fresh leaves|% dry mass|
|P_Leaf_Mn|Mean phosphorus content of fresh leaves|% dry mass|
|P_Leaf_Med|Median phosphorus content of fresh leaves|% dry mass|
|NtoP_Leaf_Mn|Mean N:P of fresh leaves|Mass ratio|
|NtoP_Leaf_Med|Median N:P of fresh leaves|Mass ratio|
|Thick_Mn|Mean fresh leaf thickness|mm|    
|Thick_Med|Median fresh leaf thickness|mm| 
|C_Leaf_Mn|Mean carbon content of fresh leaves|% dry mass|
|C_Leaf_Med|Median carbon content of fresh leaves|% dry mass|
|CtoN_Leaf_Mn|Mean C:N of fresh leaves|Mass ratio|
|CtoN_Leaf_Med|Median C:N of fresh leaves|Mass ratio|
|Ca_Leaf_Mn|Mean calcium content of fresh leaves|% dry mass|
|Ca_Leaf_Med|Median calcium content of fresh leaves|% dry mass|


### `Litter_traits_review.csv`

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


### `litter_processed.csv`

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


### `traits.csv`

Derived dataset averaging and merging leaf and litter traits by genus. Parameter definitions and units are the same as `TRY_traits.csv` and `Litter_traits_review.csv`. All senesced litter traits are reported as genus-level means. 

