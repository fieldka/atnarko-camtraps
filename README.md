---
# Title of dataset: Influence of ecotourism on grizzly bear activity depends on salmon abundance in the Atnarko River corridor, Nuxalk Territory
---

R scripts reproduce the results described in 'Influence of ecotourism on grizzly bear activity depends on salmon abundance in the Atnarko River corridor, Nuxalk Territory.' Data generated during this study (i.e., data files described below) are not publicly available due to confidentiality agreement with Nuxalk Nation, but available from the corresponding author and Nuxalk Nation upon reasonable request.

Article citation: Field, K. A., Short,M. L., Moody, J. E., Artelle, K. A., Bourbonnais,M. L., Paquet, P. C., & Darimont, C. T. (2024). Influence of ecotourism on grizzly bear activity depends on salmon abundance in the Atnarko River corridor, Nuxalk Territory. Conservation Science and Practice, e13097.

## Description of the data and scripts:
R Scripts are meant to be run sequentially. Scripts 1 and 2 produce dataframes in the global environment that scripts 3 and 4 call in. 

1 - Data exploration and preparation.R:
Data preparation and exploration; includes script to subset study period of interest, identify camera activity gaps, and calculate diel period metric. 

2 - Co-variate data.R:
Joins co-variate data: Salmon, berries, land-humans, water level.
Produces SI figures. 

3 - GLMM.R:
Weekly Detection analysis, including models with and without imputed salmon data. 

4 - Multinomial.R:
Multinomial regression analysis

Camera trap metadata.xlsx:
Metadata for camera trap variables

atnarko_deployment_data.csv:
Camera trap deployment data, including deployment begin and end dates, as well as install dates. Camera trap settings and models can be found in the SI associated with this manuscript (see also camera trap metadata for full description of variables).

atnarko_detection_data.csv:
Camera trap detection data, including (see camera trap metadata for full description of variables). 

atnarko_station_data.csv:
Camera trap station data, including deployment location ID, install date, station coordinates, treatment (see camera trap metadata for full description of variables). 

berry_surveys.csv: 
Raw berry survey data.

berries_summarized.csv:
Data that summarize the number of shrub species with ripe berries for each berry survey. 

humans.csv:
Human detection data per treatment-week. 

hydro.csv:
Weekly hydrometric data for the Atnarko river near the mouth (see source information below). 

salmon drift counts 2020-2021.csv:
Salmon counts by species.

temporally imputed_salmon drift counts 2020-2021.csv:
Temporally imputed salmon counts as described in methods and in script titled '2 - Co-variate data.R'. 

## Data used from other sources
Salmon escapement data were accessed from the Department of Fisheries and Oceans NuSEDS-New Salmon Escapement Database System. https://open.canada.ca/data/en/dataset/c48669a3-045b-400d-b730-48aafe8c5ee6

Real-Time Hydrometric Data for ATNARKO RIVER NEAR THE MOUTH (08FB006) [BC] were accessed from Government of Canada.
https://wateroffice.ec.gc.ca/report/real_time_e.html?stn=08FB006 

## Code/Software
Analyses were conducted in R.
R Core Team (2023). R: A language and environment for statistical
computing (Version 4.3.1). R Foundation for Statistical Computing, Vienna, Austria.
URL https://www.R-project.org/.
