
<!-- README.md is generated from README.Rmd. Please edit that file -->

# S2MET

## Description

Project repository for the spring two-row (**S2**) barley
multi-environment trial (**MET**) project. This project will look at the
impact of genotype-by-environment interaction when performing genomewide
prediction. The project includes a panel of 183 training population
lines, 50 validation population lines, and 9 repeated checks.

## Script Order

This is the order in which scripts should be executed to replicate the
results. All scripts are located in the
[Scripts](https://github.com/neyhartj/S2MET_Predictions/tree/master/Scripts)
directory.

Phenotypic data  
1\. `phenotype_data_adjustment.R` - Calculate genotypic means and
calculate heritability in each environment.  
2\. `phenotype_data_summary.R` - Summarize phenotypic data across all
environments, including distributions and heritability.

Environmental variables  
1\. `collect_environmental_variables.R` - Collect environmental
covariate data from NOAA and Soil Survey databases.  
2\. `manipulate_analyze_environmental_variables.R` - Calculate summary
variables, including growing degree days (GDD), accumulated GDD (AGDD),
and photothermality.  
3\. `model_environmental_variables.R` - Find the variables that are
highly correlated with the environmental mean and create covariance and
distance matrices.

Clustering  
1\. `environmental_clustering.R` - Calculate great circle distance and
phenotypic distance, then create hierarchical clusters based on the
distance objects.

Prediction 1.

## Software version

*R* version 3.5.3 was used to complete analyses in this project.

The following packages (with version information) were used to complete
this project:

<table>

<thead>

<tr>

<th style="text-align:left;">

package

</th>

<th style="text-align:left;">

version

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

boot

</td>

<td style="text-align:left;">

1.3-20

</td>

</tr>

<tr>

<td style="text-align:left;">

neyhart

</td>

<td style="text-align:left;">

0.0.0.9000

</td>

</tr>

<tr>

<td style="text-align:left;">

car

</td>

<td style="text-align:left;">

3.0-2

</td>

</tr>

<tr>

<td style="text-align:left;">

effects

</td>

<td style="text-align:left;">

4.1-0

</td>

</tr>

<tr>

<td style="text-align:left;">

carData

</td>

<td style="text-align:left;">

3.0-2

</td>

</tr>

<tr>

<td style="text-align:left;">

data.table

</td>

<td style="text-align:left;">

1.12.2

</td>

</tr>

<tr>

<td style="text-align:left;">

ggrepel

</td>

<td style="text-align:left;">

0.8.0

</td>

</tr>

<tr>

<td style="text-align:left;">

lme4qtl

</td>

<td style="text-align:left;">

0.1.10

</td>

</tr>

<tr>

<td style="text-align:left;">

rrBLUP

</td>

<td style="text-align:left;">

4.6

</td>

</tr>

<tr>

<td style="text-align:left;">

ggridges

</td>

<td style="text-align:left;">

0.5.1

</td>

</tr>

<tr>

<td style="text-align:left;">

readxl

</td>

<td style="text-align:left;">

1.3.1

</td>

</tr>

<tr>

<td style="text-align:left;">

gridExtra

</td>

<td style="text-align:left;">

2.3

</td>

</tr>

<tr>

<td style="text-align:left;">

flextable

</td>

<td style="text-align:left;">

0.5.2

</td>

</tr>

<tr>

<td style="text-align:left;">

officer

</td>

<td style="text-align:left;">

0.3.3

</td>

</tr>

<tr>

<td style="text-align:left;">

kableExtra

</td>

<td style="text-align:left;">

1.1.0

</td>

</tr>

<tr>

<td style="text-align:left;">

pls

</td>

<td style="text-align:left;">

2.7-1

</td>

</tr>

<tr>

<td style="text-align:left;">

optimx

</td>

<td style="text-align:left;">

2018-7.10

</td>

</tr>

<tr>

<td style="text-align:left;">

lmerTest

</td>

<td style="text-align:left;">

3.1-0

</td>

</tr>

<tr>

<td style="text-align:left;">

measurements

</td>

<td style="text-align:left;">

1.3.0

</td>

</tr>

<tr>

<td style="text-align:left;">

rvest

</td>

<td style="text-align:left;">

0.3.3

</td>

</tr>

<tr>

<td style="text-align:left;">

xml2

</td>

<td style="text-align:left;">

1.2.0

</td>

</tr>

<tr>

<td style="text-align:left;">

raster

</td>

<td style="text-align:left;">

2.8-19

</td>

</tr>

<tr>

<td style="text-align:left;">

sp

</td>

<td style="text-align:left;">

1.3-1

</td>

</tr>

<tr>

<td style="text-align:left;">

foreign

</td>

<td style="text-align:left;">

0.8-71

</td>

</tr>

<tr>

<td style="text-align:left;">

soilDB

</td>

<td style="text-align:left;">

2.3.5

</td>

</tr>

<tr>

<td style="text-align:left;">

aqp

</td>

<td style="text-align:left;">

1.17

</td>

</tr>

<tr>

<td style="text-align:left;">

geosphere

</td>

<td style="text-align:left;">

1.5-7

</td>

</tr>

<tr>

<td style="text-align:left;">

rnoaa

</td>

<td style="text-align:left;">

0.8.4

</td>

</tr>

<tr>

<td style="text-align:left;">

modelr

</td>

<td style="text-align:left;">

0.1.4

</td>

</tr>

<tr>

<td style="text-align:left;">

cowplot

</td>

<td style="text-align:left;">

0.9.4

</td>

</tr>

<tr>

<td style="text-align:left;">

ggdendro

</td>

<td style="text-align:left;">

0.1-20

</td>

</tr>

<tr>

<td style="text-align:left;">

agridat

</td>

<td style="text-align:left;">

1.16

</td>

</tr>

<tr>

<td style="text-align:left;">

pbr

</td>

<td style="text-align:left;">

0.1.0

</td>

</tr>

<tr>

<td style="text-align:left;">

cluster

</td>

<td style="text-align:left;">

2.0.8

</td>

</tr>

<tr>

<td style="text-align:left;">

mclust

</td>

<td style="text-align:left;">

5.4.3

</td>

</tr>

<tr>

<td style="text-align:left;">

lubridate

</td>

<td style="text-align:left;">

1.7.4

</td>

</tr>

<tr>

<td style="text-align:left;">

broom

</td>

<td style="text-align:left;">

0.5.2

</td>

</tr>

<tr>

<td style="text-align:left;">

lme4

</td>

<td style="text-align:left;">

1.1-21

</td>

</tr>

<tr>

<td style="text-align:left;">

Matrix

</td>

<td style="text-align:left;">

1.2-17

</td>

</tr>

<tr>

<td style="text-align:left;">

forcats

</td>

<td style="text-align:left;">

0.4.0

</td>

</tr>

<tr>

<td style="text-align:left;">

stringr

</td>

<td style="text-align:left;">

1.4.0

</td>

</tr>

<tr>

<td style="text-align:left;">

dplyr

</td>

<td style="text-align:left;">

0.8.0.1

</td>

</tr>

<tr>

<td style="text-align:left;">

purrr

</td>

<td style="text-align:left;">

0.3.2

</td>

</tr>

<tr>

<td style="text-align:left;">

readr

</td>

<td style="text-align:left;">

1.3.1

</td>

</tr>

<tr>

<td style="text-align:left;">

tidyr

</td>

<td style="text-align:left;">

0.8.3

</td>

</tr>

<tr>

<td style="text-align:left;">

tibble

</td>

<td style="text-align:left;">

2.1.1

</td>

</tr>

<tr>

<td style="text-align:left;">

ggplot2

</td>

<td style="text-align:left;">

3.1.1

</td>

</tr>

<tr>

<td style="text-align:left;">

tidyverse

</td>

<td style="text-align:left;">

1.2.1

</td>

</tr>

<tr>

<td style="text-align:left;">

RevoUtils

</td>

<td style="text-align:left;">

11.0.3

</td>

</tr>

<tr>

<td style="text-align:left;">

RevoUtilsMath

</td>

<td style="text-align:left;">

11.0.0

</td>

</tr>

</tbody>

</table>
