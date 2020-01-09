
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Scripts

Below are descriptions of each script (organized under the subfolders
within this folder). This is the intended order of executing these
scripts.

## PhenotypeData

1.  `phenotype_data_adjustment.R`: Calculate adjusted genotype means
    within each experimental trial according to the specific
    experimental design of that trial.
2.  `phenotype_data_summary.R`: Calculates summary statistics and fits
    linear mixed models to analyze the phenotypic data. Also included
    scripts for producing geographic maps of the trial locations.

## EnvironmentalVariables

1.  `collect_environmental_variables.R`: Queries NOAA and SoilWeb
    databases to gather daily weather variables (temperature,
    precipitation) and location-specific soil variables.
2.  `predict_growth_stage.R`: Uses temperature data gathered from NOAA
    to predict three barley growth stages (vegetative, flowering, grain
    filling) based on growing degree days and a simple barley growth
    model.
3.  `manipulate_analyze_environmental_variables.R`: Calculates summaries
    of environmental variables at different growth stages or across the
    entire growing season (i.e. isothermality, temperature range, etc.)
4.  `model_environmental_variables.R`: Determines the environmental
    variables that are associated with the environmental mean or
    environmental interaction PCA scores.

## Clustering

1.  `environmental_clustering.R`: Uses model-based clustering to group
    environments together. Analyzes data by calculating heritability
    with or without clustering. Ranks environments based on distance
    measures.
2.  `cluster_and_rank_analysis.R`: Analyze the similarity measures,
    environmental ranks, and clusters.

## Predictions

1.  `distance_rank_predictions.R`: Uses distances measures to add
    environments to a training set to make predictions
2.  `all_data_cluster_predictions.R`: Uses assigned environmental
    clusters to make predictions.

Note that scripts 1 and 2 are written to run on the Minnesota
Supercomputing Insititute.

3.  `prediction_analysis.R`: Analyze the output from the prediction
    scripts.

## Other scripts

1.  `figures.R`: Code to produce various other figures.
