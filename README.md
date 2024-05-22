
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11244560.svg)](https://doi.org/10.5281/zenodo.11244560)

# Spatial and temporal analysis on the impact of ultra-low volume indoor insecticide spraying on Aedes aegypti household density

This repository contains the data and code for our paper:

> Anna B. Kawiecki, Amy C. Morrison, Christopher M. Barker (2024).
> *`Spatial and temporal analysis on the impact of ultra-low volume indoor insecticide spraying on Aedes aegypti household density`*.
> Parasites & Vectors <https://doi.org/10.1186/s13071-024-06308-3>

## Contents

The **analysis** directory contains:

[:file_folder: data](/analysis/data): Data used in the analysis.

- [:file_folder: raw_data](/analysis/data/derived_data): The data used
  for this project can only be shared upon request. We provide a sample
  data set with simulated values that allows the code to be tested and
  explored. The provided sample data consists of 20 points randomly
  selected inside the spray zone and 10 points randomly selected in the
  buffer zone for a total of 30 household locations in each study. The
  results do not match those of the paper.

- [:file_folder: derived_data](/analysis/data/derived_data): Analysis of
  the study data. Because there were 2 separate studies in 2013 and
  2014, the are separate analysis for each year.

  - [:file_folder: variables](/analysis/data/derived_data/variables):
    Preparation of variables that are included in the models

  - [:file_folder: model](/analysis/data/derived_data/model): Model
    execution

  - [:file_folder: results](/analysis/data/derived_data/results): Data
    analysis of model outputs and data preparation for plots and
    supplementary material.

[:file_folder: figures](/analysis/figures): Plots

[:file_folder:supplementary-materials](/analysis/supplementary-materials):
Supplementary materials included in the paper.

## How to use the code

The original data and data outputs are not included in this repository,
but are available upon request to the corresponding author according to
the program funder data access and sharing policy. A sample data set of
randomized points in Iquitos, Peru, is provided to facilitate code
comprehension. Below are the contents of the repository and in what
order they should be used.

### Data needed for analysis

**1. Locate sample data sets**

Under [:file_folder: raw_data](/analysis/data/raw_data)

[data.sample.rds](/analysis/data/raw_data/data.sample.rds) : randomized
sample of study households and the dates of spray events and adult
surveys

Variables in the data set:

- “location_code”: de-identified code for each household,

- “date”: adult survey date

- “cycle\_\[1-6\]”: date when 1st through 6th spray events occurred in a
  household, if they occurred

- “exp”: observation in study year 2013 or 2014

- “e_circuit”: the circuit of entomological collection (adult survey)

- “moment”: the moment in relation to the spray events that the adult
  survey occurred ( BASELINE, PRE, FUM, POST\[1-3\])

- “org”: spray performed by study staff or MoH (only in L-2014)

- “zone”: household located in the spray zone or buffer zone

- “AA_TOTAL”: total number of *Ae.aegypti* collected in the household

- “geometry”: geolocation of the household point

[gis.buff.sample.13.rds](/analysis/data/raw_data/gis.buff.sample.13.rds),
[gis.buff.sample.14.rds](/analysis/data/raw_data/gis.buff.sample.14.rds)
Randomized households selected within 1000m of the study households.

[plot.gis.buff.sample.13.jpg](/analysis/data/raw_data/gis.buff.sample.13.rds),
[plot.gis.buff.sample.14.jpg](/analysis/data/raw_data/gis.buff.sample.14.rds):
Visual representations of the sample points in the study area and sample
points within a 1000m buffer of every point in the study area

**2. Create functions that extract INLA model summaries**

[model.fx.Rmd](/analysis/data/derived_data/model/model.fx.Rmd)

*Outputs*:

- fx.fix.eff: outputs the estimated fixed effects for each model

- fx.waic: outputs the WAIC and dWAIC for each model, where the dWAIC is
  the difference in WAIC between the WAIC of the current model and the
  WAIC of the model with the lowest WAIC

- fx.waic.b.g20\]: outputs the WAIC and dWAIC for each model, where the
  dWAIC is the difference in WAIC between the WAIC of the each model and
  the WAIC of the model selected in the within-household analysis
  ($m_{best}$).

### Within-household spray effects (*a*)

**1. Calculate of variables that measure spray effects**

[vbles.time.13.R](/analysis/data/derived_data/variables/variables_2013/vbles.time.13.R),
[vbles.time.14.R](/analysis/data/derived_data/variables/variables_2014/vbles.time.14.R)

*Inputs*:“data_raw.rds”

*Outputs*:“vbles.time.13.rds”,“vbles.time.14.rds”

**2. Run random effect models**

[ulv_model_random_13.R](/analysis/data/derived_data/model/model_2013/ulv_model_random_13.R),
[ulv_model_random_14.R](/analysis/data/derived_data/model/model_2014/ulv_model_random_14.R)

*Inputs*: “vbles.time.13.rds”, “vbles.time.14.rds”

*Outputs*: “m.random.13.rds”, “m.random.14.rds”

**3. Compare random effect models**

In order to decide the random effects to include in the models moving
forward

[model.output.Rmd](/analysis/data/derived_data/model/model.output.Rmd)

*Inputs*:

- Model lists: “m.random.13.rds”, “m.random.14.rds”

- Functions: [fx.waic](/analysis/data/derived_data/model/fx.waic.rds)

*Outputs*: choose to move forward using all 3 random effects.

**3. Run candidate models for within-household spray effect (*a*)**

Recommended to run on multiple cores to reduce computation time.

[ulv_model_time_13.R](/analysis/data/derived_data/model/model_2013/ulv_model_time_13.R),
[ulv_model_time_14.R](/analysis/data/derived_data/model/model_2014/ulv_model_time_14.R)

*Inputs*: “vbles.time.13.rds”, “vbles.time.14.rds”,

*Outputs*: “m.time.13.rds”, “m.time.14.rds”

**4. Compare candidate models measuring within-household spray effects
with those that only include random effects**

Extract WAIC and use these to compare candidate models

[model.output.Rmd](/analysis/data/derived_data/model/model.output.Rmd)

*Inputs*:

- Model lists: “m.random.13.rds”, “m.random.14.rds”, “m.time.13.rds”,
  “m.time.14.rds”

- Functions:

[fx.waic](/analysis/data/derived_data/model/fx.waic.rds),
[fx.fix.eff](/analysis/data/derived_data/model/fx.fix.eff.rds)

*Outputs*: “waic.time.13.rds”, “f.e.time.13.rds”, “waic.time.14.rds”,
“f.e.time.14.rds”

Select the best fitting model or $m_{best}$

### Effects of sprays in neighboring households (*b*)

**1. Calculate of variables that measure spray effects**

Recommended to run on multiple cores to reduce computation time

- <u>1.1 For *b*: decay rates of the spray effect over space:</u>

[vbles.space.decay.13.R](/analysis/data/derived_data/variables/variables_2013/vbles.space.decay.13.R),
[vbles.space.decay.14.R](/analysis/data/derived_data/variables/variables_2013/vbles.space.decay.14.R)

*Inputs*: “vbles.time.13.rds”, “vbles.time.14.rds”

*Outputs*: “vbles.space.decay.13.rds”,“vbles.space.decay.14.rds”

- <u>1.2 For *b*: the proportion houses sprayed within a given distance
  ring:</u>

[vbles.space.ring.13.R](/analysis/data/derived_data/variables/variables_2013/vbles.space.ring.13.R),
[vbles.space.ring.14.R](/analysis/data/derived_data/variables/variables_2014/vbles.space.ring.14.R)

*Inputs*: “vbles.time.13.rds”,
“gis.buff.sample.13.rds”,“vbles.time.14.rds”, “gis.buff.sample.14.rds”

*Outputs*:

- Ring distances based on biological plausibility:
  “vbles.space.ring.bio.13.rds”,“vbles.space.ring.bio.14.rds”

- Ring distances based on even spacing:
  “vbles.space.ring.100m.13.rds”,“vbles.space.ring.100m.14.rds”,

**2. Run candidate models for the effect of sprays in neighboring
households (*b*)**

Recommended to run on multiple cores to reduce computation time.

- <u>Models where *b* variables represent different possible decay rates
  of the spray effect over space:</u>

[ulv_model_space_decay_13.R](/analysis/data/derived_data/model/model_2013/ulv_model_space_decay_13.R),
[ulv_model_space_decay_14.R](/analysis/data/derived_data/model/model_2014/ulv_model_space_decay_14.R)

*Inputs*: “vbles.space.decay.13.rds”, “vbles.space.decay.14.rds”

*Outputs*: model lists: “m.space.decay.13.rds”, “m.space.decay.14.rds”

- <u>Models where *b* variables represent the proportion houses sprayed
  within a given distance ring:</u>

[ulv_model_space_ring_13.R](/analysis/data/derived_data/model/model_2013/ulv_model_space_ring_13.R),
[ulv_model_space_ring_14.R](/analysis/data/derived_data/model/model_2014/ulv_model_space_ring_14.R)

*Inputs*: “vbles.space.ring.bio.13.rds”,“vbles.space.ring.100m.13.rds”,
“vbles.space.ring.bio.14.rds”, “vbles.space.ring.100m.14.rds”,

*Outputs*: model lists: “m.space.rings.13.rds”, “m.space.rings.14.rds”

**3. Compare models measuring effects of sprays in neighboring
households to $m_{best}$**

[model.output.Rmd](/analysis/data/derived_data/model/model.output.Rmd)

*Inputs*:

- Model lists: “m.time.13.rds”, “m.time.14.rds”,“m.space.decay.13.rds”,
  “m.space.decay.14.rds”,“m.space.rings.13.rds”, “m.space.rings.14.rds”

- Function:
  [fx.waic.b.g20](/analysis/data/derived_data/model/fx.waic.b.g20.rds)

*Outputs*:

- Model lists: “m.space.13.rds”, “m.space.14.rds”

- Model summaries: “waic.space.13.rds”, “waic.space.14.rds”

**4. Perform prior sensitivity analysis**

[model.prior.sensitivity.analysis.13.Rmd](/analysis/data/derived_data/model/model_2013/model.prior.sensitivity.analysis.13.Rmd),
[model.prior.sensitivity.analysis.14.Rmd](/analysis/data/derived_data/model/model_2014/model.prior.sensitivity.analysis.14.Rmd)

*Inputs*: “vbles.space.decay.13.rds”, “vbles.space.decay.14.rds”

*Outputs*:  
- Model lists: “m.priors.13.rds”, “m.priors.14.rds”

- Model summaries: “waic.priors.13.rds”, “fe.priors.13.rds”,
  “waic.priors.14.rds”, “fe.priors.14.rds”

### Develop result data sets

Under [:file_folder: results](/analysis/data/derived_data/results/)

**1. Wrangle model summary data for figures and supplementary
materials**

[results.data.Rmd](/analysis/data/derived_data/results/results.data.Rmd)

*Inputs*: “data.sample.rds”, “waic.time.13.rds”,“f.e.time.13.rds”,
“waic.time.14.rds”,“f.e.time.14.rds”, “waic.space.13.rds”,
“waic.space.14.rds”, “vbles.time.14.rds”, “gis.buff.sample.14.rds”,
“waic.priors.13.rds”, “fe.priors.13.rds”, “waic.priors.14.rds”,
“fe.priors.14.rds”

*Outputs*:

- Decay curve data frames to visualize the decay weighting schemes for
  time and space: “decay.time.rds”, “decay.distance.rds”

- Model output (WAIC and fixed effect estimates) for both 2013 and 2014
  formatted for figure development : “waic.time.rds”, “f.e.time.rds”,
  “waic.space.rds”,

- Data frame that combines model WAIC, dWAIC and the decay curves that
  are used to weight time-since-spray: “curves.model.time.rds”

- Data frame with descriptions of all the candidate variables used in
  this analysis:“effects.rds”

- Example distance rings around a sample house: “i.bio.rds”,
  “i.100m.rds”, “i.buffer.bio.rds”, “i.buffer.100m.rds”

- Data frame that allow for estimation of the % reduction for *Ae.
  Aegypti* abundance following sprays: “mg20.14.red.rds”

- Model WAIC and fixed effects for prior sensitivity analysis:
  “fe.priors.rds”, “waic.priors.rds”

### Develop figures from the model results

Under [:file_folder: figures](/analysis/figures)

[ulv_figures.Rmd](/analysis/figures/ulv_figures.Rmd)

*Inputs*: “gis.buff.sample.14.rds”,“i.bio.rds”, “i.100m.rds”,
“i.buffer.bio.rds”, “i.buffer.100m.rds”
“curves.model.time.rds”,“waic.time.rds”,“mg20.14.red.rds”

*Outputs*: “\*.jpg”

### Develop additional file

Under
[:file_folder:supplementary-materials](/analysis/supplementary-materials)

[additional.file.1.Rmd](/analysis/supplementary-materials/additional.file.1.Rmd)

*Inputs*: “decay.time.rds”, “decay.distance.rds”,
“curves.model.time.rds” , “waic.time.rds”, “f.e.time.rds”,
“waic.space.rds”, “fe.priors.rds”, “waic.priors.rds”

*Outputs*:
[Additional.file.1.pdf](/analysis/supplementary-materials/Additional.file.1.pdf)

## How to cite

Please cite this compendium as:

> Kawiecki, (2024). *Compendium of R code and data for
> `Spatial and temporal analysis on the impact of ultra-low volume indoor insecticide spraying on Aedes aegypti household density`*.
> Accessed 22 May 2024. Online at
> <https://doi.org/10.5281/zenodo.11244560>

## How to run in your browser or download and run locally

This research compendium has been developed using the statistical
programming language R. To work with the compendium, you will need
installed on your computer the [R
software](https://cloud.r-project.org/) itself and optionally [RStudio
Desktop](https://rstudio.com/products/rstudio/download/).

### Licenses

**Text and figures :**
[CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

**Code :** See the [DESCRIPTION](DESCRIPTION) file

**Data :** De-identified data sets are available upon request to the
corresponding author.
