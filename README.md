# 731 Final Project

This repository contains code and data for a simulation-based study of small area estimation (SAE) methods, motivated by Demographic and Health Survey (DHS) data from Zambia.

---

## Repository Structure

```
├── source/
│   ├── k0_5_pop.r
│   ├── helper.r
│   ├── clusterModel_new.R
│   ├── directEST_national.R
│   └── fhModel_new.R
│
├── data/
│   ├── Zambia/
│   │   ├── country_shp_analysis.rds
│   │   ├── zmb_frame_ea.rds
│   │   ├── zmb_sample_ea.rds
│   │   ├── zmb_ppp_2010_1km_Aggregated_UNadj.tif
│   │   └── zmb_ppp_2018_1km_Aggregated_UNadj.tif
│   │
│   └── subpop/
│       ├── zmb_f_0_2018_1km.tif
│       ├── zmb_f_1_2018_1km.tif
│       ├── zmb_m_0_2018_1km.tif
│       ├── zmb_m_1_2018_1km.tif
│       └── zmb_k0_5_2018_1km.tif
│
├── analysis/
│   └── sim_results.R
│
├── figures/
│
├── 731_final.Rproj
└── README.md
```


---

## Code Overview

### `source/`
Core functions for estimation methods and data construction:

- `directEST_national.R`: design-based direct estimator
- `fhModel_new.R`: Fay–Herriot area-level model
- `clusterModel_new.R`: unit-level hierarchical model
- `helper.r`: supporting functions
- `k0_5_pop.r`: downloading and creating for the under-5 (k0–5) population from worldPop.

---

### `analysis/`

- `sim_results.R`: main simulation script, runs all models and produces results

---

## Data

### `data/Zambia/`

Contains DHS-based inputs and spatial data:

- `country_shp_analysis.rds`: administrative boundaries  
- `zmb_frame_ea.rds`: sampling frame  
- `zmb_sample_ea.rds`: sampled EA data  
- WorldPop raster data  

---

### `data/subpop`
Contains 1km resolution population rasters:

- female/male population (ages 0–5)
- total population (`zmb_k0_5_2018_1km.tif`)


---

## Reproducibility

To reproduce results:

1. Prepare population data
```
source("source/k0_5_pop.r")
```

2. Run simulation
```
source("analysis/sim_results.R")
```

3. Outputs:
- Results stored in R objects
- Figures saved in `figures/`

---

## Notes

- RStudio project structure
- Relative paths assumed
- Large `.RData` files excluded

