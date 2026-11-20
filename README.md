# Environmental TCE Exposure, Social Vulnerability, and Cognitive Difficulty in Pennsylvania

**BMIN5030 Final Project** · Michael Condiff

## Overview

This project asks whether environmental exposure to trichloroethylene (TCE), an industrial solvent and common groundwater contaminant, is associated with cognitive difficulty across Pennsylvania census tracts after accounting for social vulnerability.

Using only publicly available data, I built three tract-level TCE exposure measures, each reflecting a different exposure pathway:

- **Water:** TCE detections at monitoring sites from the EPA/USGS Water Quality Portal (2000-2020)
- **Superfund:** distance from each tract's population-weighted center to the nearest Superfund site with TCE as a contaminant of concern (EPA SEMS, via the Envirofacts API)
- **Air:** modeled TCE exposure concentrations from EPA's 2020 AirToxScreen

The outcome is the percentage of adults aged 18-64 reporting a cognitive difficulty in the American Community Survey (2018-2022). Analyses include descriptive comparisons of exposed and unexposed tracts, nested linear regression adjusting for social vulnerability and urbanicity, county fixed effects, and spatial error models.

**Main finding:** none of the three exposure measures was associated with higher cognitive difficulty after adjustment. Crude associations were explained by social vulnerability and by differences between counties. Exposure itself was unevenly distributed: tracts with TCE detected in water were more urban, more racially and ethnically diverse, and had more poverty. A planned next step is a phenome-wide association study (PheWAS) linking these exposure measures to Penn Medicine BioBank participants.

## Repository contents

| File or folder | Description |
|---|---|
| `final_project_draft6.qmd` | The full report: data retrieval, cleaning, exposure construction, analysis, and write-up. All code is in this file. |
| `final_project_draft6.html` | Rendered report |
| `final_project.R` | Earlier standalone script for the Superfund and AirToxScreen exposures. Its code is now in the `.qmd`; kept for reference only. |
| `data/raw/` | Downloaded source data and cached query results (created on first render; see below) |
| `data/processed/` | Tract-level outputs: exposure measures, the analytic dataset, and the AirToxScreen facility table |
| `output/figures/` | Saved figures |
| `README.md` | This file |

## Data sources

| Source | What it provides |
|---|---|
| [Water Quality Portal](https://www.waterqualitydata.us/) | TCE sample results and monitoring site locations, Pennsylvania, 2000-2020 |
| [EPA Envirofacts SEMS tables](https://www.epa.gov/enviro/envirofacts-data-service-api) | Superfund sites and their contaminants of concern |
| [EPA AirToxScreen 2020](https://www.epa.gov/AirToxScreen/2020-airtoxscreen-assessment-results) | Block-level modeled TCE exposure concentrations (EPA Region 3) and the facility emissions inventory |
| [U.S. Census Bureau ACS 2018-2022](https://www.census.gov/programs-surveys/acs) | Cognitive difficulty (table B18104) and social vulnerability measures, via `tidycensus` |
| [2020 Decennial Census](https://www.census.gov/programs-surveys/decennial-census/decade/2020/2020-census-main.html) | Block populations (for weighting) and tract centers of population |
| [TIGER/Line shapefiles](https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html) | 2020 tract and state boundaries, via `tigris` |

No individual-level or protected health data are included in this repository.

## Reproducing the analysis

1. **Install R packages** (once, in the R console):

   ```r
   install.packages(c("tidyverse", "sf", "tigris", "tidycensus", "dataRetrieval",
                      "ggspatial", "readxl", "janitor", "gtsummary", "sandwich",
                      "lmtest", "broom", "spdep", "spatialreg", "knitr"),
                    type = "binary")
   ```

2. **Set up a Census API key** (free, from <https://api.census.gov/data/key_signup.html>):

   ```r
   tidycensus::census_api_key("YOUR_KEY", install = TRUE)
   ```

   Restart R afterward. The key is stored in your `.Renviron` file and should never be committed to the repository.

3. **Open the RStudio project** and render `final_project_draft5.qmd` (Ctrl+Shift+K).

The first render downloads all source data and may take several minutes. Every download and query result is saved under `data/raw/` and `data/processed/`, and later renders read the saved copies instead of downloading again. To refresh a data source, delete its saved file and render again.

## Notes on data files and Git

Some downloaded files are large (the unzipped AirToxScreen Region 3 file in particular) and may exceed GitHub's file size limit. Because the report re-creates everything in `data/raw/` automatically, that folder can be excluded from version control with a `.gitignore` entry:

```
data/raw/
```

The smaller tract-level files in `data/processed/` can be committed so readers can inspect the results without rendering.

## Acknowledgments

Thanks to Dr. Joe Romano (BMIN), who helped focus the project on TCE and Superfund sites, and Colleen Morse (Penn Medicine BioBank), who advised on linking tract-level exposures to BioBank participants without sharing location data.
