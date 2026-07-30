# rdavis_vibrio

MSc thesis project (LSHTM, Environmental Epidemiology) studying non-cholera *Vibrio* infections in Cádiz, Charente-Maritime, Arcachon, and Stockholm (2010–2025), and their association with sea/air environmental exposures (sea surface temperature, salinity, air temperature, precipitation, radiation) using a time-stratified case-crossover design.

All scripts are R and are written to be run from a working directory containing the raw case data (`Outcome Data Vibrio.csv`), an assembled `exposure.csv`, and the relevant `.nc` (NetCDF) climate/ocean data files. Several scripts have hard-coded `setwd()` calls pointing at the original author's local folders — update these before rerunning.

## Files

### Data preparation — NetCDF to CSV
- **`Nc file to csv.R`** — Loops over multiple years and three NetCDF variables (temperature, precipitation, radiation) for a given site, computes daily spatial averages, converts units (Kelvin→Celsius, metres→mm), and compiles everything into one CSV (e.g. `Cadiz_daily_averages_2010_2025.csv`). Site-specific: working directory and file prefixes need editing per location.
- **`SST_NC_to_CSV.R`** — Reads a sea surface temperature NetCDF file (`SST_daily_Stockholm_2010_2025.nc`), computes the daily spatial mean, exports to CSV, and plots the daily SST timeline.
- **`NC to CSV_SSS_Data.R`** — Same as above but for sea surface salinity (`SSS_France_2010_2025.nc`), exporting `France_daily_average_2010_2025_SSS.csv` and a timeline plot.

### Descriptive / exploratory analysis
- **`Vibrio_analysis_descriptive_for_CC_study.R`** — Early/simple descriptive analysis of the case data: total case counts, date range, case counts by age group and by sex.
- **`epi_tables.R`** — Expanded descriptive analysis plus formatted summary tables (via `gtsummary`/`gt`): a pooled demographic table and a table stratified by study area, exported as CSV/Word documents.

### Epidemic curves
- **`epicurve_vibrio_all_years.R`** — Builds a faceted epicurve (area × year grid) of monthly case counts with seasonal background shading; saves `vibrio_epi_curves_seasonal_grid.png`.
- **`epicurve_vibrio_combined_months.R`** — Same idea but pools all years together into one seasonal profile per area (area stacked rows); saves `vibrio_combined_seasonal_curves.png`.

### Main analysis
- **`Logistic_Regression_code`** — The core analysis script (plain text, no `.R` extension, ~1000 lines). End-to-end pipeline:
  1. Loads and cleans case data; builds descriptive stats and demographic tables (pooled + stratified by area).
  2. Loads exposure data (`exposure.csv`) and builds a **time-stratified case-crossover dataset**: for each case, matches control days with the same day-of-week/month/year within the study period, then joins environmental exposures by date and area.
  3. Categorizes SST/SSS exposures and explores lag structures (0–3 days).
  4. Fits **conditional logistic regression** (`clogit`) models: univariate screening across exposures/lags, an SST×SSS interaction model, area-stratified SST+SSS models, precipitation-only models, and a fully adjusted SST+SSS+Pmean model.
  5. Computes a Spearman correlation matrix between exposure variables (checks multicollinearity).
  6. Produces forest plots of odds ratios and exports all results tables as Word documents (`.docx`) via `flextable`/`officer`.
  7. Generates shaded seasonal timeline plots of cases vs. environmental variables per area.

## Notes
- Several scripts duplicate logic (e.g. descriptive stats and epi tables appear in `Vibrio_analysis_descriptive_for_CC_study.R`, `epi_tables.R`, and again inside `Logistic_Regression_code`) — `Logistic_Regression_code` is the most complete/current version and supersedes the smaller scripts.
- Outputs (CSVs, `.docx` tables, `.png` figures) are written to the working directory and are not stored in this repo.
