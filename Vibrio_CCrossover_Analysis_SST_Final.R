################################################################################
#MSc Vibrio Study
#Code based on 1301 ENVIRONMENTAL EPIDEMIOLOGY - CHERIE PART & MAX EYRE
################################################################################

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#1 Set up and package load
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#set working directory
setwd("~/Documents/LSHTM/Thesis/Environmental data")

# load packages.
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

library(tidyverse)

pacman::p_load(dplyr, lubridate, gnm, Epi, survival, dlnm, ggplot2, gtsummary, gt, flextable)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#2 Prepare and clean the case data
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#Import and look at case data
cases <- read.csv("Outcome Data Vibrio.csv")
head(cases)                             # show first 6 rows
str(cases)                              # structure of the data

#Convert the date from character to date class
cases$date <- as.Date(cases$date, format = "%d/%m/%Y")

#Make sure data is ordered by date chronologically from 2010 to 2025
cases <- cases %>% arrange(date)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#3 Descriptive Analysis
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

### A. How many non-cholera vibrio infections were diagnosed and over what time period?
sum(cases$isb)                          # total N cases of vibrio cases
range(cases$date)                       # date range

### B. Which age group had the highest number of vibrio cases?
cases %>%
  group_by(mage) %>%
  summarise(n=sum(isb),
            prop = n/sum(cases$isb)*100)


### C. What percentage of vibrio cases occurred in women and men
cases %>%
  group_by(sex) %>%
  summarise(n=sum(isb),
            prop = n/sum(cases$isb)*100)

### C. What percentage of vibrio cases occurred in what area
cases %>%
  group_by(area) %>%
  summarise(n=sum(isb),
            prop = n/sum(cases$isb)*100)

### C. What percentage of vibrio cases occurred with what exposure type
cases %>%
  group_by(exposure_type) %>%
  summarise(n=sum(isb),
            prop = n/sum(cases$isb)*100)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#4 Epidemiological tables 1 & 2 (pooled and stratified by area)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

############ Table 1: Make pooled epi table word doc ##############
epi_table <- cases %>%
  filter(id <= 249) %>%
  mutate(across(where(is.factor), droplevels)) %>%
  select(sex, mage, exposure_type, vibrio_species, area) %>%
  tbl_summary(
    label = list(sex ~ "Sex", mage ~ "Age", exposure_type ~ "Exposure type", vibrio_species ~ "Vibrio species", area ~ "Region"),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = all_categorical() ~ 1,
    missing = "no"
  ) %>%
  modify_header(label = "**Demographic data**", stat_0 = "**Cases (n)**")

epi_table
epi_table %>% as_flex_table() %>% save_as_docx(path = "Epi_Table.docx")

############ Table 2: Epi_Table_Stratified_by_Area word doc ##############

processed_cases <- cases %>%
  filter(!is.na(area) & area != "" & area != " ") %>%
  mutate(
    vibrio_species = as.factor(vibrio_species),
    sex = as.factor(sex),
    area = as.factor(area),
    exposure_type = as.factor(exposure_type),
    mage = forcats::fct_relevel(mage, c("<20", "20-34", "35-49", "50-64", "65-79", "80+"))
  )

area_stratified_table <- processed_cases %>%
  select(sex, mage, exposure_type, vibrio_species, area) %>%
  tbl_summary(
    by = area,
    type = list(mage ~ "categorical"),
    label = list(sex ~ "Sex", mage ~ "Age", exposure_type ~ "Exposure Type", vibrio_species ~ "Vibrio species"),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = all_categorical() ~ 1,
    missing = "no"
  ) %>%
  modify_header(label = "**Demographic data**")

area_stratified_table
area_stratified_table %>% as_flex_table() %>% save_as_docx(path = "Epi_Table_by_Area.docx")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 5. Exposure Data import and cleaning
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#Open and read exposure data
exposure <- read.csv("exposure.csv")    # import data
head(exposure)                          # look at first 6 rows
str(exposure)                           # check structure
exposure$date <- as.Date(exposure$date, "%d/%m/%Y")    # make dates consistent format
str(exposure)                           # check

# make sure data are ordered by date
exposure <- exposure %>% arrange(date)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 6. Time stratfied logistic regression. #Pooled. Case crossover design
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# 1. Load libraries
library(tidyverse)
library(lubridate)
library(survival)  # For clogit - logistic regression
library(dlnm)      # For onebasis and crosspred

n_distinct(cases$id) == nrow(cases) # check all case ids are unique

# 1. Standardise date formats in both exposure and cases data for case crossover
cases$date            <- as.Date(parse_date_time(cases$date, orders = c("ymd", "dmy", "mdy")))
exposure_cleaned$date <- as.Date(parse_date_time(exposure_cleaned$date, orders = c("ymd", "dmy", "mdy")))

# 2. Ignore spacing in cells, and join text characters for 'area'
cases$area            <- trimws(tolower(as.character(cases$area)))
exposure_cleaned$area <- trimws(tolower(as.character(exposure_cleaned$area)))

# Sort files chronologically to maintain consistency
cases            <- cases %>% arrange(date)
exposure_cleaned <- exposure_cleaned %>% arrange(date)

# include 249 cases, add dow, month, and year for matching
cases_cleaned <- cases %>%
  filter(id <= 249) %>%
  mutate(case_date = date,
         dow   = wday(date, week_start = 1),
         month = month(date),
         year  = year(date))

# create a dataframe of all days in the study period
start_date <- min(cases_cleaned$date, na.rm = TRUE)
end_date   <- max(cases_cleaned$date, na.rm = TRUE)

alldays <- data.frame(date = seq(start_date, end_date, by = "day")) %>%
  mutate(dow   = wday(date, week_start = 1),
         month = month(date),
         year  = year(date))

# add control days: attach all calendar dates that fall on the same dow,
# month, and year as the case_date.
cc_data <- cases_cleaned %>%
  select(id, case_date, dow, month, year, mage, area, sex) %>%
  inner_join(alldays, by = c("dow","month","year"), relationship = "many-to-many") %>%
  mutate(isb = as.integer(date == case_date))

# Prevent duplicated rows
cc_data <- cc_data %>% distinct(id, date, .keep_all = TRUE)

# Join exposure data by date and area
cc_data <- cc_data %>%
  left_join(exposure_cleaned, by = c("date" = "date", "area" = "area"))


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 7 calculate SST specific lags
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

cc_data <- cc_data %>%
  arrange(id, date) %>%
  group_by(id) %>%
  mutate(
    hi18_lag0 = pmax(SSTmean - 18, 0),
    hi18_lag1 = pmax(dplyr::lag(SSTmean, 1) - 18, 0),
    hi18_lag2 = pmax(dplyr::lag(SSTmean, 2) - 18, 0),
    hi18_lag3 = pmax(dplyr::lag(SSTmean, 3) - 18, 0)
  ) %>%
  ungroup()

# show no missing rows in the data to confirm matching
print(paste("Total Observations Matched:", nrow(cc_data)))
print(paste("Remaining SSTmean missing items:", sum(is.na(cc_data$SSTmean))))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 9. PART A: Pooled analysis. Conditional logistic regression. Case crossover design
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Create a natural cubic spline basis for sea surface temperature (lag 0)
ncstemp <- onebasis(cc_data$SSTmean, "ns", df = 3)

# Fit the conditional logistic regression model
modclog <- clogit(isb ~ ncstemp + strata(interaction(id, area)), data = cc_data)
summary(modclog) # The non-linear conditional logistic regression model confirms a significant overall association

# Predict the exposure-response curve and plot the spline
predclog <- crosspred(ncstemp, modclog, by = 1)
plot(predclog, "overall", xlab = "Daily mean temperature (°C), lag 0", ylab = "Odds Ratio", lwd = 2, ci.arg = list(density = 100, col = grey(0.5)))

# set 18°C threshold baseline for SST, quantify the heat effect and stat significance
modlthres <- clogit(isb ~ hi18_lag0 + strata(interaction(id, area)), data = cc_data)
summary(modlthres)
# for every 1°C increase in daily mean sea surface temperature above the 18°C threshold, the odds of a diagnosed Vibrio infection increase significantly by 39.9%
#(OR = 1.399; 95% CI: [1.195, 1.636]; (p < 0.001)).

#run logistic regression for SSTmean Pmean, SSSmean and Rmean to see if the lagged exposure variables
#are associated with outcome.

# 1. Run all the separate models
mod_SST   <- clogit(isb ~ hi18_lag0 + strata(interaction(id, area)), data = cc_data)
mod_P     <- clogit(isb ~ Pmean     + strata(interaction(id, area)), data = cc_data)
mod_SSS   <- clogit(isb ~ SSSmean   + strata(interaction(id, area)), data = cc_data)
mod_Rmean <- clogit(isb ~ Rmean     + strata(interaction(id, area)), data = cc_data)

# 2. Show standard row coefficients from each model
row_SST   <- summary(mod_SST)$coefficients["hi18_lag0", ]
row_P     <- summary(mod_P)$coefficients["Pmean", ]
row_SSS   <- summary(mod_SSS)$coefficients["SSSmean", ]
row_Rmean <- summary(mod_Rmean)$coefficients["Rmean", ]

# 3. Combine exposure variables and associations into one table
results_table <- rbind(
  SSTmean_Above18 = row_SST,
  Pmean           = row_P,
  SSSmean         = row_SSS,
  Rmean           = row_Rmean
)

# 4. Clean the columns and show table
results_table <- results_table[, c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")]
round(results_table, 4)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 10. #Unadjusted distributed lag sea surface temperature model
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

exposure_cleaned <- exposure_cleaned %>%
  arrange(area, date) %>%
  group_by(area) %>%
  mutate(
    hi18_lag0 = pmax(SSTmean - 18, 0),
    hi18_lag1 = pmax(dplyr::lag(SSTmean, 1) - 18, 0),
    hi18_lag2 = pmax(dplyr::lag(SSTmean, 2) - 18, 0),
    hi18_lag3 = pmax(dplyr::lag(SSTmean, 3) - 18, 0)
  ) %>%
  ungroup()

exposure_cleaned <- exposure_cleaned %>%
  mutate(across(c(hi18_lag0, hi18_lag1, hi18_lag2, hi18_lag3), ~replace_na(.x, 0)))

cc_data <- cases_cleaned %>%
  select(id, case_date, dow, month, year, mage, area, sex) %>%
  inner_join(alldays, by = c("dow","month","year"), relationship = "many-to-many") %>%
  mutate(isb = as.integer(date == case_date)) %>%
  distinct(id, date, .keep_all = TRUE) %>%
  left_join(exposure_cleaned, by = c("date" = "date", "area" = "area"))

modlt03_raw <- clogit(isb ~ hi18_lag0 + hi18_lag1 + hi18_lag2 + hi18_lag3 +
                        strata(interaction(id, area)), data = cc_data)
summary(modlt03_raw)

# 3. Extract and pool the raw coefficients
lags    <- c("hi18_lag0","hi18_lag1","hi18_lag2","hi18_lag3")
b_raw   <- coef(modlt03_raw)[lags]
V_raw   <- vcov(modlt03_raw)[lags, lags]
L       <- rep(1, length(lags))
est_raw <- sum(L * b_raw)
se_raw  <- sqrt(as.numeric(t(L) %*% V_raw %*% L))
OR_raw  <- exp(est_raw)
LCL_raw <- exp(est_raw - 1.96 * se_raw)
UCL_raw <- exp(est_raw + 1.96 * se_raw)

# show cumulative unadjusted OR
cat("\n--- POOLED CUMULATIVE UNADJUSTED OVERALL RE-CALCULATED EFFECT ---\n")
print(c(Unadjusted_OR = OR_raw, Lower_CI = LCL_raw, Upper_CI = UCL_raw))


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#11. Intermediate single variable (rainfall) adjusted model
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
modlt03_Pmean <- clogit(isb ~ hi18_lag0 + hi18_lag1 + hi18_lag2 + hi18_lag3 + Pmean + strata(interaction(id, area)), data = cc_data)
summary(modlt03_Pmean)

# Extract and pool the rainfall-adjusted coefficients
b_Pmean   <- coef(modlt03_Pmean)[lags]
V_Pmean   <- vcov(modlt03_Pmean)[lags, lags]
est_Pmean <- sum(L * b_Pmean)
se_Pmean  <- sqrt(as.numeric(t(L) %*% V_Pmean %*% L))
OR_Pmean  <- exp(est_Pmean)
LCL_Pmean <- exp(est_Pmean - 1.96 * se_Pmean)
UCL_Pmean <- exp(est_Pmean + 1.96 * se_Pmean)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 12. Fully adjusted multivariate model
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

modlt03 <- clogit(isb ~ hi18_lag0 + hi18_lag1 + hi18_lag2 + hi18_lag3 + Pmean + SSSmean + Rmean + strata(interaction(id, area)), data = cc_data)
summary(modlt03)

# Extract and pool fully adjusted coefficients
b_adj   <- coef(modlt03)[lags]
V_adj   <- vcov(modlt03)[lags, lags]
est_adj <- sum(L * b_adj)
se_adj  <- sqrt(as.numeric(t(L) %*% V_adj %*% L))
OR_adj  <- exp(est_adj)
LCL_adj <- exp(est_adj - 1.96 * se_adj)
UCL_adj <- exp(est_adj + 1.96 * se_adj)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 13. Model tests
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
anova(modlt03_Pmean, modlt03, test = "Chisq")

cat("\n--- CUMULATIVE ODDS RATIO COMPARISONS ---\n")
print(c(Unadjusted_OR = OR_raw, Lower_CI = LCL_raw, Upper_CI = UCL_raw))
print(c(Rain_Adjusted_OR = OR_Pmean, Lower_CI = LCL_Pmean, Upper_CI = UCL_Pmean))
print(c(Fully_Adjusted_OR = OR_adj, Lower_CI = LCL_adj, Upper_CI = UCL_adj))

#after adjusting for potential confounders the relationship between isb and
#SSTmean is unaffected

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 14. PART B: Area stratified logistic regression SST, adjusting for other variables
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

regional_results_full <- cc_data %>%
  group_by(area) %>%
  group_modify(~{

    # Clean the subset: Keep only patient IDs that have BOTH case and control days
    clean_sub <- .x %>%
      group_by(id) %>%
      filter(any(isb == 1) & any(isb == 0)) %>%
      ungroup()

    tryCatch({
      mod <- clogit(isb ~ SSTmean + Pmean + SSSmean + Rmean + strata(id), data = clean_sub)

      # shows all exposure rows for this region
      as.data.frame(summary(mod)$coefficients) %>%
        tibble::rownames_to_column("term")

    }, error = function(e) {
      data.frame(term = "Model Convergence Failure", coef = NA, `exp(coef)` = NA, `se(coef)` = NA, z = NA, `Pr(>|z|)` = NA, check.names = FALSE)
    })
  }) %>%
  ungroup()

# View
print(regional_results_full, n = 30)

