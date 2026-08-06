################################################################################
#MSc Vibrio Study
# Code based on 1301 ENVIRONMENTAL EPIDEMIOLOGY
# Made by MAX EYRE - MARCH 2026
################################################################################

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#1 Set up and package load
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#set working directory
setwd("~/Documents/LSHTM/Thesis/Environmental data")

# load packages. dplyr for cleaning data; gnm & survival for regression;
#dlnm for distributed lag non-linear models
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, lubridate, gnm, Epi, survival, dlnm, ggplot2)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#2 Prepare case data
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# import and inspect case data (Outcome Data Vibrio.csv file)
cases <- read.csv("Outcome Data Vibrio.csv") #import data
head(cases)                             # first 6 rows
str(cases)                              # structure of the data, variables

# convert 'date' from character to proper date format
cases$date <- as.Date(cases$date, format = "%d/%m/%Y")
# convert remaining characters (sex or area, etc) to factors for comparison
cases <- cases %>% mutate(across(where(is.character), as.factor))
#   and check format
str(cases)

# order data by date chronologically
cases <- cases %>% arrange(date)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#3 Basic descriptive Analysis
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

### A. How many non-cholera vibrio infections were recorded and over what time period?
sum(cases$isb)                          # total N cases of vibrio
range(cases$date)                       # date range

### B. Which age category had the highest number of vibrio cases?
cases %>%
  group_by(mage) %>%
  summarise(n=sum(isb),                 # sum number
            prop = n/nrow(cases)*100)   # calculate proportion (%)

#Bar chart all areas

#Bar chart, by area


### C. What percentage of vibrio cases occurred in women and men
cases %>%
  group_by(sex) %>%
  summarise(n=sum(isb),
            prop = n/nrow(cases)*100)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#4 Creating plots and exporting epidemiological summary tables
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Summary of a single variable, sex
summary(cases$sex)

#Cross table, cases by area and sex
table(cases$area, cases$sex)

# Group by many variables - summarises demographic data
cases %>%
  group_by(area, sex, mage,exposure_type, vibrio_species) %>%
  summarise(cases_n = n(), .groups = "drop")

# Install the package if not already: install.packages("gtsummary")
library(gtsummary)
library(dplyr)

#A make Epi table
epi_table <- cases %>%
  select(sex, mage, exposure_type, vibrio_species, area) %>%
  tbl_summary(
    label = list(
      sex ~ "Sex",
      mage ~ "Age",
      exposure_type ~ "Exposure type",
      vibrio_species ~ "Vibrio species",
      area ~ "Region"
    ),
    statistic = list(
      all_continuous() ~ "{mean} (y)\n{sd} (y)\n{min}–{max}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no"
  ) %>%
  # FIX: Use '=' instead of '~' inside modify_header
  modify_header(
    label = "**Demographic data**",
    stat_0 = "**Cases (n)**"
  )

#view epitable 1
epi_table

# Add this to the end of your table pipeline:
write.csv(as.data.frame(epi_table), "Epi_Table.csv", row.names = FALSE)

#word doc save as
gtsummary::as_gt(epi_table) %>% gt::gtsave("Epi_Table.docx")

#--------

#B Epi_Table_Stratified_by_Area

library(dplyr)
library(gtsummary)
library(gt)

#. Process data, ensure chronological age order
processed_cases <- cases %>%
  mutate(
    mage = forcats::fct_relevel(mage, c("<20", "20 - 39", "40 - 59", "60 - 79", "80+"))
  )


#. Build the stratified table by area

area_stratified_table <- processed_cases %>%
  select(sex, mage, exposure_type, vibrio_species, area) %>%
  tbl_summary(
    by = area, # Group side-by-side columns by area
    type = list(mage ~ "categorical"),
    label = list(
      sex ~ "Sex",
      mage ~ "Age",
      exposure_type ~ "Exposure Type",
      vibrio_species ~ "Vibrio species"
    ),
    statistic = list(
      all_categorical() ~ "{n} ({p}%)" #counts and percentage for categories
    ),
    missing = "no"
  ) %>%
  modify_header(
    label = "**Demographic data**"
  )


#View are stratified table
area_stratified_table

#. Export directly to your Word Document path
gtsummary::as_gt(area_stratified_table) %>% gt::gtsave("Epi_Table_by_Area.docx")
