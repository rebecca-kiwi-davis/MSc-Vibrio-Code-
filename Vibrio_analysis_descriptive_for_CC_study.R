################################
##### DESCRIPTIVE ANALYSIS #####
################################

#set working directory
setwd("~/Documents/LSHTM/Thesis/Environmental data")

# load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, lubridate, gnm, Epi, survival, dlnm, ggplot2)

# import and inspect the case data (Outcome Data Vibrio.csv file) 
cases <- read.csv("Outcome Data Vibrio.csv")
head(cases)                             # first 6 rows
str(cases)                              # structure of the data

# convert 'date' from character to date class
cases$date <- as.Date(cases$date, format = "%d/%m/%Y")
# convert remaining characters to factors
cases <- cases %>% mutate(across(where(is.character), as.factor))
#   and check
str(cases)

# order data by date 
cases <- cases %>% arrange(date)

### 1. How many non-cholera vibrio infections were recorded and over what time period?
sum(cases$isb)                          # total N cases of vibrio
range(cases$date)                       # date range

### 2. Which age category had the highest number of vibrio cases?
cases %>% 
  group_by(mage) %>%
  summarise(n=sum(isb),                 # sum number
            prop = n/nrow(cases)*100)   # calculate proportion (%)

### 3. What percentage of vibrio cases occurred in women and men
cases %>% 
  group_by(sex) %>%
  summarise(n=sum(isb), 
            prop = n/nrow(cases)*100)
