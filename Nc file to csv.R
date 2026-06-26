# 1. Install and load required packages to convert nc.file to csv file
if(!require("terra")) install.packages("terra")
if(!require("dplyr")) install.packages("dplyr")
if(!require("purrr")) install.packages("purrr") # for multiple year analysis

library(terra)
library(dplyr)
library(purrr)

# 2. Define the temporal range and saved file name - change based on country
years <- 2010:2025
output_file <- "Cadiz_daily_averages_2010_2025.csv"

# 3. Include and process each year using a loop/map function to make one table
compiled_data <- map_df(years, function(yr) {

  # Make file names for temp, precip, and radiation.
  temp_file   <- paste0("T_Cadiz_", yr, ".nc")
  precip_file <- paste0("P_Cadiz_", yr, ".nc")
  rad_file    <- paste0("R_Cadiz_", yr, ".nc")

  # Check if all three files exist before processing to avoid errors
  if(!file.exists(temp_file) | !file.exists(precip_file) | !file.exists(rad_file)) {
    warning(paste("Missing files for year:", yr, "- Skipping."))
    return(NULL)
  }

  # Load the raster data
  r_temp   <- rast(temp_file)
  r_precip <- rast(precip_file)
  r_rad    <- rast(rad_file)

  # Extract daily means from the hourly readings taken from selected 0.25 degree grid squares. Using unweighted mean.
  year_data <- data.frame(
    date_time   = time(r_temp),
    temperature = global(r_temp, fun = "mean", na.rm = TRUE)[, 1],
    precip      = global(r_precip, fun = "mean", na.rm = TRUE)[, 1],
    radiation   = global(r_rad, fun = "mean", na.rm = TRUE)[, 1]
  ) %>%
    # Convert timestamps to dates and aggregate by day (remove hours) in table.
    mutate(day = as.Date(date_time)) %>%
    group_by(day) %>%
    summarise(
      year            = yr,
      avg_temp_C      = mean(temperature, na.rm = TRUE) - 273.15, # Convert Kelvin to Celsius
      avg_rad_Wm2     = mean(radiation, na.rm = TRUE),            # Mean radiation
      total_precip_mm = sum(precip, na.rm = TRUE) * 1000,         # Metres convert to daily mm
      .groups         = "drop"
    )

  return(year_data)
})

# 4. View the first 6 rows of final table structure to see if correct
print(head(compiled_data))
print(tail(compiled_data))

# 5. Save the complete dataset to laptop
write.csv(compiled_data, output_file, row.names = FALSE)
