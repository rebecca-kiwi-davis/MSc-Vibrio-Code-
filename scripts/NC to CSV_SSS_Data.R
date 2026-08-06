# 1. Load the libraries 
library(terra) #Nc data specific
library(tidyverse) #clean data and plot
library(tidyterra) #map spatial data with ggplot

# 2. Load the Multi Observation Global Ocean Sea Surface Salinity and Sea Surface Density dataset
nc_data <- rast("SSS_France_2010_2025.nc")

# 4. Calculate the daily spatial average across the defined coastal area, ignore empty data points including land
daily_averages <- global(nc_data, fun = "mean", na.rm = TRUE)

# 5. Extract dates and pair them with Sea Surface Salinity (SSS) averages
daily_summary <- tibble(
  date  = time(nc_data),
  SSS_average = daily_averages$mean
)

# 6. View the dataset 
print(daily_summary)

# 7. Save to CSV file
write_csv(daily_summary, "France_daily_average_2010_2025_SSS.csv")

# 8. Plot the daily SSS timeline over the years
ggplot(daily_summary, aes(x = date, y = SSS_average)) +
  geom_line(color = "steelblue", linewidth = 0.5) +
  theme_minimal() +
  labs(
    title = "France Daily Regional SSS Average (2010-2024)",
    x = "Date",
    y = "Average Sea Surface Salinity (SSS)"
  )
