#1. Load the libraries 
library(terra) #Nc data specific
library(tidyverse) #clean data and plot
library(tidyterra) #map spatial data with ggplot

# 3. Load the Global Ocean Physics Reanalysis Nc file dataset
nc_data <- rast("SST_daily_Stockholm_2010_2025.nc")

# 4. Calculate the daily spatial average across the defined coastal area, ignore empty data points including land
daily_averages <- global(nc_data, fun = "mean", na.rm = TRUE)

# 5. Align dates with Potential Temperature (thetao) averages
daily_summary <- tibble(
  date  = time(nc_data),
  Temp_average = daily_averages$mean
)

# 6. View the final dataset structure by showing first 10 rows
print(daily_summary)

# 7. Save to CSV under new file name
write_csv(daily_summary, "Stockholm_daily_average_2010_2025_SST.csv")

# 8. Plot the daily temperature timeline over the years
ggplot(daily_summary, aes(x = date, y = Temp_average)) +
  geom_line(color = "darkred", linewidth = 0.5) +
  theme_minimal() +
  labs(
    title = "Daily Average Sea Water Potential Temperature",
    x = "Date",
    y = "Average Potential Temperature (°C)"
  )
