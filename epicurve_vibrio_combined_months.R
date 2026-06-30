#loads dplyr for cleaning and ggplot2 for plotting
library(tidyverse)

#Read case data
Outcome.Data.Vibrio <- read.csv("Outcome Data Vibrio.csv", stringsAsFactors = FALSE)

# drop any row with missing year, month, or area data.
cleaned_data <- Outcome.Data.Vibrio %>%
  filter(!is.na(year) & !is.na(month) & !is.na(area)) %>%
  mutate(
    #Convert columns safely to numbers to remove any typos and close modification function
    Year_num  = as.numeric(as.character(year)),
    Month_num = as.numeric(as.character(month))
  ) %>%
  filter(Month_num >= 1 & Month_num <= 12) %>%

  # Use a year (2020) to create a one shared 12-month timeline axis
  mutate(plot_date = make_date(year = 2020, month = Month_num, day = 1)) %>%
  filter(!is.na(plot_date))


# combines all years together for each month, by each area
plot_counts <- cleaned_data %>%
  count(area, plot_date, name = "case_counts")


#epicurve creation

#Create the plot with month on X axis and cases on Y axis
epi_curve_plot <- ggplot(plot_counts, aes(x = plot_date, y = case_counts)) +

  # Seasonal background panels with different grey shades
  geom_rect(aes(xmin = make_date(2019, 12, 15), xmax = make_date(2020, 3, 31), ymin = -Inf, ymax = Inf), fill = "gray96") +
  geom_rect(aes(xmin = make_date(2020, 4, 1), xmax = make_date(2020, 6, 30), ymin = -Inf, ymax = Inf), fill = "gray86") +
  geom_rect(aes(xmin = make_date(2020, 7, 1), xmax = make_date(2020, 9, 30), ymin = -Inf, ymax = Inf), fill = "gray70") + #summer
  geom_rect(aes(xmin = make_date(2020, 10, 1), xmax = make_date(2020, 12, 31), ymin = -Inf, ymax = Inf), fill = "gray86") +


  # Monthly timeline data bars (combined multi-year totals) coloured blue
  geom_col(fill = "steelblue", width = 25) +

  # Stacks 3 wide rows vertically (one row for each city)
  # Removes the year columns and gives the months more room
  facet_wrap(~ area, ncol = 1, scales = "free_y") +

  # Filters labels
  # This expands the timeline to full year edges even if your CSV completely lacks Jan or Dec rows, with limits.
  scale_x_date(
    breaks = make_date(2020, c(1, 4, 7, 10), 1),
    labels = c("Jan", "Apr", "Jul", "Oct"),
    limits = c(make_date(2019, 12, 15), make_date(2020, 12, 31)), # FIX: Prevents January/December clipping
    expand = c(0, 0)
  ) +


  scale_y_continuous(breaks = scales::pretty_breaks()) +

  #layout detail
  theme_bw() +

  # Formatting settings for clean, wide rows, with a tilt for readability
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),

    # Formats the 3 city titles at the top of each wide bar panel
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    strip.background = element_rect(fill = "gray90", color = "gray70"),

    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),

    plot.caption = element_text(hjust = 0, size = 10, face = "plain", margin = margin(t = 15))
  ) +

  # Titles and Figure Caption
  labs(
    title = "Seasonal Distribution of locally acquired Non-Cholera Vibrio Cases (2010-2025 Combined), in Cádiz, Charente-maritime, and Stockholm",
    x = "Month of Year",
    y = "Total Combined Case Counts",
    caption = "Figure 1: Combined monthly distribution of non-cholera Vibrio case admissions notified between 2010 and 2025.\nData is combined across all study years to show the regional seasonality profiles. Shading panels define four seasons (Dark grey field = Summer risk peak)."
  )

# Display the plot
print(epi_curve_plot)

# Export with a proportional layout (wider than tall)
ggsave("vibrio_combined_seasonal_curves.png", width = 12, height = 9, dpi = 300)
