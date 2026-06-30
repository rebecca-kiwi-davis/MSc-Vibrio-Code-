library(tidyverse)

# 2. Read CSV file
Outcome.Data.Vibrio <- read.csv("Outcome Data Vibrio.csv", stringsAsFactors = FALSE)

#clean data
cleaned_data <- Outcome.Data.Vibrio %>%

   # Remove rows with missing crucial variables (year, month, area)
  filter(!is.na(Year) & !is.na(Month) & !is.na(area)) %>%

  # Convert columns safely to numbers to remove any typos and close modification function
  mutate(
    Year_num  = as.numeric(as.character(Year)),
    Month_num = as.numeric(as.character(Month))
  ) %>%

  # Align data to valid calendar months (1 to 12)
  filter(Month_num >= 1 & Month_num <= 12) %>%

  # Create a valid monthly timeline date object for plotting using a year
  mutate(plot_date = make_date(year = 2020, month = Month_num, day = 1)) %>%

#Build Epicurve

# Calculate fresh case counts from the new cleaned data (2026 will be gone)
plot_counts <- cleaned_data %>%
  count(area, Year_num, plot_date, name = "case_counts")

# Create the plot with month on X axis and combined cases on Y axis, by area
epi_curve_plot <- ggplot(plot_counts, aes(x = plot_date, y = case_counts)) +

  # Vary the shading of panels depending on month. Light grey (jan-may), light grey (april-june), darker grey (july- sept), light grey (oct-dec)
  geom_rect(aes(xmin = make_date(2020, 1, 1), xmax = make_date(2020, 3, 31), ymin = -Inf, ymax = Inf), fill = "gray96") +
  geom_rect(aes(xmin = make_date(2020, 4, 1), xmax = make_date(2020, 6, 30), ymin = -Inf, ymax = Inf), fill = "gray86") +
  geom_rect(aes(xmin = make_date(2020, 7, 1), xmax = make_date(2020, 9, 30), ymin = -Inf, ymax = Inf), fill = "gray70") + # Summer Focus
  geom_rect(aes(xmin = make_date(2020, 10, 1), xmax = make_date(2020, 12, 31), ymin = -Inf, ymax = Inf), fill = "gray86") + # Autumn

  # Monthly timeline data bars in blue to stand out
  geom_col(fill = "steelblue", width = 25) +

  # 15 panels grid matrix. Splits panels into regions as rows, and years as columns.
  facet_grid(rows = vars(area), cols = vars(Year_num), scales = "free_y") +

  # Filters labels to seasons - Jan, Apr, Jul, Oct
  scale_x_date(
    breaks = make_date(2020, c(1, 4, 7, 10), 1),
    labels = c("Jan", "Apr", "Jul", "Oct"),
    expand = c(0, 0)
  ) +

  # layout style with white background
  theme_bw() +

  # Rotate text to see all and adjust size
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),

    # Formats the titles (cities) on the side to look horizontal, large, and prominent
    strip.text.x = element_text(face = "bold", size = 10),
    strip.text.y = element_text(face = "bold", size = 12, color = "black", angle = 0),
    strip.background.y = element_rect(fill = "gray90", color = "gray70"),

    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),

    # Formats the caption underneath to be left-aligned and legible
    plot.caption = element_text(hjust = 0, size = 10, face = "plain", margin = margin(t = 15))
  ) +

  # Labels and Caption updated to reflect a 2010–2025 timeline
  labs(
    title = "Locally Acquired Cases of Non-Cholera Vibrio by Month of Admission, Cádiz, Charente-Maritime, and Stockholm (2010–2025)",
    x = "Month of Onset",
    y = "Monthly Case Counts",
    caption = "Figure 1: Date of non-cholera Vibrio case admission (mostly hospital) notified between 2010 and 2025.\nShading of background panels define four seasons (Dark grey field = July–September, summer risk peak)."
  )

# Display final plot
print(epi_curve_plot)

# Export as your file
ggsave("vibrio_epi_curves_seasonal_grid.png", width = 22, height = 9, dpi = 300)
