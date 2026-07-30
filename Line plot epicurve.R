# 1. Automated Package Check & Load
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("viridis")) install.packages("viridis")

library(tidyverse)
library(viridis)

# 2. Read CSV file
Outcome.Data.Vibrio <- read.csv("Outcome Data Vibrio.csv", stringsAsFactors = FALSE)

# 3. Clean data
cleaned_data <- Outcome.Data.Vibrio %>%
  filter(!is.na(year) & !is.na(month) & !is.na(area)) %>%
  mutate(
    Year_num  = as.numeric(as.character(year)),
    Month_num = as.numeric(as.character(month))
  ) %>%
  filter(Month_num >= 1 & Month_num <= 12)

# 4. Generate a complete grid
complete_grid <- expand_grid(
  area = unique(cleaned_data$area),
  Year_num = min(cleaned_data$Year_num):max(cleaned_data$Year_num),
  Month_num = 1:12
)

# 5. Calculate monthly case counts per year
plot_counts <- cleaned_data %>%
  count(area, Year_num, Month_num, name = "case_counts") %>%
  right_join(complete_grid, by = c("area", "Year_num", "Month_num")) %>%
  mutate(
    case_counts = replace_na(case_counts, 0),
    Year_factor = as.factor(Year_num)
  )

# extend edges of plot
edge_left <- plot_counts %>% filter(Month_num == 1) %>% mutate(Month_num = 0.5, case_counts = 0)
edge_right <- plot_counts %>% filter(Month_num == 12) %>% mutate(Month_num = 12.5, case_counts = 0)

plot_counts <- bind_rows(plot_counts, edge_left, edge_right) %>%
  arrange(area, Year_num, Month_num)

# 6. colours of lines
unique_years <- sort(unique(plot_counts$Year_num))

highlight_colors <- case_when(
  unique_years == 2025 ~ "#990000",
  unique_years == 2024 ~ "#1D3557",
  unique_years == 2023 ~ "#2A9D8F",
  unique_years == 2022 ~ "#D9A752",
  unique_years == 2021 ~ "#9E9D58",
  unique_years == 2020 ~ "#567D65",
  unique_years == 2019 ~ "#437A80",
  unique_years == 2018 ~ "#466987",
  unique_years == 2017 ~ "#585A8A",
  unique_years == 2016 ~ "#70507D",
  unique_years == 2015 ~ "#80486C",
  unique_years == 2014 ~ "#8A4556",
  unique_years == 2013 ~ "#7A4D4F",
  unique_years == 2012 ~ "#5C5657",
  unique_years == 2011 ~ "#6E7A8A",
  unique_years == 2010 ~ "#8C9BA5",
  TRUE                 ~ "#A9A9A9"
)

#Season key
uk_seasons_df <- data.frame(
  xmin = c(0.5,  2.5,  5.5,  8.5,  11.5),
  xmax = c(2.5,  5.5,  8.5,  11.5, 12.5),
  Season = factor(c("Winter", "Spring", "Summer", "Autumn", "Winter"),
                  levels = c("Winter", "Spring", "Summer", "Autumn"))
)

# 7. Month of Admission on X axis
epi_curve_plot <- ggplot(plot_counts, aes(x = Month_num, y = case_counts, group = Year_factor, color = Year_factor)) +

  # Panel shading for seasons
  geom_rect(data = uk_seasons_df, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Season),
            color = NA, inherit.aes = FALSE) +

  # Trendlines
  geom_line(linewidth = 0.95) +
  geom_point(data = filter(plot_counts, Month_num %% 1 == 0), aes(size = ifelse(Year_num %in% c(2023, 2024, 2025), 2.0, 1.0))) +

  # 8. 4 area rows with proportional heights
  facet_grid(area ~ ., space = "free_y", scales = "free_y") +

  # 9.
  scale_x_continuous(
    breaks = c(1, 4, 7, 10),
    labels = c("Jan", "Apr", "Jul", "Oct"),
    limits = c(0.5, 12.5),
    expand = c(0, 0)
  ) +

  # 10. Y-axis
  scale_y_continuous(
    breaks = scales::breaks_width(1),
    expand = expansion(mult = c(0.03, 0.15))
  ) +

  # 11. scaling vectors
  scale_color_manual(values = highlight_colors) +
  scale_size_identity() +

  # linking colors directly to the season key
  scale_fill_manual(
    values = c("Winter" = "gray96", "Spring" = "gray88", "Summer" = "gray65", "Autumn" = "gray78"),
    name = "Season Key"
  ) +

  # 12. Visual layout styling
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),

    # FIXED: Places region names horizontally and centered on the right margin strips
    strip.text.y = element_text(angle = 0, face = "bold", size = 11, color = "black", hjust = 0.5),
    strip.background = element_rect(fill = "gray90", color = "gray70"),

    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),

    # Legend settings
    legend.position = "right",
    legend.box = "vertical",
    legend.key.width = unit(1.2, "cm"),
    legend.title = element_text(face = "bold", size = 10),

    # Centers the plot title and caption text smoothly
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 10), hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, size = 10, face = "plain", margin = margin(t = 15))
  ) +

  # 13. Wrapped title and caption text strings prevent edge clipping
  labs(
    title = str_wrap("Locally acquired cases of non-cholera vibrio by month of admission and year (2010–2025)", width = 75),
    x = "Month of Admission",
    y = "Monthly Case Counts",
    color = "Year",
    caption = str_wrap("Figure 1: Comparison of non-cholera Vibrio trend lines by calendar year across geographic areas. All individual years are tracked using a professional chronological jewel-tone color profile. Shading panels define UK meteorological seasons.", width = 110)
  )

# show final plot
print(epi_curve_plot)

# Export
ggsave("vibrio_epi_curves_final_professional.png", width = 13, height = 12, dpi = 300)
