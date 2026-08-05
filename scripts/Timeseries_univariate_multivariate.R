
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 18. Time series
# regional standardization (1 SD Z-scores ) to ensure
# fair comparison across areas.
# temporal Lags of 0-3 are used to match the Case-Crossover model used prior
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


library(dplyr)
library(tidyr)
library(splines)
library(flextable)
library(officer)
library(ggplot2)

#include all areas
ts_areas   <- c("cadiz", "charente-maritime", "stockholm", "arcachon")
ts_lag_max <- 3

format_p_ts <- function(p) {
  if (is.na(p)) return("-")
  if (p < 0.001) return("< 0.001")
  return(sprintf("%.3f", round(p, 3)))
}

# timelines
build_region_ts <- function(current_area) {

  daily_cases <- cases %>%
    mutate(clean_area = trimws(tolower(as.character(area)))) %>%
    filter(isb == 1, clean_area == current_area) %>%
    count(date, name = "Daily_Cases")

  daily_weather <- exposure %>%
    mutate(clean_area = trimws(tolower(as.character(area)))) %>%
    filter(clean_area == current_area) %>%
    distinct(date, SSTmean, SSSmean, Pmean) %>%
    arrange(date)

  full_calendar <- data.frame(date = seq(min(daily_weather$date), max(daily_weather$date), by = "day"))

  ts_data <- full_calendar %>%
    left_join(daily_cases, by = "date") %>%
    mutate(Daily_Cases = ifelse(is.na(Daily_Cases), 0, Daily_Cases)) %>%
    left_join(daily_weather, by = "date") %>%
    arrange(date) %>%
    filter(!is.na(SSTmean), !is.na(SSSmean), !is.na(Pmean)) %>%
    mutate(
      dow        = factor(weekdays(date)),
      time_index = as.numeric(date - min(date)),
      doy        = as.numeric(format(date, "%j")),
      sin365     = sin(2 * pi * doy / 365.25),
      cos365     = cos(2 * pi * doy / 365.25)
    )

  # make lags
  for (l in 0:ts_lag_max) {
    ts_data[[paste0("SST_lag", l)]]   <- dplyr::lag(ts_data$SSTmean, l)
    ts_data[[paste0("SSS_lag", l)]]   <- dplyr::lag(ts_data$SSSmean, l)
    ts_data[[paste0("Pmean_lag", l)]] <- dplyr::lag(ts_data$Pmean, l)
  }

  lag_cols <- c(paste0("SST_lag", ts_lag_max), paste0("SSS_lag", ts_lag_max), paste0("Pmean_lag", ts_lag_max))
  ts_data  <- ts_data[complete.cases(ts_data[, lag_cols]), ]

  n_years   <- as.numeric(diff(range(ts_data$date))) / 365.25
  trend_df  <- max(3, round(n_years / 2))

  list(data = ts_data, trend_df = trend_df)
}

region_ts_objects <- setNames(lapply(ts_areas, build_region_ts), ts_areas)

# 2. scalining for each area
cat("\n=== REGIONAL BASELINE CONTEXTS FOR ANOMALY SCALING ===\n")
for (a in ts_areas) {
  d_context <- region_ts_objects[[a]]$data
  cat(sprintf("[%s] -> Data Rows: %d | Cases: %d\n", toupper(a), nrow(d_context), sum(d_context$Daily_Cases)))
}

# naming
ts_screen_configs <- list(
  list(prefix = "SST_lag",   label = "SST (>18°C Threshold, per 1°C increase)", type = "threshold"),
  list(prefix = "SSS_lag",   label = "Sea Surface Salinity (Continuous 1 SD increase)", type = "standardized"),
  list(prefix = "Pmean_lag", label = "Precipitation (Continuous 1 SD increase)", type = "standardized")
)

ts_screen_raw <- list()

# 1. Run independent regressions across individual areas 
for (current_area in ts_areas) {
  d        <- region_ts_objects[[current_area]]$data
  trend_df <- region_ts_objects[[current_area]]$trend_df

  mu_sss   <- mean(d$SSSmean, na.rm = TRUE); sd_sss   <- sd(d$SSSmean, na.rm = TRUE)
  mu_pmean <- mean(d$Pmean, na.rm = TRUE);   sd_pmean <- sd(d$Pmean, na.rm = TRUE)

  for (config in ts_screen_configs) {
    for (current_lag in 0:ts_lag_max) {

      raw_col <- paste0(config$prefix, current_lag)

      #  local regional standardisation
      if (config$type == "threshold") {
        d$screen_x <- pmax(d[[raw_col]] - 18, 0)
      } else if (config$prefix == "SSS_lag") {
        d$screen_x <- as.numeric((d[[raw_col]] - mu_sss) / sd_sss)
      } else {
        d$screen_x <- as.numeric((d[[raw_col]] - mu_pmean) / sd_pmean)
      }

      fit <- tryCatch(
        glm(Daily_Cases ~ screen_x + ns(time_index, df = trend_df) + sin365 + cos365 + dow,
            family = quasipoisson(link = "log"), data = d),
        error = function(e) NULL
      )
      if (is.null(fit)) next

      co <- summary(fit)$coefficients
      if (!"screen_x" %in% rownames(co)) next

      ts_screen_raw[[length(ts_screen_raw) + 1]] <- data.frame(
        Region = toupper(current_area), Lag = paste0("Lag ", current_lag),
        Cases = sum(d$Daily_Cases), Exposure = config$label,
        b = co["screen_x", "Estimate"], se = co["screen_x", "Std. Error"], p_raw = co["screen_x", "Pr(>|t|)"],
        stringsAsFactors = FALSE
      )
    }
  }
}

# 2. apply true global FDR corrections
ts_screen_results <- do.call(rbind, ts_screen_raw) %>%
  mutate(
    p_adjusted = p.adjust(p_raw, method = "BH"),
    sig        = p_adjusted < 0.05,
    RR_CI      = paste0(sprintf("%.2f", exp(b)), " (", sprintf("%.2f", exp(b - 1.96*se)), ", ", sprintf("%.2f", exp(b + 1.96*se)), ")"),
    P_Value    = sapply(p_adjusted, format_p_ts)
  )

# 3. Structural layout (48 rows across the 4 regions)
ft_ts_screen <- flextable(ts_screen_results %>% select(Region, Lag, Cases, Exposure, RR_CI, P_Value)) %>%
  font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>% align(align = "center", j = c("Lag", "Cases", "RR_CI", "P_Value"), part = "all") %>%
  align(align = "left", j = c("Region", "Exposure"), part = "all") %>% merge_v(j = c("Region")) %>% border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 1.5)) %>% hline_bottom(border = fp_border(color = "black", width = 1.5)) %>%
  hline(i = seq(12, nrow(ts_screen_results), by = 12), border = fp_border(color = "black", width = 1.0)) %>%
  set_header_labels(Region = "Study Region", Lag = "Exposure Lag", Cases = "Total Cases (N)",
                    Exposure = "Environmental Predictor (Single)", RR_CI = "Crude Rate Ratio (95% CI)", P_Value = "p-value") %>%
  bold(i = which(ts_screen_results$sig), j = c("RR_CI", "P_Value")) %>% autofit()

# 4. word doc
doc_ts_screen <- read_docx() %>%
  body_add_par("Table TS.1: Time Series Screening Associations (Single-Exposure, Single-Lag Harmonized Quasi-Poisson Models)", style = "heading 1") %>%
  body_add_flextable(ft_ts_screen) %>%
  body_add_par("\nNote: Bold entries indicate statistical significance at the False Discovery Rate (BH) adjusted p < 0.05 level. Each parameter row represents a completely standalone model. SST is tracked above an absolute 18°C threshold. Salinity and precipitation are locally standardized per region to continuous 1 SD Z-score scales before regression. The analysis spans Lags 0-3.", style = "Normal")

print(doc_ts_screen, target = "MSc_Thesis_Table_TimeSeries_univariate.docx")
cat("\n[SUCCESS]: Original stratified table format successfully restored and exported with Arcachon included!\n")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# STEP 3: Area-Specific, Single-Lag Mutually Adjusted Time Series Model
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(dplyr)
library(splines)
library(flextable)
library(officer)

ts_adjusted_raw   <- list()
ts_model_objects  <- list()

for (current_area in ts_areas) {
  d        <- region_ts_objects[[current_area]]$data
  trend_df <- region_ts_objects[[current_area]]$trend_df

  ts_model_objects[[current_area]] <- list()

  sd_sss   <- sd(d$SSSmean, na.rm = TRUE)
  sd_pmean <- sd(d$Pmean, na.rm = TRUE)

  for (current_lag in 0:ts_lag_max) {

    d_lag_fit <- d %>%
      mutate(
        sst_raw   = dplyr::lag(SSTmean, current_lag),
        sss_raw   = dplyr::lag(SSSmean, current_lag),
        pmean_raw = dplyr::lag(Pmean,   current_lag),

        sst_trans   = pmax(sst_raw - 18, 0),
        sss_trans   = as.numeric((sss_raw - mean(d$SSSmean, na.rm = TRUE)) / sd_sss),
        pmean_trans = as.numeric((pmean_raw - mean(d$Pmean, na.rm = TRUE)) / sd_pmean)
      ) %>%
      filter(complete.cases(sst_trans, strans = sss_trans, pmean_trans))

    fit <- glm(Daily_Cases ~ sst_trans + sss_trans + pmean_trans +
                 ns(time_index, df = trend_df) + sin365 + cos365 + dow,
               family = quasipoisson(link = "log"), data = d_lag_fit)

    model_label <- paste0("lag_", current_lag)
    ts_model_objects[[current_area]][[model_label]] <- fit

    predictors <- list(
      list(coef_name = "sst_trans",   label = "SST (>18°C Threshold, per 1°C increase)"),
      list(coef_name = "sss_trans",   label = "Sea Surface Salinity (Continuous 1 SD increase)"),
      list(coef_name = "pmean_trans", label = "Precipitation (Continuous 1 SD increase)")
    )
    for (pred in predictors) {
      if (pred$coef_name %in% names(coef(fit))) {
        logRR   <- coef(fit)[pred$coef_name]
        seLogRR <- sqrt(vcov(fit)[pred$coef_name, pred$coef_name])

        ts_adjusted_raw[[length(ts_adjusted_raw) + 1]] <- data.frame(
          Region   = toupper(current_area),
          Lag      = paste0("Lag ", current_lag),
          Cases    = sum(d_lag_fit$Daily_Cases),
          Exposure = pred$label,
          b        = as.numeric(logRR),
          se       = as.numeric(seLogRR),
          p_raw    = as.numeric(2 * pnorm(-abs(logRR / seLogRR))),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

#  geographic and variable factor boundaries
ts_adjusted_results <- do.call(rbind, ts_adjusted_raw) %>%
  mutate(
    Region   = factor(Region, levels = c("CADIZ", "CHARENTE-MARITIME", "STOCKHOLM", "ARCACHON")),
    Exposure = factor(Exposure, levels = c(
      "SST (>18°C Threshold, per 1°C increase)",
      "Sea Surface Salinity (Continuous 1 SD increase)",
      "Precipitation (Continuous 1 SD increase)"
    ))
  ) %>%
  arrange(Region, Exposure, Lag) %>%
  mutate(
    p_adjusted = p.adjust(p_raw, method = "BH"),
    sig        = p_adjusted < 0.05,
    RR_CI      = paste0(sprintf("%.2f", exp(b)), " (", sprintf("%.2f", exp(b - 1.96*se)), ", ", sprintf("%.2f", exp(b + 1.96*se)), ")"),
    P_Value    = sapply(p_adjusted, format_p_ts)
  )

ft_ts_adjusted <- flextable(ts_adjusted_results %>% select(Region, Lag, Cases, Exposure, RR_CI, P_Value)) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  align(align = "center", j = c("Lag", "Cases", "RR_CI", "P_Value"), part = "all") %>%
  align(align = "left", j = c("Region", "Exposure"), part = "all") %>%
  merge_v(j = c("Region")) %>%
  border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 1.5)) %>%
  hline_bottom(border = fp_border(color = "black", width = 1.5)) %>%
  hline(i = seq(12, nrow(ts_adjusted_results), by = 12), border = fp_border(color = "black", width = 1.0)) %>%
  set_header_labels(Region = "Study Region", Lag = "Exposure Lag", Cases = "Total Cases (N)",
                    Exposure = "Environmental Predictor (Adjusted)", RR_CI = "Fully Adjusted Rate Ratio (95% CI)", P_Value = "p-value") %>%
  bold(i = which(ts_adjusted_results$sig), j = c("RR_CI", "P_Value")) %>%
  autofit()

doc_ts_adjusted <- read_docx() %>%
  body_add_par("Table TS.2: Mutually Adjusted Area-Specific Time Series Model with Standardized Single Lags (Lags 0-3)", style = "heading 1") %>%
  body_add_flextable(ft_ts_adjusted) %>%
  body_add_par("\nNote: Bold entries indicate statistical significance at the Benjamini-Hochberg False Discovery Rate (BH) adjusted p < 0.05 level. Each single lag distance (Lag 0 through Lag 3) is evaluated within its own independent, area-specific model to prevent multi-collinearity inflation. Models jointly adjust for all three exposures on that specified lag day, along with long-term trend, annual seasonality, and day-of-week via quasi-Poisson regression. Sea Surface Temperature is modeled as an absolute linear threshold effect above 18°C. Sea Surface Salinity and Precipitation are standardized locally using Z-scores to represent risk change per 1 regional Standard Deviation (SD) increase.", style = "Normal")

print(doc_ts_adjusted, target = "MSc_Thesis_Table_TimeSeries_Adjusted.docx")
cat("\n[SUCCESS]: Adjusted table finalized with customized SST -> SSS -> P sorting layout!\n")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# STEP 4: Single-Lag-response plots (per area, per exposure)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


plot_single_lags <- function(region_models, lag_max, transformed_var_name, main_title, y_label) {
  lags  <- 0:lag_max
  RR    <- numeric(length(lags))
  lower <- numeric(length(lags))
  upper <- numeric(length(lags))

  for (i in seq_along(lags)) {
    l <- lags[i]
    model_label <- paste0("lag_", l)
    fit <- region_models[[model_label]]

    # Maps 
    coef_name <- paste0(transformed_var_name, "_trans")

    if (!is.null(fit) && coef_name %in% names(coef(fit))) {
      logRR   <- coef(fit)[coef_name]
      se_log  <- sqrt(vcov(fit)[coef_name, coef_name])

      RR[i]    <- exp(logRR)
      lower[i] <- exp(logRR - 1.96 * se_log)
      upper[i] <- exp(logRR + 1.96 * se_log)
    } else {
      RR[i]    <- 1; lower[i] <- 1; upper[i] <- 1
    }
  }

  plot(lags, RR, type = "b", pch = 19, ylim = range(c(lower, upper, 1)),
       xlab = "Lag (days)", ylab = y_label, main = main_title, xaxt = "n")
  axis(1, at = 0:lag_max, labels = 0:lag_max)
  arrows(lags, lower, lags, upper, angle = 90, code = 3, length = 0.05)
  abline(h = 1, lty = 2, col = "gray50")
}


for (current_area in ts_areas) {
  region_models <- ts_model_objects[[current_area]]

  # --- 1. Generate the Multipanel Single-Lag Response Curve Figures ---
  png(paste0("MSc_Thesis_Figure_TimeSeries_LagResponse_", current_area, ".png"), width = 2400, height = 800, res = 200)
  par(mfrow = c(1, 3), family = "Times New Roman")

  plot_single_lags(region_models, ts_lag_max, "sst",
                   paste0(toupper(current_area), ": SST (>18°C)"), "RR (per 1°C above 18°C)")
  plot_single_lags(region_models, ts_lag_max, "sss",
                   paste0(toupper(current_area), ": SSS (1 SD Anomaly)"), "RR (per 1 Local SD increase)")
 
  plot_single_lags(region_models, ts_lag_max, "pmean",
                   paste0(toupper(current_area), ": Precipitation (1 SD Anomaly)"), "RR (per 1 Local SD increase)")

  dev.off()
}
cat("\n[SUCCESS]: Multipanel single-lag figures cleanly generated for all 4 ordered study regions!\n")





