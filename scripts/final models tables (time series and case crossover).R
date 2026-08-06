#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# 11. Case-Crossover Univariable Models
# FDR Adjustment for p-value
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(dplyr)
library(survival)
library(officer)
library(flextable)

cc_clean_univar <- cc_data_lags
target_lags     <- 0:3
unique_areas    <- c("cadiz", "charente-maritime", "stockholm", "arcachon")
cc_univar_raw   <- list()

format_p <- function(p) {
  if (is.na(p)) return("-")
  if (p < 0.001) return("< 0.001")
  return(sprintf("%.3f", p))
}

# Include 5 parameters
cc_univar_configs <- list(
  list(prefix = "hi18_lag",   label = "SST (>18°C Threshold, per 1°C increase)", type = "threshold"),
  list(prefix = "Tmean_lag_", label = "Air Temperature (Continuous 1 SD increase)", type = "standardized"),
  list(prefix = "sal_lag",    label = "Sea Surface Salinity (Continuous 1 SD increase)", type = "threshold"),
  list(prefix = "Pmean_lag_", label = "Precipitation (Continuous 1 SD increase)", type = "standardized"),
  list(prefix = "Rmean_lag_", label = "Relative Humidity (Continuous 1 SD increase)", type = "standardized")
)

# 1. Apply independent univariate regressions across area strata
for (current_area in unique_areas) {
  for (config in cc_univar_configs) {
    for (current_lag in target_lags) {

      target_var <- paste0(config$prefix, current_lag)

      cc_sub <- cc_clean_univar %>%
        filter(area == current_area & !is.na(.data[[target_var]]) & !is.na(id)) %>%
        group_by(id) %>%
        filter(any(isb == 1) & any(isb == 0)) %>%
        ungroup()

      if (nrow(cc_sub) < 10) next
      actual_case_count <- sum(cc_sub$isb == 1, na.rm = TRUE)

      if (config$type == "threshold") {
        cc_sub$transformed_x <- cc_sub[[target_var]]
      } else {
        cc_sub$transformed_x <- as.numeric(scale(cc_sub[[target_var]]))
      }

      fit <- tryCatch(
        clogit(isb ~ transformed_x + strata(id), data = cc_sub),
        error = function(e) NULL
      )
      if (is.null(fit)) next

      co <- summary(fit)$coefficients
      if (!"transformed_x" %in% rownames(co)) next

      cc_univar_raw[[length(cc_univar_raw) + 1]] <- data.frame(
        Region   = toupper(current_area),
        Lag      = paste0("Lag ", current_lag),
        Cases    = actual_case_count,
        Variable = config$label,
        b        = co["transformed_x", "coef"],
        se       = co["transformed_x", "se(coef)"],
        p_raw    = co["transformed_x", "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
    }
  }
}

# 2. Put results in table, p-value adjustment
cc_univar_results <- do.call(rbind, cc_univar_raw) %>%
  mutate(
    Region   = factor(Region, levels = c("CADIZ", "CHARENTE-MARITIME", "STOCKHOLM", "ARCACHON")),
    Variable = factor(Variable, levels = c(
      "SST (>18°C Threshold, per 1°C increase)",
      "Air Temperature (Continuous 1 SD increase)",
      "Sea Surface Salinity (Continuous 1 SD increase)",
      "Precipitation (Continuous 1 SD increase)",
      "Relative Humidity (Continuous 1 SD increase)"
    ))
  ) %>%
  arrange(Region, Variable, Lag) %>%
  group_by(Region) %>%                                          # <-- CHANGED: region-level FDR family
  mutate(p_adjusted = p.adjust(p_raw, method = "BH")) %>%
  ungroup() %>%
  mutate(
    OR_CI         = paste0(sprintf("%.2f", exp(b)), " (", sprintf("%.2f", exp(b - 1.96*se)), ", ", sprintf("%.2f", exp(b + 1.96*se)), ")"),
    P_Raw_Display = sapply(p_raw, format_p),
    P_FDR_Display = sapply(p_adjusted, format_p),
    sig           = p_adjusted < 0.05
  )

# make two tables; primary variables, and additional (tmean, rmean)

# Table A: Primary Environmental/Marine Exposures (SST, Salinity, Rainfall)
cc_results_primary <- cc_univar_results %>%
  filter(Variable %in% c("SST (>18°C Threshold, per 1°C increase)",
                         "Sea Surface Salinity (Continuous 1 SD increase)",
                         "Precipitation (Continuous 1 SD increase)")) %>%
  mutate(Variable = droplevels(Variable))

# Table B: Secondary Atmospheric/Weather Exposures (Air Temp, Relative Humidity)
cc_results_secondary <- cc_univar_results %>%
  filter(Variable %in% c("Air Temperature (Continuous 1 SD increase)",
                         "Relative Humidity (Continuous 1 SD increase)")) %>%
  mutate(Variable = droplevels(Variable))


#make table of primary variables
ft_cc_primary <- flextable(cc_results_primary %>% select(Region, Lag, Cases, Variable, OR_CI, P_Raw_Display, P_FDR_Display)) %>%
  font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 9.5, part = "all") %>%
  bold(part = "header") %>% align(align = "center", j = c("Lag", "Cases", "OR_CI", "P_Raw_Display", "P_FDR_Display"), part = "all") %>%
  align(align = "left", j = c("Region", "Variable"), part = "all") %>% merge_v(j = c("Region")) %>% border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 1.5)) %>% hline_bottom(border = fp_border(color = "black", width = 1.5)) %>%
  hline(i = seq(12, nrow(cc_results_primary), by = 12), border = fp_border(color = "black", width = 1.0)) %>%
  set_header_labels(Region = "Study Region", Lag = "Exposure Lag", Cases = "Cases (N)", Variable = "Primary Environmental Predictor",
                    OR_CI = "Crude Odds Ratio (95% CI)", P_Raw_Display = "Raw p-value", P_FDR_Display = "FDR-Adjusted p-value") %>%
  bold(i = which(cc_results_primary$sig), j = c("OR_CI", "P_FDR_Display")) %>% autofit()


# make table of secondary variables
ft_cc_secondary <- flextable(cc_results_secondary %>% select(Region, Lag, Cases, Variable, OR_CI, P_Raw_Display, P_FDR_Display)) %>%
  font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 9.5, part = "all") %>%
  bold(part = "header") %>% align(align = "center", j = c("Lag", "Cases", "OR_CI", "P_Raw_Display", "P_FDR_Display"), part = "all") %>%
  align(align = "left", j = c("Region", "Variable"), part = "all") %>% merge_v(j = c("Region")) %>% border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 1.5)) %>% hline_bottom(border = fp_border(color = "black", width = 1.5)) %>%
  hline(i = seq(8, nrow(cc_results_secondary), by = 8), border = fp_border(color = "black", width = 1.0)) %>%
  set_header_labels(Region = "Study Region", Lag = "Exposure Lag", Cases = "Cases (N)", Variable = "Secondary Atmospheric Predictor",
                    OR_CI = "Crude Odds Ratio (95% CI)", P_Raw_Display = "Raw p-value", P_FDR_Display = "FDR-Adjusted p-value") %>%
  bold(i = which(cc_results_secondary$sig), j = c("OR_CI", "P_FDR_Display")) %>% autofit()


#make word doc
doc_combined_univar <- read_docx() %>%
  body_add_par("Table CC.1-A: Stratified Single-Variable Case-Crossover Odds Ratios for Primary Environmental Predictors Across Lags", style = "heading 1") %>%
  body_add_flextable(ft_cc_primary) %>%
  body_add_par("\nNote: Bold entries indicate statistical significance at the False Discovery Rate (BH) adjusted p < 0.05 level calculated locally per distinct variable baseline array within each specific region block.\n\n", style = "Normal") %>%

  body_add_par("Table CC.1-B: Stratified Single-Variable Case-Crossover Odds Ratios for Secondary Atmospheric Predictors Across Lags", style = "heading 1") %>%
  body_add_flextable(ft_cc_secondary) %>%
  body_add_par("\nNote: Bold entries indicate statistical significance at the Benjamini-Hochberg FDR-adjusted p < 0.05 level calculated locally per regional parameter family framework.", style = "Normal")

print(doc_combined_univar, target = "MSc_Thesis_Table_CaseCrossover_Univariate.docx")
cat("\n[SUCCESS]: Combined Univar Word Document containing both Table CC.1-A and Table CC.1-B generated successfully!\n")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# MODEL 2: Case-Crossover Multivariable Model
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
cc_clean_triple     <- cc_data_lags
triple_model_output <- data.frame()

for (current_area in unique_areas) {
  for (current_lag in target_lags) {
    v_sst   = paste0("hi18_lag", current_lag)
    v_sss   = paste0("sal_lag", current_lag)
    v_pmean = paste0("Pmean_lag_", current_lag)

    triple_subset = cc_clean_triple %>%
      filter(area == current_area & !is.na(.data[[v_sst]]) & !is.na(.data[[v_sss]]) & !is.na(.data[[v_pmean]]) & !is.na(id)) %>%
      group_by(id) %>% filter(any(isb == 1) & any(isb == 0)) %>% ungroup()

    if (nrow(triple_subset) < 10) next
    actual_case_count = sum(triple_subset$isb == 1, na.rm = TRUE)

    triple_subset$SST_threshold_var = triple_subset[[v_sst]]
    triple_subset$sss_scaled_vector = as.numeric(scale(triple_subset[[v_sss]]))
    triple_subset$pmean_scaled_1sd  = as.numeric(scale(triple_subset[[v_pmean]]))

    mod_triple = clogit(isb ~ SST_threshold_var + sss_scaled_vector + pmean_scaled_1sd + strata(id), data = triple_subset)
    sum_triple = summary(mod_triple)$coefficients

    exposures = c("SST_threshold_var", "sss_scaled_vector", "pmean_scaled_1sd")
    labels    = c("SST (>18°C Threshold, per 1°C increase)", "Sea Surface Salinity (Continuous 1 SD increase)", "Precipitation (Continuous 1 SD increase)")

    for (k in seq_along(exposures)) {
      exp_name = exposures[k]
      b        = sum_triple[exp_name, "coef"]
      se       = sum_triple[exp_name, "se(coef)"]
      p_raw    = sum_triple[exp_name, "Pr(>|z|)"]
      or_txt   = paste0(round(exp(b), 2), " (", round(exp(b - 1.96*se), 2), ", ", round(exp(b + 1.96*se), 2), ")")

      triple_model_output = rbind(triple_model_output, data.frame(
        Region = toupper(current_area), Lag = paste0("Lag ", current_lag), Cases = actual_case_count,
        Variable = labels[k], OR = or_txt, p_raw = p_raw, stringsAsFactors = FALSE
      ))
    }
  }
}

#define areas and adjustment of p-value
triple_model_output = triple_model_output %>%
  mutate(
    Region   = factor(Region, levels = c("CADIZ", "CHARENTE-MARITIME", "STOCKHOLM", "ARCACHON")),
    Variable = factor(Variable, levels = c(
      "SST (>18°C Threshold, per 1°C increase)", "Sea Surface Salinity (Continuous 1 SD increase)", "Precipitation (Continuous 1 SD increase)"
    ))
  ) %>%
  arrange(Region, Variable, Lag) %>%
  group_by(Region) %>%                                          # <-- CHANGED: region-level FDR family
  mutate(p_adjusted = p.adjust(p_raw, method = "BH")) %>%
  ungroup() %>%
  mutate(
    P_Raw_Display = sapply(p_raw, format_p),
    P_FDR_Display = sapply(p_adjusted, format_p),
    sig = p_adjusted < 0.05
  )

ft_triple_table = flextable(triple_model_output %>% select(Region, Lag, Cases, Variable, OR, P_Raw_Display, P_FDR_Display)) %>%
  font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 9.5, part = "all") %>%
  bold(part = "header") %>% align(align = "center", j = c("Lag", "Cases", "OR", "P_Raw_Display", "P_FDR_Display"), part = "all") %>%
  align(align = "left", j = c("Region", "Variable"), part = "all") %>% merge_v(j = c("Region")) %>% border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 1.5)) %>% hline_bottom(border = fp_border(color = "black", width = 1.5)) %>%
  hline(i = seq(12, nrow(triple_model_output), by = 12), border = fp_border(color = "black", width = 1.0)) %>%
  set_header_labels(Region = "Study Region", Lag = "Exposure Lag", Cases = "Cases (N)", Variable = "Environmental Predictor (Adjusted)",
                    OR = "Fully Adjusted Odds Ratio (95% CI)", P_Raw_Display = "Raw p-value", P_FDR_Display = "FDR-Adjusted p-value") %>%
  bold(i = which(triple_model_output$sig), j = c("OR", "P_FDR_Display")) %>% autofit()

#make word doc
doc_output_cc2 = read_docx() %>%
  body_add_par("Table CC.2 REVISED: Fully Adjusted Multi-Exposure Case-Crossover Stratified Sensitivity Matrix", style = "heading 1") %>%
  body_add_flextable(ft_triple_table) %>%
  body_add_par("\nNote: Bold entries indicate significance at the FDR-adjusted p < 0.05 level calculated locally per region-exposure category block.", style = "Normal")
print(doc_output_cc2, target = "MSc_Thesis_Table_CaseCrossover_multivariate_LOCAL_FDR.docx")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# MODEL 3: Time Series univariate model
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
library(dplyr)
library(tidyr)
library(splines)
library(flextable)
library(officer)

#define areas
ts_areas   <- c("cadiz", "charente-maritime", "stockholm", "arcachon")
ts_lag_max <- 3

format_p_ts <- function(p) {
  if (is.na(p)) return("-")
  if (p < 0.001) return("< 0.001")
  return(sprintf("%.3f", round(p, 3)))
}

ts_screen_configs <- list(
  list(prefix = "SST_lag",   label = "SST (>18°C Threshold, per 1°C increase)", type = "threshold"),
  list(prefix = "SSS_lag",   label = "Sea Surface Salinity (Continuous 1 SD increase)", type = "standardized"),
  list(prefix = "Pmean_lag", label = "Precipitation (Continuous 1 SD increase)", type = "standardized")
)
ts_screen_raw <- list()

for (current_area in ts_areas) {
  d        <- region_ts_objects[[current_area]]$data
  trend_df <- region_ts_objects[[current_area]]$trend_df
  mu_sss   <- mean(d$SSSmean, na.rm = TRUE); sd_sss   <- sd(d$SSSmean, na.rm = TRUE)
  mu_pmean <- mean(d$Pmean, na.rm = TRUE);   sd_pmean <- sd(d$Pmean, na.rm = TRUE)

  for (config in ts_screen_configs) {
    for (current_lag in 0:ts_lag_max) {
      raw_col <- paste0(config$prefix, current_lag)
      if (config$type == "threshold") {
        d$screen_x <- pmax(d[[raw_col]] - 18, 0)
      } else if (config$prefix == "SSS_lag") {
        d$screen_x <- as.numeric((d[[raw_col]] - mu_sss) / sd_sss)
      } else {
        d$screen_x <- as.numeric((d[[raw_col]] - mu_pmean) / sd_pmean)
      }

      fit <- tryCatch(glm(Daily_Cases ~ screen_x + ns(time_index, df = trend_df) + sin365 + cos365 + dow,
                          family = quasipoisson(link = "log"), data = d), error = function(e) NULL)
      if (is.null(fit)) next
      co <- summary(fit)$coefficients
      if (!"screen_x" %in% rownames(co)) next

      ts_screen_raw[[length(ts_screen_raw) + 1]] <- data.frame(
        Region = toupper(current_area), Lag = paste0("Lag ", current_lag), Cases = sum(d$Daily_Cases),
        Exposure = config$label, b = co["screen_x", "Estimate"], se = co["screen_x", "Std. Error"],
        p_raw = co["screen_x", "Pr(>|t|)"], stringsAsFactors = FALSE
      )
    }
  }
}

#p-value correction
ts_screen_results <- do.call(rbind, ts_screen_raw) %>%
  mutate(
    Region   = factor(Region, levels = c("CADIZ", "CHARENTE-MARITIME", "STOCKHOLM", "ARCACHON")),
    Exposure = factor(Exposure, levels = c(
      "SST (>18°C Threshold, per 1°C increase)", "Sea Surface Salinity (Continuous 1 SD increase)", "Precipitation (Continuous 1 SD increase)"
    ))
  ) %>%
  arrange(Region, Exposure, Lag) %>%
  group_by(Region) %>%                                          # <-- CHANGED: region-level FDR family
  mutate(p_adjusted = p.adjust(p_raw, method = "BH")) %>%
  ungroup() %>%
  mutate(
    RR_CI         = paste0(sprintf("%.2f", exp(b)), " (", sprintf("%.2f", exp(b - 1.96*se)), ", ", sprintf("%.2f", exp(b + 1.96*se)), ")"),
    P_Raw_Display = sapply(p_raw, format_p_ts),
    P_FDR_Display = sapply(p_adjusted, format_p_ts),
    sig           = p_adjusted < 0.05
  )

ft_ts_screen <- flextable(ts_screen_results %>% select(Region, Lag, Cases, Exposure, RR_CI, P_Raw_Display, P_FDR_Display)) %>%
  font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 9.5, part = "all") %>%
  bold(part = "header") %>% align(align = "center", j = c("Lag", "Cases", "RR_CI", "P_Raw_Display", "P_FDR_Display"), part = "all") %>%
  align(align = "left", j = c("Region", "Exposure"), part = "all") %>% merge_v(j = c("Region")) %>% border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 1.5)) %>% hline_bottom(border = fp_border(color = "black", width = 1.5)) %>%
  hline(i = seq(12, nrow(ts_screen_results), by = 12), border = fp_border(color = "black", width = 1.0)) %>%
  set_header_labels(Region = "Study Region", Lag = "Exposure Lag", Cases = "Total Cases (N)", Exposure = "Environmental Predictor (Single)",
                    RR_CI = "Crude Rate Ratio (95% CI)", P_Raw_Display = "Raw p-value", P_FDR_Display = "FDR-Adjusted p-value") %>%
  bold(i = which(ts_screen_results$sig), j = c("RR_CI", "P_FDR_Display")) %>% autofit()

#word doc
doc_ts_screen <- read_docx() %>%
  body_add_par("Table TS.1 REVISED: Time Series Screening Associations (Single-Exposure, Single-Lag Harmonized Quasi-Poisson Models)", style = "heading 1") %>%
  body_add_flextable(ft_ts_screen) %>%
  body_add_par("\nNote: Bold entries indicate statistical significance at the False Discovery Rate (BH) adjusted p < 0.05 level calculated locally within distinct region-exposure parameters.", style = "Normal")
print(doc_ts_screen, target = "MSc_Thesis_Table_TimeSeries_univariate_LOCAL_FDR.docx")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# MODEL 4: Mutually Adjusted Area-Specific Time Series Model
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
      filter(complete.cases(sst_trans, sss_trans, pmean_trans))

    fit <- glm(Daily_Cases ~ sst_trans + sss_trans + pmean_trans + ns(time_index, df = trend_df) + sin365 + cos365 + dow,
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
          Region   = toupper(current_area), Lag      = paste0("Lag ", current_lag), Cases    = sum(d_lag_fit$Daily_Cases),
          Exposure = pred$label, b        = as.numeric(logRR), se       = as.numeric(seLogRR),
          p_raw    = as.numeric(2 * pnorm(-abs(logRR / seLogRR))), stringsAsFactors = FALSE
        )
      }
    }
  }
}

#p-value adjusted
ts_adjusted_results <- do.call(rbind, ts_adjusted_raw) %>%
  mutate(
    Region   = factor(Region, levels = c("CADIZ", "CHARENTE-MARITIME", "STOCKHOLM", "ARCACHON")),
    Exposure = factor(Exposure, levels = c(
      "SST (>18°C Threshold, per 1°C increase)", "Sea Surface Salinity (Continuous 1 SD increase)", "Precipitation (Continuous 1 SD increase)"
    ))
  ) %>%
  arrange(Region, Exposure, Lag) %>%
  group_by(Region) %>%                                          # <-- CHANGED: region-level FDR family
  mutate(p_adjusted = p.adjust(p_raw, method = "BH")) %>%
  ungroup() %>%
  mutate(
    RR_CI         = paste0(sprintf("%.2f", exp(b)), " (", sprintf("%.2f", exp(b - 1.96*se)), ", ", sprintf("%.2f", exp(b + 1.96*se)), ")"),
    P_Raw_Display = sapply(p_raw, format_p_ts),
    P_FDR_Display = sapply(p_adjusted, format_p_ts),
    sig           = p_adjusted < 0.05
  )

ft_ts_adjusted <- flextable(ts_adjusted_results %>% select(Region, Lag, Cases, Exposure, RR_CI, P_Raw_Display, P_FDR_Display)) %>%
  font(fontname = "Times New Roman", part = "all") %>% fontsize(size = 9.5, part = "all") %>%
  bold(part = "header") %>% align(align = "center", j = c("Lag", "Cases", "RR_CI", "P_Raw_Display", "P_FDR_Display"), part = "all") %>%
  align(align = "left", j = c("Region", "Exposure"), part = "all") %>% merge_v(j = c("Region")) %>% border_remove() %>%
  hline_top(border = fp_border(color = "black", width = 1.5)) %>% hline_bottom(border = fp_border(color = "black", width = 1.5)) %>%
  hline(i = seq(12, nrow(ts_adjusted_results), by = 12), border = fp_border(color = "black", width = 1.0)) %>%
  set_header_labels(Region = "Study Region", Lag = "Exposure Lag", Cases = "Total Cases (N)", Exposure = "Environmental Predictor (Adjusted)",
                    RR_CI = "Fully Adjusted Rate Ratio (95% CI)", P_Raw_Display = "Raw p-value", P_FDR_Display = "FDR-Adjusted p-value") %>%
  bold(i = which(ts_adjusted_results$sig), j = c("RR_CI", "P_FDR_Display")) %>% autofit()

#word doc
doc_ts_adjusted <- read_docx() %>%
  body_add_par("Table TS.2 REVISED: Mutually Adjusted Area-Specific Time Series Model with Standardized Single Lags (Lags 0-3)", style = "heading 1") %>%
  body_add_flextable(ft_ts_adjusted) %>%
  body_add_par("\nNote: Each lag is evaluated within its own independent model with FDR corrections run locally within specific region-by-exposure metrics.", style = "Normal")
print(doc_ts_adjusted, target = "MSc_Thesis_Table_TimeSeries_Adjusted_LOCAL_FDR.docx")
cat("\n[SUCCESS]: All 4 updated model blocks have completed execution successfully!\n")

#lag plots
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# PLOT: Mutually Adjusted Area-Specific Lag-Response Curves (Lags 0-3)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

library(dplyr)
library(ggplot2)
library(gridExtra)

# 1. Reuse the exact clean dataframe generated during your Model 4 run
# We harvest the data directly from the 'ts_adjusted_results' output table
if (!exists("ts_adjusted_results")) {
  stop("Please run the Model 4 script first to generate 'ts_adjusted_results'.")
}

# Ensure correct factor levels for clean plotting layout
plot_data <- ts_adjusted_results %>%
  mutate(
    # Strip "Lag " string to get numeric values for the X-axis
    Lag_Numeric = as.numeric(gsub("Lag ", "", Lag)),
    # Extract lower and upper 95% CI bounds mathematically from b and se
    RR       = exp(b),
    RR_Lower = exp(b - 1.96 * se),
    RR_Upper = exp(b + 1.96 * se),
    # Prettify labels for plot facets
    Exposure_Short = case_when(
      Exposure == "SST (>18°C Threshold, per 1°C increase)" ~ "SST (>18°C)",
      Exposure == "Sea Surface Salinity (Continuous 1 SD increase)" ~ "Salinity (1 SD)",
      Exposure == "Precipitation (Continuous 1 SD increase)" ~ "Precipitation (1 SD)",
      TRUE ~ as.character(Exposure)
    )
  )

# 2. Plotting Function for Consistent Styling Across Regions
generate_lag_curve <- function(data_subset, region_title) {
  ggplot(data_subset, aes(x = Lag_Numeric, y = RR, group = Exposure_Short, color = Exposure_Short)) +
    # Reference line at RR = 1.0 (No Risk Change)
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray50", size = 0.6) +
    # Confidence intervals represented by shaded ribbons
    geom_ribbon(aes(ymin = RR_Lower, ymax = RR_Upper, fill = Exposure_Short), alpha = 0.15, color = NA) +
    # Main continuous trend lines
    geom_line(size = 1.1) +
    # Exact point estimates per discrete lag
    geom_point(size = 2) +
    # Dynamic grid layout splitting the three mutually adjusted predictors
    facet_wrap(~Exposure_Short, scales = "free_y", ncol = 3) +
    # Aesthetic adjustments
    scale_x_continuous(breaks = 0:3, labels = c("Lag 0", "Lag 1", "Lag 2", "Lag 3")) +
    labs(
      title = paste("Mutually Adjusted Lag-Response Curves:", region_title),
      x = "Exposure Lag Time",
      y = "Rate Ratio (95% CI)"
    ) +
    theme_minimal(base_family = "Times New Roman", base_size = 11) +
    theme(
      plot.title        = element_text(face = "bold", size = 12, hjust = 0.5),
      strip.background  = element_rect(fill = "gray95", color = "gray80"),
      strip.text        = element_text(face = "bold", color = "black"),
      legend.position   = "none", # Legend omitted since facets are self-explanatory
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(color = "gray80", fill = NA, size = 0.8)
    )
}

# 3. Generate individual plots for the 4 distinct regions
plots_list <- list()
for (current_area in unique(plot_data$Region)) {
  area_df <- plot_data %>% filter(Region == current_area)
  plots_list[[current_area]] <- generate_lag_curve(area_df, current_area)
}

# 4. Save each plot as a high-resolution image for your thesis layout
# Saves a unified layout configuration per area block
for (current_area in names(plots_list)) {
  file_name <- paste0("MSc_Thesis_LagResponse_Mutually_Adjusted_", current_area, ".png")
  ggsave(
    filename = file_name,
    plot = plots_list[[current_area]],
    width = 8.5,
    height = 3.5,
    dpi = 300
  )
}

# 5. Display the plots sequentially in your R studio window
grid.arrange(grobs = plots_list, ncol = 1)
cat("\n[SUCCESS]: Mutually adjusted lag-response plots saved for all 4 study regions!\n")

