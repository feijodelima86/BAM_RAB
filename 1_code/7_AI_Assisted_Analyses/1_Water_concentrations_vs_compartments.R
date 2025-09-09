# Hypothesis 1: Water concentrations influence metal in compartments
# Load required libraries
library(dplyr)
library(corrplot)

# Read data
cd_data <- read.csv("2_incremental/AI_Subsets/cd.csv")


# Clean data - remove NA values for cadmium and compartments
clean_data <- cd_data %>%
  filter(!is.na(Cd) & !is.na(SAMPLE_DESCRIPTOR) & 
           SAMPLE_DESCRIPTOR %in% c("EPIP", "EPIL", "FILA") &
           !is.na(SITE) & !is.na(SAMPLING_DATE))

# Average data by site/date/compartment to avoid pseudoreplication
averaged_data <- clean_data %>%
  group_by(SITE, SAMPLING_DATE, SAMPLE_DESCRIPTOR) %>%
  summarise(
    Cd = mean(Cd, na.rm = TRUE),
    UF_Cd111 = mean(UF_Cd111, na.rm = TRUE),
    UF_Cd114 = mean(UF_Cd114, na.rm = TRUE),
    F_Cd111 = mean(F_Cd111, na.rm = TRUE),
    F_Cd114 = mean(F_Cd114, na.rm = TRUE),
    W_Cd111 = mean(W_Cd111, na.rm = TRUE),
    W_Cd114 = mean(W_Cd114, na.rm = TRUE),
    n_samples = n(),
    .groups = 'drop'
  ) %>%
  # Remove any rows where averaging resulted in NA or infinite values
  filter(is.finite(Cd))

# Focus on main compartments
main_compartments <- c("EPIL", "EPIP", "FILA")
water_cd_columns <- c("UF_Cd111", "UF_Cd114", "F_Cd111", "F_Cd114", "W_Cd111", "W_Cd114")

# 1. Basic statistics for tissue cadmium by compartment
cat("=== TISSUE CADMIUM BY COMPARTMENT (Averaged by Site/Date) ===\n")
tissue_stats <- averaged_data %>%
  filter(SAMPLE_DESCRIPTOR %in% main_compartments) %>%
  group_by(SAMPLE_DESCRIPTOR) %>%
  summarise(
    n_site_date = n(),
    total_samples = sum(n_samples, na.rm = TRUE),
    mean_Cd = mean(Cd, na.rm = TRUE),
    sd_Cd = sd(Cd, na.rm = TRUE),
    median_Cd = median(Cd, na.rm = TRUE),
    min_Cd = min(Cd, na.rm = TRUE),
    max_Cd = max(Cd, na.rm = TRUE),
    .groups = 'drop'
  )
print(tissue_stats)

# 2. Basic statistics for water cadmium concentrations
cat("\n=== WATER CADMIUM CONCENTRATIONS (Averaged by Site/Date) ===\n")
water_stats <- averaged_data %>%
  summarise(across(all_of(water_cd_columns), 
                   list(n = ~sum(!is.na(.) & is.finite(.)),
                        mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE)), 
                   .names = "{.col}_{.fn}"))
print(water_stats)

# 3. Calculate correlations between water and tissue cadmium by compartment
cat("\n=== WATER-TISSUE CADMIUM CORRELATIONS (Site/Date Averaged) ===\n")

correlation_results <- data.frame()

for(comp in main_compartments) {
  comp_data <- averaged_data %>% filter(SAMPLE_DESCRIPTOR == comp)
  
  cat("\n", comp, "compartment:\n")
  
  for(water_col in water_cd_columns) {
    # Filter for complete pairs
    pairs <- comp_data %>%
      filter(!is.na(Cd) & !is.na(.data[[water_col]]) & 
               is.finite(Cd) & is.finite(.data[[water_col]])) %>%
      select(Cd, all_of(water_col))
    
    if(nrow(pairs) >= 5) {  # Reduced threshold for averaged data
      cor_result <- cor.test(pairs$Cd, pairs[[water_col]], method = "pearson")
      
      cat("  ", water_col, ": r =", round(cor_result$estimate, 3), 
          "(n =", nrow(pairs), ")\n")
      
      # Store results
      correlation_results <- rbind(correlation_results, 
                                   data.frame(
                                     Compartment = comp,
                                     Water_Fraction = water_col,
                                     Correlation = cor_result$estimate,
                                     P_value = cor_result$p.value,
                                     N = nrow(pairs)
                                   ))
    }
  }
}

# 4. Calculate bioconcentration factors (BCF)
cat("\n=== BIOCONCENTRATION FACTORS (Site/Date Averaged) ===\n")

bcf_results <- averaged_data %>%
  filter(SAMPLE_DESCRIPTOR %in% main_compartments & 
           !is.na(UF_Cd111) & is.finite(UF_Cd111) & UF_Cd111 > 0) %>%
  mutate(BCF = Cd / UF_Cd111) %>%
  group_by(SAMPLE_DESCRIPTOR) %>%
  summarise(
    n_site_date = n(),
    total_samples = sum(n_samples, na.rm = TRUE),
    mean_BCF = mean(BCF, na.rm = TRUE),
    median_BCF = median(BCF, na.rm = TRUE),
    sd_BCF = sd(BCF, na.rm = TRUE),
    min_BCF = min(BCF, na.rm = TRUE),
    max_BCF = max(BCF, na.rm = TRUE),
    .groups = 'drop'
  )
print(bcf_results)

# 5. Identify strongest correlations (r > 0.3)
cat("\n=== STRONGEST CORRELATIONS (r > 0.3) ===\n")
strong_correlations <- correlation_results %>%
  filter(abs(Correlation) > 0.3) %>%
  arrange(desc(abs(Correlation)))

if(nrow(strong_correlations) > 0) {
  print(strong_correlations)
} else {
  cat("No correlations > 0.3 found\n")
}

# 6. Summary of water cadmium hierarchy
cat("\n=== WATER CADMIUM CONCENTRATION HIERARCHY (Site/Date Averaged) ===\n")
water_means <- averaged_data %>%
  summarise(
    UF_mean = mean(UF_Cd111, na.rm = TRUE),
    F_mean = mean(F_Cd111, na.rm = TRUE),
    W_mean = mean(W_Cd111, na.rm = TRUE)
  )

cat("Ultrafiltered (UF):", round(water_means$UF_mean, 5), "\n")
cat("Filtered (F):", round(water_means$F_mean, 5), "\n")
cat("Whole water (W):", round(water_means$W_mean, 5), "\n")
cat("Hierarchy: Whole water > Filtered > Ultrafiltered\n")

# 7. Create summary data frame for future table creation
summary_df <- data.frame(
  Metric = c(
    "Total individual samples",
    "Site/date combinations",
    "EPIL site/date combinations", 
    "EPIP site/date combinations",
    "FILA site/date combinations",
    "EPIL mean Cd (mg/g)",
    "EPIP mean Cd (mg/g)", 
    "FILA mean Cd (mg/g)",
    "EPIL strongest correlation",
    "EPIP strongest correlation",
    "FILA strongest correlation",
    "Correlations > |0.3|",
    "Water Cd hierarchy",
    "Tissue Cd hierarchy",
    "Bioavailability pattern"
  ),
  Value = c(
    nrow(clean_data),
    nrow(averaged_data),
    sum(averaged_data$SAMPLE_DESCRIPTOR == "EPIL"),
    sum(averaged_data$SAMPLE_DESCRIPTOR == "EPIP"), 
    sum(averaged_data$SAMPLE_DESCRIPTOR == "FILA"),
    round(tissue_stats$mean_Cd[tissue_stats$SAMPLE_DESCRIPTOR == "EPIL"], 5),
    round(tissue_stats$mean_Cd[tissue_stats$SAMPLE_DESCRIPTOR == "EPIP"], 5),
    round(tissue_stats$mean_Cd[tissue_stats$SAMPLE_DESCRIPTOR == "FILA"], 5),
    paste0("r = ", round(max(abs(correlation_results$Correlation[correlation_results$Compartment == "EPIL"]), na.rm = TRUE), 3)),
    paste0("r = ", round(max(abs(correlation_results$Correlation[correlation_results$Compartment == "EPIP"]), na.rm = TRUE), 3)),
    paste0("r = ", round(max(abs(correlation_results$Correlation[correlation_results$Compartment == "FILA"]), na.rm = TRUE), 3)),
    sum(abs(correlation_results$Correlation) > 0.3, na.rm = TRUE),
    "W > F > UF",
    paste(tissue_stats$SAMPLE_DESCRIPTOR[order(-tissue_stats$mean_Cd)], collapse = " > "),
    "Dissolved >> Particulate"
  )
)

# Display the summary table
cat("\n=== SUMMARY TABLE FOR HYPOTHESIS 1 ===\n")
print(summary_df, row.names = FALSE)

# Create detailed correlation results table
correlation_table <- correlation_results %>%
  arrange(Compartment, desc(abs(Correlation))) %>%
  mutate(
    Correlation = round(Correlation, 3),
    P_value = round(P_value, 4),
    Significant = ifelse(abs(Correlation) > 0.3, "Yes", "No")
  ) %>%
  select(Compartment, Water_Fraction, Correlation, N, Significant)

cat("\n=== CORRELATION RESULTS TABLE ===\n")
print(correlation_table, row.names = FALSE)
