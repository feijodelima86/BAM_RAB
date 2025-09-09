# Step 1: Data Quality & Basic Transport Patterns for Cadmium Mass Balance
# Goal: Understand the physical transport framework before isolating biological effects

library(tidyverse)
library(lubridate)
library(viridis)
library(patchwork)
library(scales)
library(ggplot2)

# ============================================================================
# PART 1: DATA PREPARATION AND QUALITY CHECK
# ============================================================================

# Load your data (adjust path as needed)
df <- read.csv("2_incremental/AI_Subsets/integrated_parsed.csv")

# Create ordered site factor based on distance
site_order <- df %>%
  select(Site, Distance_km) %>%
  distinct() %>%
  arrange(Distance_km) %>%
  pull(Site)

df <- df %>%
  mutate(Site = factor(Site, levels = site_order))

# Quick data summary
cat("=== DATA SUMMARY ===\n")
cat("Date range:", as.character(min(df$Date)), "to", as.character(max(df$Date)), "\n")
cat("Number of sites:", n_distinct(df$Site), "\n")
cat("Number of sampling dates:", n_distinct(df$Date), "\n")
cat("Total observations:", nrow(df), "\n")
cat("Sites in order:", paste(site_order, collapse = ", "), "\n")
cat("Distance range:", min(df$Distance_km), "to", max(df$Distance_km), "km\n\n")

# Check for data gaps
cat("=== DATA QUALITY CHECK ===\n")
gap_summary <- df %>%
  summarise(
    mean_nutrient_gap = mean(nutrient_days_diff, na.rm = TRUE),
    max_nutrient_gap = max(nutrient_days_diff, na.rm = TRUE),
    mean_doc_gap = mean(doc_days_diff, na.rm = TRUE),
    max_doc_gap = max(doc_days_diff, na.rm = TRUE),
    mean_alkalinity_gap = mean(alkalinity_days_diff, na.rm = TRUE),
    max_alkalinity_gap = max(alkalinity_days_diff, na.rm = TRUE)
  )
print(gap_summary)

# ============================================================================
# PART 2: DISCHARGE PATTERNS - IDENTIFYING GAINS AND LOSSES
# ============================================================================

cat("\n=== DISCHARGE ANALYSIS ===\n")

# Calculate discharge changes between consecutive sites
discharge_changes <- df %>%
  arrange(Date, Distance_km) %>%
  group_by(Date) %>%
  mutate(
    discharge_change = discharge_m3s - lag(discharge_m3s),
    discharge_pct_change = 100 * discharge_change / lag(discharge_m3s),
    upstream_site = lag(Site),
    upstream_distance = lag(Distance_km),
    distance_between = Distance_km - upstream_distance,
    reach = paste(upstream_site, "->", Site),
    # Flag tributary inputs vs abstractions
    flow_type = case_when(
      discharge_change > 0.1 ~ "Tributary Input",
      discharge_change < -0.1 ~ "Abstraction/Loss",
      TRUE ~ "Minimal Change"
    )
  ) %>%
  filter(!is.na(discharge_change)) %>%
  ungroup()

# Summary statistics of discharge changes by reach
discharge_summary <- discharge_changes %>%
  group_by(reach, upstream_distance, Distance_km) %>%
  summarise(
    reach_length_km = first(distance_between),
    mean_discharge_change = mean(discharge_change, na.rm = TRUE),
    median_discharge_change = median(discharge_change, na.rm = TRUE),
    sd_discharge_change = sd(discharge_change, na.rm = TRUE),
    mean_pct_change = mean(discharge_pct_change, na.rm = TRUE),
    n_observations = n(),
    n_gains = sum(discharge_change > 0.1),
    n_losses = sum(discharge_change < -0.1),
    dominant_flow = ifelse(n_gains > n_losses, "Gaining", 
                           ifelse(n_losses > n_gains, "Losing", "Variable")),
    .groups = 'drop'
  ) %>%
  arrange(upstream_distance)

cat("\nDischarge changes by reach:\n")
print(discharge_summary)

# Identify reaches with consistent water losses (potential abstraction)
abstraction_reaches <- discharge_summary %>%
  filter(dominant_flow == "Losing" & mean_discharge_change < -0.1) %>%
  pull(reach)

if(length(abstraction_reaches) > 0) {
  cat("\nReaches with likely water abstraction:\n")
  print(abstraction_reaches)
} else {
  cat("\nNo reaches with consistent water abstraction identified.\n")
}

# Identify reaches with consistent tributary inputs
tributary_reaches <- discharge_summary %>%
  filter(dominant_flow == "Gaining" & mean_discharge_change > 0.1) %>%
  pull(reach)

if(length(tributary_reaches) > 0) {
  cat("\nReaches with likely tributary inputs:\n")
  print(tributary_reaches)
}

# Plot 1: Discharge patterns along river
p1 <- ggplot(df, aes(x = Distance_km, y = discharge_m3s)) +
  geom_point(aes(color = as.factor(month(Date))), alpha = 0.6, size = 3) +
  geom_line(aes(group = Date, color = as.factor(month(Date))), alpha = 0.3) +
  scale_color_viridis_d(name = "Month") +
  labs(
    title = "Discharge Patterns Along River",
    subtitle = "Lines connect same sampling date",
    x = "Distance from headwater (km)",
    y = expression("Discharge (m"^3*"/s)")
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Plot 2: Discharge changes between sites
p2 <- discharge_changes %>%
  ggplot(aes(x = reach, y = discharge_change)) +
  geom_boxplot(aes(fill = flow_type), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_viridis_d() +
  labs(
    title = "Discharge Changes Between Sites",
    subtitle = "Positive = tributary input, Negative = abstraction/loss",
    x = "River Reach",
    y = expression("Discharge Change (m"^3*"/s)")
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================================
# PART 3: IDENTIFY MINIMUM CD CONCENTRATION FOR TRIBUTARY ASSUMPTION
# ============================================================================

cat("\n=== CADMIUM CONCENTRATION ANALYSIS ===\n")

# Find minimum Cd concentrations across entire dataset
min_concentrations <- df %>%
  summarise(
    min_W_Cd = min(W_Cd, na.rm = TRUE),
    min_UF_Cd = min(UF_Cd, na.rm = TRUE),
    min_F_Cd = min(F_Cd, na.rm = TRUE),
    percentile_5_W_Cd = quantile(W_Cd, 0.05, na.rm = TRUE),
    percentile_10_W_Cd = quantile(W_Cd, 0.10, na.rm = TRUE)
  )

cat("Minimum and low percentile Cd concentrations:\n")
print(min_concentrations)

# Find where and when minimum concentrations occur
min_W_Cd_info <- df %>%
  filter(W_Cd == min(W_Cd, na.rm = TRUE)) %>%
  select(Site, Date, Distance_km, W_Cd, discharge_m3s, Position)

cat("\nMinimum W_Cd occurred at:\n")
print(min_W_Cd_info)

# Distribution of Cd concentrations by site
cd_by_site <- df %>%
  group_by(Site, Distance_km, Position) %>%
  summarise(
    n_samples = n(),
    mean_W_Cd = mean(W_Cd, na.rm = TRUE),
    median_W_Cd = median(W_Cd, na.rm = TRUE),
    min_W_Cd = min(W_Cd, na.rm = TRUE),
    max_W_Cd = max(W_Cd, na.rm = TRUE),
    cv_W_Cd = sd(W_Cd, na.rm = TRUE) / mean(W_Cd, na.rm = TRUE),
    mean_discharge = mean(discharge_m3s, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Distance_km)

cat("\nCadmium statistics by site:\n")
print(cd_by_site)

# Identify potential clean tributary sites (consistently low Cd)
clean_sites <- cd_by_site %>%
  filter(mean_W_Cd < quantile(cd_by_site$mean_W_Cd, 0.25)) %>%
  pull(Site)

cat("\nPotential clean water sites (bottom 25% mean Cd):\n")
print(as.character(clean_sites))

# Check if headwater is contaminated as expected
headwater_stats <- cd_by_site %>%
  filter(Distance_km == 0)

cat("\nHeadwater (Distance = 0) statistics:\n")
print(headwater_stats)

# ============================================================================
# PART 4: SPATIAL-TEMPORAL PATTERNS IN CD CONCENTRATIONS AND FRACTIONS
# ============================================================================

# Plot 3: Spatial pattern of Cd concentrations
p3 <- ggplot(df, aes(x = Distance_km, y = W_Cd)) +
  geom_point(aes(color = discharge_m3s), size = 3, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2, color = "red") +
  scale_color_viridis_c(name = expression("Discharge (m"^3*"/s)")) +
  labs(
    title = "Spatial Pattern of Total Cadmium",
    subtitle = "LOESS smooth shows general trend",
    x = "Distance from headwater (km)",
    y = "W_Cd concentration"
  ) +
  theme_minimal()

# Plot 4: Temporal patterns
p4 <- df %>%
  group_by(Date) %>%
  summarise(
    mean_W_Cd = mean(W_Cd, na.rm = TRUE),
    se_W_Cd = sd(W_Cd, na.rm = TRUE) / sqrt(n()),
    mean_discharge = mean(discharge_m3s, na.rm = TRUE),
    mean_temp = mean(daily_temp, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = Date)) +
  geom_ribbon(aes(ymin = mean_W_Cd - se_W_Cd, 
                  ymax = mean_W_Cd + se_W_Cd), 
              alpha = 0.3, fill = "blue") +
  geom_line(aes(y = mean_W_Cd), color = "blue", size = 1) +
  geom_point(aes(y = mean_W_Cd, size = mean_discharge), 
             color = "blue", alpha = 0.7) +
  scale_size_continuous(name = expression("Mean Q (m"^3*"/s)")) +
  labs(
    title = "Temporal Pattern of Cadmium",
    subtitle = "System-wide average ± SE",
    x = "Date",
    y = "Mean W_Cd concentration"
  ) +
  theme_minimal()

# Plot 5: Particle fraction patterns
p5 <- df %>%
  select(Site, Date, Distance_km, 
         Dissolved_pct, Colloidal_pct, Suspended_pct) %>%
  pivot_longer(cols = c(Dissolved_pct, Colloidal_pct, Suspended_pct),
               names_to = "Fraction",
               values_to = "Percentage") %>%
  mutate(Fraction = str_remove(Fraction, "_pct")) %>%
  group_by(Distance_km, Fraction) %>%
  summarise(mean_pct = mean(Percentage, na.rm = TRUE), .groups = 'drop') %>%
  ggplot(aes(x = Distance_km, y = mean_pct, fill = Fraction)) +
  geom_area(alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(
    title = "Average Cadmium Partitioning Along River",
    x = "Distance from headwater (km)",
    y = "Mean Percentage (%)",
    fill = "Fraction"
  ) +
  theme_minimal()

# Plot 6: Discharge-Concentration Relationship
p6 <- ggplot(df, aes(x = discharge_m3s, y = W_Cd)) +
  geom_point(aes(color = Distance_km), alpha = 0.6, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_viridis_c(name = "Distance (km)") +
  labs(
    title = "Discharge-Concentration Relationship",
    subtitle = "Log-log scale; dashed line = linear fit",
    x = expression("Discharge (m"^3*"/s) [log scale]"),
    y = "W_Cd concentration [log scale]"
  ) +
  annotation_logticks() +
  theme_minimal()

# Calculate dilution power law exponent
lm_result <- lm(log10(W_Cd) ~ log10(discharge_m3s), data = df)
slope <- coef(lm_result)[2]
cat("\n=== DILUTION ANALYSIS ===\n")
cat("Power law exponent (C = aQ^b): b =", round(slope, 3), "\n")
cat("If b = -1, perfect dilution; b > -1, additional sources; b < -1, removal processes\n")

# ============================================================================
# PART 5: CALCULATE CD LOADS AND IDENTIFY TRANSPORT ANOMALIES
# ============================================================================

# Calculate Cd loads and changes between sites
load_changes <- df %>%
  arrange(Date, Distance_km) %>%
  group_by(Date) %>%
  mutate(
    # Load changes
    W_Cd_load_change = W_Cd_load - lag(W_Cd_load),
    UF_Cd_load_change = UF_Cd_load - lag(UF_Cd_load),
    F_Cd_load_change = F_Cd_load - lag(F_Cd_load),
    # Expected from simple dilution (assuming tributary has min Cd)
    tributary_water = discharge_m3s - lag(discharge_m3s),
    expected_trib_Cd_load = tributary_water * min(W_Cd, na.rm = TRUE),
    expected_total_load = lag(W_Cd_load) + expected_trib_Cd_load,
    # Deviation from expectation
    load_residual = W_Cd_load - expected_total_load,
    residual_pct = 100 * load_residual / expected_total_load,
    # Site info
    upstream_site = lag(Site),
    reach = paste(upstream_site, "->", Site)
  ) %>%
  filter(!is.na(W_Cd_load_change)) %>%
  ungroup()

# Summary of load changes by reach
load_summary <- load_changes %>%
  group_by(reach) %>%
  summarise(
    n_obs = n(),
    mean_load_change = mean(W_Cd_load_change, na.rm = TRUE),
    sd_load_change = sd(W_Cd_load_change, na.rm = TRUE),
    mean_residual = mean(load_residual, na.rm = TRUE),
    mean_residual_pct = mean(residual_pct, na.rm = TRUE),
    n_sources = sum(load_residual > 0, na.rm = TRUE),
    n_sinks = sum(load_residual < 0, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(abs(mean_residual)))

cat("\n=== CADMIUM LOAD ANALYSIS ===\n")
cat("Reaches with largest deviations from simple transport:\n")
print(load_summary)

# Identify hotspots (large positive residuals = unexplained sources)
hotspots <- load_summary %>%
  filter(mean_residual_pct > 20) %>%
  pull(reach)

if(length(hotspots) > 0) {
  cat("\nPotential Cd source hotspots (>20% above expected):\n")
  print(hotspots)
}

# Plot 7: Load residuals by reach
p7 <- load_changes %>%
  ggplot(aes(x = reach, y = residual_pct)) +
  geom_boxplot(aes(fill = residual_pct > 0), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(
    values = c("TRUE" = "darkgreen", "FALSE" = "darkred"),
    labels = c("TRUE" = "Source", "FALSE" = "Sink"),
    name = "Process"
  ) +
  labs(
    title = "Deviations from Simple Transport Model",
    subtitle = "% difference from expected based on dilution alone",
    x = "River Reach",
    y = "Residual (%)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================================
# PART 6: ENVIRONMENTAL CONTEXT
# ============================================================================

# Create summary of environmental conditions by site
env_summary <- df %>%
  group_by(Site, Distance_km, Position) %>%
  summarise(
    mean_pH = mean(pH, na.rm = TRUE),
    range_pH = max(pH, na.rm = TRUE) - min(pH, na.rm = TRUE),
    mean_DOC = mean(DOC, na.rm = TRUE),
    mean_temp = mean(daily_temp, na.rm = TRUE),
    mean_DO = mean(daily_DO, na.rm = TRUE),
    mean_alkalinity = mean(Alkalinity, na.rm = TRUE),
    mean_NO3N = mean(NO3N, na.rm = TRUE),
    mean_NH4N = mean(NH4N, na.rm = TRUE),
    mean_SRP = mean(SRP, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Distance_km)

# Plot 8: Key environmental gradients
env_plot_data <- env_summary %>%
  select(Distance_km, mean_pH, mean_DOC, mean_temp, mean_alkalinity) %>%
  pivot_longer(cols = starts_with("mean_"),
               names_to = "Parameter",
               values_to = "Value") %>%
  mutate(Parameter = case_when(
    Parameter == "mean_pH" ~ "pH",
    Parameter == "mean_DOC" ~ "DOC (mg/L)",
    Parameter == "mean_temp" ~ "Temperature (°C)",
    Parameter == "mean_alkalinity" ~ "Alkalinity"
  ))

p8 <- ggplot(env_plot_data, aes(x = Distance_km, y = Value)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 3) +
  facet_wrap(~Parameter, scales = "free_y", ncol = 2) +
  labs(
    title = "Environmental Gradients Along River",
    x = "Distance from headwater (km)",
    y = "Mean Value"
  ) +
  theme_minimal()

# ============================================================================
# COMBINE AND DISPLAY PLOTS
# ============================================================================

# Create combined plot layout
combined_plot <- (p1 + p2) / (p3 + p6) / (p4 + p5) / (p7 + p8)

print(combined_plot)

# ============================================================================
# CREATE SUMMARY OUTPUT FOR STEP 2
# ============================================================================

# Compile key findings for mass balance calculations
transport_framework <- list(
  # Minimum Cd for tributary assumption
  tributary_Cd = min_concentrations$min_W_Cd,
  tributary_Cd_5pct = min_concentrations$percentile_5_W_Cd,
  
  # Flow characteristics
  abstraction_reaches = abstraction_reaches,
  tributary_reaches = tributary_reaches,
  discharge_summary = discharge_summary,
  
  # Cd patterns
  clean_sites = clean_sites,
  contaminated_headwater = headwater_stats,
  site_Cd_stats = cd_by_site,
  
  # Transport anomalies
  hotspot_reaches = hotspots,
  load_residuals = load_summary,
  
  # Environmental context
  env_gradients = env_summary,
  
  # Dilution characteristics
  dilution_exponent = slope
)

# Save prepared datasets for Step 2
step2_data <- list(
  main_data = df,
  discharge_changes = discharge_changes,
  load_changes = load_changes,
  transport_framework = transport_framework
)

# Save for use in next step
# saveRDS(step2_data, "step1_output_for_step2.rds")

cat("\n=== STEP 1 ANALYSIS COMPLETE ===\n")
cat("Minimum Cd for tributary assumption:", round(min_concentrations$min_W_Cd, 6), "\n")
cat("Alternative (5th percentile):", round(min_concentrations$percentile_5_W_Cd, 6), "\n")
cat("Dilution power law exponent:", round(slope, 3), "\n")
cat("Number of hotspot reaches identified:", length(hotspots), "\n")
cat("\nData prepared for Step 2: Mass Balance Calculations\n")
cat("Access results in 'transport_framework' object\n")