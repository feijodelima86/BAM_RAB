# Full Spectrum Analysis: Biological Compartments vs Water Column Cd Fractions
# Distance-weighted, temporally resolved, with interpolation

library(dplyr)
library(tidyr)
library(lubridate)

# Assume previous data processing steps are completed
# (algal_data, water_data, site_info loaded and processed)

# =================== BIOLOGICAL DATA BY COMPARTMENT ===================

# Get individual compartment data with temporal changes
compartment_data <- algal_clean %>%
  group_by(Date, Site, Distance_km, Compartment) %>%
  summarise(
    Mean_biomass_g_m2 = mean(Biomass_g_m2, na.rm = TRUE),
    Mean_Cd_standing_crop_mg_m2 = mean(Cd_standing_crop_mg_m2, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Site, Compartment, Date)

# Calculate biological deltas for each compartment - SEPARATED EFFECTS
compartment_deltas <- compartment_data %>%
  group_by(Site, Distance_km, Compartment) %>%
  arrange(Date) %>%
  mutate(
    Prev_date = lag(Date),
    Delta_days = as.numeric(Date - Prev_date),
    
    # Previous values
    Prev_biomass = lag(Mean_biomass_g_m2),
    Prev_Cd_crop = lag(Mean_Cd_standing_crop_mg_m2),
    
    # Calculate tissue Cd concentration (mg Cd per g biomass)
    Cd_concentration_mg_g = ifelse(Mean_biomass_g_m2 > 0, Mean_Cd_standing_crop_mg_m2 / Mean_biomass_g_m2, NA),
    Prev_Cd_concentration = ifelse(Prev_biomass > 0, Prev_Cd_crop / Prev_biomass, NA),
    
    # SEPARATED EFFECTS
    # 1. Pure biomass effect (growth/death at constant Cd concentration)
    Delta_biomass_g_m2 = Mean_biomass_g_m2 - Prev_biomass,
    Biomass_change_rate_per_day = ifelse(Delta_days > 0, Delta_biomass_g_m2 / Delta_days, NA),
    Pure_biomass_Cd_effect = ifelse(!is.na(Prev_Cd_concentration), 
                                    Delta_biomass_g_m2 * Prev_Cd_concentration, NA),
    Pure_biomass_Cd_rate = ifelse(Delta_days > 0, Pure_biomass_Cd_effect / Delta_days, NA),
    
    # 2. Pure tissue concentration effect (uptake/release at constant biomass)
    Delta_Cd_concentration = Cd_concentration_mg_g - Prev_Cd_concentration,
    Cd_concentration_change_rate = ifelse(Delta_days > 0, Delta_Cd_concentration / Delta_days, NA),
    Pure_concentration_Cd_effect = ifelse(!is.na(Prev_biomass), 
                                          Delta_Cd_concentration * Prev_biomass, NA),
    Pure_concentration_Cd_rate = ifelse(Delta_days > 0, Pure_concentration_Cd_effect / Delta_days, NA),
    
    # 3. Combined effect (total standing crop change)
    Delta_Cd_crop_mg_m2 = Mean_Cd_standing_crop_mg_m2 - Prev_Cd_crop,
    Cd_change_rate_per_day = ifelse(Delta_days > 0, Delta_Cd_crop_mg_m2 / Delta_days, NA),
    
    # Verification: Pure biomass + Pure concentration should = Combined
    Combined_check = Pure_biomass_Cd_effect + Pure_concentration_Cd_effect,
    
    # Classify dominant mechanism
    Dominant_mechanism = case_when(
      abs(Pure_biomass_Cd_effect) > abs(Pure_concentration_Cd_effect) ~ "Growth/Death Driven",
      abs(Pure_concentration_Cd_effect) > abs(Pure_biomass_Cd_effect) ~ "Uptake/Release Driven", 
      TRUE ~ "Mixed/Equal"
    )
  ) %>%
  filter(!is.na(Prev_date)) %>%
  ungroup()

# =================== WATER CHEMISTRY DELTAS ===================

# Calculate water chemistry deltas between sampling dates
water_deltas <- water_clean %>%
  group_by(Site, Distance_km) %>%
  arrange(Date) %>%
  mutate(
    Prev_date = lag(Date),
    Delta_days_water = as.numeric(Date - Prev_date),
    
    # Previous values
    Prev_dissolved = lag(Dissolved_Cd_mg_L),
    Prev_colloidal = lag(Colloidal_Cd_mg_L), 
    Prev_suspended = lag(Suspended_Cd_mg_L),
    
    # Delta changes in each fraction
    Delta_dissolved_mg_L = Dissolved_Cd_mg_L - Prev_dissolved,
    Delta_colloidal_mg_L = Colloidal_Cd_mg_L - Prev_colloidal,
    Delta_suspended_mg_L = Suspended_Cd_mg_L - Prev_suspended,
    
    # Rate of change per day
    Dissolved_change_rate = ifelse(Delta_days_water > 0, Delta_dissolved_mg_L / Delta_days_water, NA),
    Colloidal_change_rate = ifelse(Delta_days_water > 0, Delta_colloidal_mg_L / Delta_days_water, NA),
    Suspended_change_rate = ifelse(Delta_days_water > 0, Delta_suspended_mg_L / Delta_days_water, NA)
  ) %>%
  filter(!is.na(Prev_date)) %>%
  ungroup()

# =================== DISTANCE WEIGHTING SYSTEM ===================

# Create distance weights - upstream sites influence downstream
# Weight decreases with distance, but all upstream sites contribute
all_sites <- unique(compartment_deltas$Site)
all_distances <- site_info$Distance_km[match(all_sites, site_info$Site)]
names(all_distances) <- all_sites

# For each site, calculate weights from all upstream sites
distance_weights <- expand.grid(
  Target_site = all_sites,
  Source_site = all_sites,
  stringsAsFactors = FALSE
) %>%
  mutate(
    Target_distance = all_distances[Target_site],
    Source_distance = all_distances[Source_site],
    Distance_diff = Target_distance - Source_distance,
    
    # Weight calculation: upstream sites (negative diff) get positive weights
    # Downstream sites (positive diff) get zero weight
    Weight = ifelse(Distance_diff >= 0 & Distance_diff <= 50, 
                    exp(-Distance_diff / 25),  # Exponential decay over 50km
                    ifelse(Distance_diff > 50, 0.1 * exp(-(Distance_diff-50) / 50), 0))
  ) %>%
  filter(Weight > 0.01)  # Keep only meaningful weights

# =================== CREATE COMPREHENSIVE DATA MATRIX ===================

# Get all unique dates from both datasets
all_dates <- sort(unique(c(compartment_deltas$Date, water_deltas$Date)))

# Create complete grid of Site x Date x Compartment
complete_bio_grid <- expand.grid(
  Date = all_dates,
  Site = all_sites,
  Compartment = unique(compartment_deltas$Compartment),
  stringsAsFactors = FALSE
) %>%
  left_join(site_info, by = "Site") %>%
  left_join(compartment_deltas, by = c("Date", "Site", "Distance_km", "Compartment"))

# Create complete grid for water chemistry
complete_water_grid <- expand.grid(
  Date = all_dates,
  Site = all_sites,
  stringsAsFactors = FALSE
) %>%
  left_join(site_info, by = "Site") %>%
  left_join(water_deltas, by = c("Date", "Site", "Distance_km"))

# =================== INTERPOLATION FOR MISSING DATA ===================

# Simple interpolation by site and compartment for biological data
bio_interpolated <- complete_bio_grid %>%
  group_by(Site, Compartment) %>%
  arrange(Date) %>%
  mutate(
    # Check if we have enough data points for interpolation
    n_biomass_points = sum(!is.na(Biomass_change_rate_per_day)),
    n_cd_points = sum(!is.na(Cd_change_rate_per_day)),
    n_pure_biomass_points = sum(!is.na(Pure_biomass_Cd_rate)),
    n_pure_conc_points = sum(!is.na(Pure_concentration_Cd_rate)),
    
    # Linear interpolation for key variables (only if >=2 points)
    Biomass_change_interp = ifelse(n_biomass_points >= 2,
                                   approx(x = which(!is.na(Biomass_change_rate_per_day)), 
                                          y = Biomass_change_rate_per_day[!is.na(Biomass_change_rate_per_day)], 
                                          xout = 1:n(), rule = 2)$y,
                                   Biomass_change_rate_per_day),
    
    Cd_change_interp = ifelse(n_cd_points >= 2,
                              approx(x = which(!is.na(Cd_change_rate_per_day)), 
                                     y = Cd_change_rate_per_day[!is.na(Cd_change_rate_per_day)], 
                                     xout = 1:n(), rule = 2)$y,
                              Cd_change_rate_per_day),
    
    # NEW: Interpolate separated effects
    Pure_biomass_Cd_interp = ifelse(n_pure_biomass_points >= 2,
                                    approx(x = which(!is.na(Pure_biomass_Cd_rate)), 
                                           y = Pure_biomass_Cd_rate[!is.na(Pure_biomass_Cd_rate)], 
                                           xout = 1:n(), rule = 2)$y,
                                    Pure_biomass_Cd_rate),
    
    Pure_concentration_Cd_interp = ifelse(n_pure_conc_points >= 2,
                                          approx(x = which(!is.na(Pure_concentration_Cd_rate)), 
                                                 y = Pure_concentration_Cd_rate[!is.na(Pure_concentration_Cd_rate)], 
                                                 xout = 1:n(), rule = 2)$y,
                                          Pure_concentration_Cd_rate)
  ) %>%
  select(-n_biomass_points, -n_cd_points, -n_pure_biomass_points, -n_pure_conc_points) %>%
  ungroup()

# Simple interpolation for water chemistry
water_interpolated <- complete_water_grid %>%
  group_by(Site) %>%
  arrange(Date) %>%
  mutate(
    # Check if we have enough data points for interpolation
    n_dissolved_points = sum(!is.na(Dissolved_change_rate)),
    n_colloidal_points = sum(!is.na(Colloidal_change_rate)),
    n_suspended_points = sum(!is.na(Suspended_change_rate)),
    
    # Linear interpolation for water deltas (only if >=2 points)
    Dissolved_change_interp = ifelse(n_dissolved_points >= 2,
                                     approx(x = which(!is.na(Dissolved_change_rate)), 
                                            y = Dissolved_change_rate[!is.na(Dissolved_change_rate)], 
                                            xout = 1:n(), rule = 2)$y,
                                     Dissolved_change_rate),
    
    Colloidal_change_interp = ifelse(n_colloidal_points >= 2,
                                     approx(x = which(!is.na(Colloidal_change_rate)), 
                                            y = Colloidal_change_rate[!is.na(Colloidal_change_rate)], 
                                            xout = 1:n(), rule = 2)$y,
                                     Colloidal_change_rate),
    
    Suspended_change_interp = ifelse(n_suspended_points >= 2,
                                     approx(x = which(!is.na(Suspended_change_rate)), 
                                            y = Suspended_change_rate[!is.na(Suspended_change_rate)], 
                                            xout = 1:n(), rule = 2)$y,
                                     Suspended_change_rate)
  ) %>%
  select(-n_dissolved_points, -n_colloidal_points, -n_suspended_points) %>%
  ungroup()

# =================== DISTANCE-WEIGHTED BIOLOGICAL INFLUENCE ===================

# Calculate distance-weighted biological influence for each site-date
weighted_bio_influence <- bio_interpolated %>%
  left_join(distance_weights, by = c("Site" = "Source_site")) %>%
  filter(!is.na(Target_site)) %>%
  group_by(Date, Target_site, Compartment) %>%
  summarise(
    Weighted_biomass_change = sum(Biomass_change_interp * Weight, na.rm = TRUE),
    Weighted_Cd_change = sum(Cd_change_interp * Weight, na.rm = TRUE),
    
    # NEW: Weighted separated effects
    Weighted_pure_biomass_effect = sum(Pure_biomass_Cd_interp * Weight, na.rm = TRUE),
    Weighted_pure_concentration_effect = sum(Pure_concentration_Cd_interp * Weight, na.rm = TRUE),
    
    Total_weight = sum(Weight, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  rename(Site = Target_site) %>%
  left_join(site_info, by = "Site")

# =================== BIOLOGICAL FLUX vs DOWNSTREAM LOADS COMPARISON ===================

# Calculate total biological flux by site and date (across all compartments)
total_bio_flux <- weighted_bio_influence %>%
  group_by(Date, Site, Distance_km, Width_m) %>%
  summarise(
    # Total effects across all compartments
    Total_biomass_effect_mg_m2_day = sum(Weighted_pure_biomass_effect, na.rm = TRUE),
    Total_concentration_effect_mg_m2_day = sum(Weighted_pure_concentration_effect, na.rm = TRUE),
    Total_combined_effect_mg_m2_day = sum(Weighted_Cd_change, na.rm = TRUE),
    
    # Convert to reach-scale flux (mg/day for entire reach width)
    Reach_biomass_flux_mg_day = Total_biomass_effect_mg_m2_day * Width_m * 1000,  # 1000m reach length
    Reach_concentration_flux_mg_day = Total_concentration_effect_mg_m2_day * Width_m * 1000,
    Reach_total_bio_flux_mg_day = Total_combined_effect_mg_m2_day * Width_m * 1000,
    
    # Classify biological activity
    Bio_mechanism = case_when(
      abs(Total_biomass_effect_mg_m2_day) > abs(Total_concentration_effect_mg_m2_day) ~ "Growth/Death Dominated",
      abs(Total_concentration_effect_mg_m2_day) > abs(Total_biomass_effect_mg_m2_day) ~ "Uptake/Release Dominated",
      TRUE ~ "Mixed Mechanisms"
    ),
    
    .groups = 'drop'
  )

# Merge with water column loads
bio_vs_loads <- total_bio_flux %>%
  left_join(water_clean, by = c("Date", "Site", "Distance_km")) %>%
  mutate(
    # Calculate daily loads for each particle size (mg/day)
    Dissolved_load_mg_day = ifelse(!is.na(Dissolved_Cd_mg_L) & !is.na(Discharge_m3_s),
                                   Dissolved_Cd_mg_L * Discharge_m3_s * 86400 * 1000, NA),
    Colloidal_load_mg_day = ifelse(!is.na(Colloidal_Cd_mg_L) & !is.na(Discharge_m3_s),
                                   Colloidal_Cd_mg_L * Discharge_m3_s * 86400 * 1000, NA),
    Suspended_load_mg_day = ifelse(!is.na(Suspended_Cd_mg_L) & !is.na(Discharge_m3_s),
                                   Suspended_Cd_mg_L * Discharge_m3_s * 86400 * 1000, NA),
    Total_water_load_mg_day = ifelse(!is.na(W_Cd_mg_L) & !is.na(Discharge_m3_s),
                                     W_Cd_mg_L * Discharge_m3_s * 86400 * 1000, NA),
    
    # Calculate biological flux as percentage of water column loads
    Bio_pct_dissolved_load = ifelse(Dissolved_load_mg_day > 0, 
                                    abs(Reach_total_bio_flux_mg_day) / Dissolved_load_mg_day * 100, NA),
    Bio_pct_colloidal_load = ifelse(Colloidal_load_mg_day > 0,
                                    abs(Reach_total_bio_flux_mg_day) / Colloidal_load_mg_day * 100, NA),
    Bio_pct_suspended_load = ifelse(Suspended_load_mg_day > 0,
                                    abs(Reach_total_bio_flux_mg_day) / Suspended_load_mg_day * 100, NA),
    Bio_pct_total_load = ifelse(Total_water_load_mg_day > 0,
                                abs(Reach_total_bio_flux_mg_day) / Total_water_load_mg_day * 100, NA),
    
    # Classify biological significance
    Bio_significance = case_when(
      Bio_pct_total_load > 10 ~ "Major Biological Effect (>10%)",
      Bio_pct_total_load > 1 ~ "Moderate Biological Effect (1-10%)",
      Bio_pct_total_load > 0.1 ~ "Minor Biological Effect (0.1-1%)",
      TRUE ~ "Negligible Biological Effect (<0.1%)"
    ),
    
    # Identify which particle size is most affected
    Most_affected_fraction = case_when(
      Bio_pct_dissolved_load == pmax(Bio_pct_dissolved_load, Bio_pct_colloidal_load, Bio_pct_suspended_load, na.rm = TRUE) ~ "Dissolved",
      Bio_pct_colloidal_load == pmax(Bio_pct_dissolved_load, Bio_pct_colloidal_load, Bio_pct_suspended_load, na.rm = TRUE) ~ "Colloidal",
      Bio_pct_suspended_load == pmax(Bio_pct_dissolved_load, Bio_pct_colloidal_load, Bio_pct_suspended_load, na.rm = TRUE) ~ "Suspended",
      TRUE ~ "Unknown"
    )
  )

# Merge weighted biological data with water chemistry
full_spectrum_data <- weighted_bio_influence %>%
  left_join(water_interpolated, by = c("Date", "Site", "Distance_km")) %>%
  filter(!is.na(Weighted_biomass_change), !is.na(Dissolved_change_interp))

# =================== COMPREHENSIVE CORRELATION ANALYSIS ===================

# Merge weighted biological data with water chemistry
full_spectrum_data <- weighted_bio_influence %>%
  left_join(water_interpolated, by = c("Date", "Site", "Distance_km")) %>%
  filter(!is.na(Weighted_biomass_change), !is.na(Dissolved_change_interp))

# Create correlation matrix: Each compartment vs each water fraction - WITH SEPARATED EFFECTS
correlation_matrix <- full_spectrum_data %>%
  group_by(Site, Distance_km, Compartment) %>%
  summarise(
    n_observations = n(),
    
    # Traditional biomass change correlations
    Bio_vs_Dissolved = cor(Weighted_biomass_change, Dissolved_change_interp, use = "complete.obs"),
    Bio_vs_Colloidal = cor(Weighted_biomass_change, Colloidal_change_interp, use = "complete.obs"),
    Bio_vs_Suspended = cor(Weighted_biomass_change, Suspended_change_interp, use = "complete.obs"),
    
    # Combined Cd change correlations  
    Cd_vs_Dissolved = cor(Weighted_Cd_change, Dissolved_change_interp, use = "complete.obs"),
    Cd_vs_Colloidal = cor(Weighted_Cd_change, Colloidal_change_interp, use = "complete.obs"),
    Cd_vs_Suspended = cor(Weighted_Cd_change, Suspended_change_interp, use = "complete.obs"),
    
    # NEW: Pure biomass effect correlations (growth/death at constant Cd concentration)
    PureBio_vs_Dissolved = cor(Weighted_pure_biomass_effect, Dissolved_change_interp, use = "complete.obs"),
    PureBio_vs_Colloidal = cor(Weighted_pure_biomass_effect, Colloidal_change_interp, use = "complete.obs"),
    PureBio_vs_Suspended = cor(Weighted_pure_biomass_effect, Suspended_change_interp, use = "complete.obs"),
    
    # NEW: Pure concentration effect correlations (uptake/release at constant biomass)
    PureConc_vs_Dissolved = cor(Weighted_pure_concentration_effect, Dissolved_change_interp, use = "complete.obs"),
    PureConc_vs_Colloidal = cor(Weighted_pure_concentration_effect, Colloidal_change_interp, use = "complete.obs"),
    PureConc_vs_Suspended = cor(Weighted_pure_concentration_effect, Suspended_change_interp, use = "complete.obs"),
    
    .groups = 'drop'
  ) %>%
  arrange(Distance_km, Compartment)

# =================== TEMPORAL HOTSPOT ANALYSIS ===================

# Monthly analysis
monthly_patterns <- full_spectrum_data %>%
  mutate(
    Month = month(Date),
    Season = case_when(
      Month %in% c(3,4,5) ~ "Spring",
      Month %in% c(6,7,8) ~ "Summer", 
      Month %in% c(9,10,11) ~ "Fall",
      TRUE ~ "Winter"
    )
  ) %>%
  group_by(Month, Season, Compartment) %>%
  summarise(
    n_obs = n(),
    Mean_bio_change = mean(Weighted_biomass_change, na.rm = TRUE),
    Mean_dissolved_change = mean(Dissolved_change_interp, na.rm = TRUE),
    Mean_suspended_change = mean(Suspended_change_interp, na.rm = TRUE),
    
    # Monthly correlations
    Monthly_bio_dissolved_cor = cor(Weighted_biomass_change, Dissolved_change_interp, use = "complete.obs"),
    Monthly_bio_suspended_cor = cor(Weighted_biomass_change, Suspended_change_interp, use = "complete.obs"),
    .groups = 'drop'
  ) %>%
  arrange(Month, Compartment)

# =================== SPATIAL HOTSPOT ANALYSIS ===================

# Distance-based analysis
spatial_patterns <- correlation_matrix %>%
  mutate(
    Distance_bin = case_when(
      Distance_km < 50 ~ "Upstream (0-50km)",
      Distance_km < 100 ~ "Middle (50-100km)",
      TRUE ~ "Downstream (100+km)"
    ),
    
    # Strongest correlation magnitude for each site-compartment
    Strongest_correlation = pmax(abs(Bio_vs_Dissolved), abs(Bio_vs_Colloidal), abs(Bio_vs_Suspended), na.rm = TRUE),
    
    Dominant_fraction = case_when(
      abs(Bio_vs_Dissolved) == Strongest_correlation ~ "Dissolved",
      abs(Bio_vs_Colloidal) == Strongest_correlation ~ "Colloidal", 
      abs(Bio_vs_Suspended) == Strongest_correlation ~ "Suspended",
      TRUE ~ "None"
    )
  )

# =================== HOTSPOT IDENTIFICATION ===================

# Identify top biological-chemical interaction hotspots
hotspot_analysis <- spatial_patterns %>%
  filter(!is.na(Strongest_correlation)) %>%
  mutate(
    Hotspot_score = Strongest_correlation * n_observations / 10,  # Weight by data quality
    Hotspot_category = case_when(
      Hotspot_score > 0.7 ~ "Strong Hotspot",
      Hotspot_score > 0.5 ~ "Moderate Hotspot", 
      Hotspot_score > 0.3 ~ "Weak Hotspot",
      TRUE ~ "No Signal"
    )
  ) %>%
  arrange(desc(Hotspot_score))

# =================== COMPREHENSIVE RESULTS ===================

cat("=== FULL SPECTRUM BIOLOGICAL-CHEMICAL ANALYSIS ===\n")
cat("=== WITH SEPARATED BIOMASS vs CONCENTRATION EFFECTS ===\n\n")

cat("=== BIOLOGICAL EFFECTS SUMMARY ===\n")
bio_effects_summary <- compartment_deltas %>%
  group_by(Compartment) %>%
  summarise(
    n_observations = n(),
    
    # Mean effects
    Mean_biomass_change = mean(Biomass_change_rate_per_day, na.rm = TRUE),
    Mean_concentration_change = mean(Cd_concentration_change_rate, na.rm = TRUE),
    Mean_combined_change = mean(Cd_change_rate_per_day, na.rm = TRUE),
    
    # Proportion of times each mechanism dominates
    Pct_growth_death_dominated = sum(Dominant_mechanism == "Growth/Death Driven", na.rm = TRUE) / n() * 100,
    Pct_uptake_release_dominated = sum(Dominant_mechanism == "Uptake/Release Driven", na.rm = TRUE) / n() * 100,
    Pct_mixed = sum(Dominant_mechanism == "Mixed/Equal", na.rm = TRUE) / n() * 100,
    
    .groups = 'drop'
  )
print(bio_effects_summary)

cat("\n=== BIOLOGICAL vs DOWNSTREAM LOADS COMPARISON ===\n")
load_comparison_summary <- bio_vs_loads %>%
  filter(!is.na(Bio_pct_total_load)) %>%
  group_by(Site, Distance_km) %>%
  summarise(
    n_observations = n(),
    
    # Average biological contribution to loads
    Mean_bio_pct_dissolved = mean(Bio_pct_dissolved_load, na.rm = TRUE),
    Mean_bio_pct_colloidal = mean(Bio_pct_colloidal_load, na.rm = TRUE),
    Mean_bio_pct_suspended = mean(Bio_pct_suspended_load, na.rm = TRUE),
    Mean_bio_pct_total = mean(Bio_pct_total_load, na.rm = TRUE),
    
    # Maximum biological contribution observed
    Max_bio_pct_total = max(Bio_pct_total_load, na.rm = TRUE),
    
    # Most common significance level
    Typical_significance = names(sort(table(Bio_significance), decreasing = TRUE))[1],
    
    # Most commonly affected fraction
    Most_affected = names(sort(table(Most_affected_fraction), decreasing = TRUE))[1],
    
    .groups = 'drop'
  ) %>%
  arrange(Distance_km)
print(load_comparison_summary)

cat("\n=== SEPARATED EFFECTS CORRELATION MATRIX ===\n")
print(correlation_matrix)

cat("\n=== MECHANISM ANALYSIS BY COMPARTMENT ===\n")
mechanism_analysis <- full_spectrum_data %>%
  group_by(Compartment) %>%
  summarise(
    n_observations = n(),
    
    # Compare pure effects
    Mean_pure_biomass_effect = mean(Weighted_pure_biomass_effect, na.rm = TRUE),
    Mean_pure_concentration_effect = mean(Weighted_pure_concentration_effect, na.rm = TRUE),
    
    # Correlations for each mechanism
    PureBio_dissolved_cor = cor(Weighted_pure_biomass_effect, Dissolved_change_interp, use = "complete.obs"),
    PureBio_suspended_cor = cor(Weighted_pure_biomass_effect, Suspended_change_interp, use = "complete.obs"),
    
    PureConc_dissolved_cor = cor(Weighted_pure_concentration_effect, Dissolved_change_interp, use = "complete.obs"),
    PureConc_suspended_cor = cor(Weighted_pure_concentration_effect, Suspended_change_interp, use = "complete.obs"),
    
    # Which mechanism shows stronger correlations?
    Stronger_mechanism = ifelse(
      abs(PureBio_dissolved_cor) + abs(PureBio_suspended_cor) > 
        abs(PureConc_dissolved_cor) + abs(PureConc_suspended_cor),
      "Growth/Death Effects", "Uptake/Release Effects"
    ),
    
    .groups = 'drop'
  )
print(mechanism_analysis)

cat("\n=== TEMPORAL HOTSPOTS WITH MECHANISMS ===\n")
temporal_mechanisms <- full_spectrum_data %>%
  mutate(
    Month = month(Date),
    Season = case_when(
      Month %in% c(3,4,5) ~ "Spring",
      Month %in% c(6,7,8) ~ "Summer", 
      Month %in% c(9,10,11) ~ "Fall",
      TRUE ~ "Winter"
    )
  ) %>%
  group_by(Season, Compartment) %>%
  summarise(
    n_obs = n(),
    
    # Seasonal averages
    Season_biomass_effect = mean(Weighted_pure_biomass_effect, na.rm = TRUE),
    Season_concentration_effect = mean(Weighted_pure_concentration_effect, na.rm = TRUE),
    
    # Seasonal correlations
    Season_bio_dissolved = cor(Weighted_pure_biomass_effect, Dissolved_change_interp, use = "complete.obs"),
    Season_conc_dissolved = cor(Weighted_pure_concentration_effect, Dissolved_change_interp, use = "complete.obs"),
    
    .groups = 'drop'
  ) %>%
  arrange(Season, Compartment)
print(temporal_mechanisms)

cat("\n=== SPATIAL HOTSPOTS FOR BIOLOGICAL SIGNIFICANCE ===\n")
spatial_significance <- bio_vs_loads %>%
  filter(!is.na(Bio_pct_total_load)) %>%
  group_by(Site, Distance_km) %>%
  summarise(
    n_observations = n(),
    
    # Biological significance metrics
    Mean_total_bio_flux = mean(abs(Reach_total_bio_flux_mg_day), na.rm = TRUE),
    Mean_pct_of_load = mean(Bio_pct_total_load, na.rm = TRUE),
    Max_pct_of_load = max(Bio_pct_total_load, na.rm = TRUE),
    
    # How often is biology significant?
    Pct_major_effect = sum(Bio_significance == "Major Biological Effect (>10%)", na.rm = TRUE) / n() * 100,
    Pct_moderate_effect = sum(Bio_significance == "Moderate Biological Effect (1-10%)", na.rm = TRUE) / n() * 100,
    
    # Dominant mechanism at this site
    Dominant_mechanism = names(sort(table(Bio_mechanism), decreasing = TRUE))[1],
    
    .groups = 'drop'
  ) %>%
  arrange(desc(Mean_pct_of_load))
print(spatial_significance)

# Export all results
write.csv(full_spectrum_data, "full_spectrum_analysis_data.csv", row.names = FALSE)
write.csv(correlation_matrix, "separated_effects_correlations.csv", row.names = FALSE)
write.csv(bio_vs_loads, "biological_vs_downstream_loads.csv", row.names = FALSE)
write.csv(bio_effects_summary, "biological_mechanisms_summary.csv", row.names = FALSE)
write.csv(mechanism_analysis, "compartment_mechanism_analysis.csv", row.names = FALSE)
write.csv(compartment_deltas, "detailed_biological_deltas.csv", row.names = FALSE)
write.csv(spatial_significance, "spatial_biological_significance.csv", row.names = FALSE)

cat("\n=== INTERPRETATION GUIDE ===\n")
cat("SEPARATED EFFECTS:\n")
cat("- Pure Biomass Effect: Cd change from growth/death at constant tissue concentration\n")
cat("- Pure Concentration Effect: Cd change from uptake/release at constant biomass\n")
cat("- Combined Effect: Total standing crop change (sum of above)\n\n")

cat("BIOLOGICAL SIGNIFICANCE:\n")
cat("- Major (>10%): Biology substantially affects downstream loads\n")
cat("- Moderate (1-10%): Biology noticeably affects loads\n")
cat("- Minor (0.1-1%): Biology slightly affects loads\n")
cat("- Negligible (<0.1%): Biology has minimal effect\n\n")

cat("MECHANISMS:\n")
cat("- Growth/Death Dominated: Changes driven by algae growing/dying\n")
cat("- Uptake/Release Dominated: Changes driven by Cd concentration in tissue\n")
cat("- Mixed: Both mechanisms equally important\n\n")

cat("KEY QUESTIONS ANSWERED:\n")
cat("1. Which compartments drive Cd cycling through growth vs uptake?\n")
cat("2. When/where do biological processes significantly affect downstream loads?\n")
cat("3. What proportion of Cd transport is biological vs physical?\n")
cat("4. Which particle sizes are most affected by biological activity?\n")