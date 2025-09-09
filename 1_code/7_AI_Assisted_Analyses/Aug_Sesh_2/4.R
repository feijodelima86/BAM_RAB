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

# Calculate biological deltas for each compartment
compartment_deltas <- compartment_data %>%
  group_by(Site, Distance_km, Compartment) %>%
  arrange(Date) %>%
  mutate(
    Prev_date = lag(Date),
    Delta_days = as.numeric(Date - Prev_date),
    
    # Biomass changes
    Prev_biomass = lag(Mean_biomass_g_m2),
    Delta_biomass_g_m2 = Mean_biomass_g_m2 - Prev_biomass,
    Biomass_change_rate_per_day = ifelse(Delta_days > 0, Delta_biomass_g_m2 / Delta_days, NA),
    
    # Cd standing crop changes  
    Prev_Cd_crop = lag(Mean_Cd_standing_crop_mg_m2),
    Delta_Cd_crop_mg_m2 = Mean_Cd_standing_crop_mg_m2 - Prev_Cd_crop,
    Cd_change_rate_per_day = ifelse(Delta_days > 0, Delta_Cd_crop_mg_m2 / Delta_days, NA)
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
                              Cd_change_rate_per_day)
  ) %>%
  select(-n_biomass_points, -n_cd_points) %>%
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
    Total_weight = sum(Weight, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  rename(Site = Target_site) %>%
  left_join(site_info, by = "Site")

# =================== COMPREHENSIVE CORRELATION ANALYSIS ===================

# Merge weighted biological data with water chemistry
full_spectrum_data <- weighted_bio_influence %>%
  left_join(water_interpolated, by = c("Date", "Site", "Distance_km")) %>%
  filter(!is.na(Weighted_biomass_change), !is.na(Dissolved_change_interp))

# Create correlation matrix: Each compartment vs each water fraction
correlation_matrix <- full_spectrum_data %>%
  group_by(Site, Distance_km, Compartment) %>%
  summarise(
    n_observations = n(),
    
    # Biomass change correlations
    Bio_vs_Dissolved = cor(Weighted_biomass_change, Dissolved_change_interp, use = "complete.obs"),
    Bio_vs_Colloidal = cor(Weighted_biomass_change, Colloidal_change_interp, use = "complete.obs"),
    Bio_vs_Suspended = cor(Weighted_biomass_change, Suspended_change_interp, use = "complete.obs"),
    
    # Cd change correlations  
    Cd_vs_Dissolved = cor(Weighted_Cd_change, Dissolved_change_interp, use = "complete.obs"),
    Cd_vs_Colloidal = cor(Weighted_Cd_change, Colloidal_change_interp, use = "complete.obs"),
    Cd_vs_Suspended = cor(Weighted_Cd_change, Suspended_change_interp, use = "complete.obs"),
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

cat("=== FULL SPECTRUM BIOLOGICAL-CHEMICAL ANALYSIS ===\n\n")

cat("=== CORRELATION MATRIX BY SITE AND COMPARTMENT ===\n")
print(correlation_matrix)

cat("\n=== TEMPORAL PATTERNS BY MONTH/SEASON ===\n") 
print(monthly_patterns)

cat("\n=== SPATIAL HOTSPOT RANKING ===\n")
print(hotspot_analysis)

cat("\n=== SUMMARY STATISTICS ===\n")
summary_stats <- full_spectrum_data %>%
  group_by(Compartment) %>%
  summarise(
    n_total_obs = n(),
    Mean_weighted_bio_change = mean(Weighted_biomass_change, na.rm = TRUE),
    Mean_dissolved_delta = mean(Dissolved_change_interp, na.rm = TRUE),
    Mean_suspended_delta = mean(Suspended_change_interp, na.rm = TRUE),
    
    Overall_bio_dissolved_cor = cor(Weighted_biomass_change, Dissolved_change_interp, use = "complete.obs"),
    Overall_bio_suspended_cor = cor(Weighted_biomass_change, Suspended_change_interp, use = "complete.obs"),
    .groups = 'drop'
  )
print(summary_stats)

# Export all results
write.csv(full_spectrum_data, "full_spectrum_analysis_data.csv", row.names = FALSE)
write.csv(correlation_matrix, "site_compartment_correlations.csv", row.names = FALSE)
write.csv(monthly_patterns, "temporal_hotspot_patterns.csv", row.names = FALSE)
write.csv(hotspot_analysis, "spatial_biological_hotspots.csv", row.names = FALSE)
write.csv(bio_interpolated, "interpolated_biological_data.csv", row.names = FALSE)
write.csv(water_interpolated, "interpolated_water_chemistry.csv", row.names = FALSE)

cat("\n=== INTERPRETATION GUIDE ===\n")
cat("HOTSPOT SCORES: Higher = stronger biological-chemical coupling\n")
cat("CORRELATIONS: Positive = synchronized changes, Negative = opposing changes\n")
cat("DISTANCE WEIGHTING: Accounts for cumulative upstream biological effects\n")
cat("COMPARTMENTS: EPIL=epilithic, EPIP=epiphytic, FILA=filamentous algae\n")
cat("FRACTIONS: Dissolved=bioavailable, Colloidal=intermediate, Suspended=particulate\n")