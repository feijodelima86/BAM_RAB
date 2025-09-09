# Biological Cadmium Flux Analysis - Delta Standing Crop Approach
# Calculate biological flux as change in Cd standing crop over time

library(dplyr)
library(tidyr)
library(lubridate)

# Read the datasets
algal_data <- read.csv("2_incremental/AI_Subsets/Cd_redux.csv")
water_data <- read.csv("2_incremental/AI_Subsets/integrated_parsed.csv")

# Site characteristics
site_info <- data.frame(
  Site = c("WS", "DL", "GR", "GC", "BG", "BN"),
  Width_m = c(11, 24, 19, 35, 27, 31),
  Distance_km = c(0, 44, 64.8, 89.2, 144.3, 167.8)
)

# Clean and prepare algal data
algal_clean <- algal_data %>%
  mutate(
    Date = as.Date(SAMPLING_DATE, format = "%m/%d/%Y"),
    Site = SITE,
    Compartment = SAMPLE_DESCRIPTOR,
    Biomass_g_m2 = OM.AREA.g.m2,
    Cd_mg_g = Cd,
    # Calculate Cd standing crop (mg/m²)
    Cd_standing_crop_mg_m2 = Biomass_g_m2 * Cd_mg_g
  ) %>%
  filter(!is.na(Date), !is.na(Site), !is.na(Compartment), 
         !is.na(Cd_standing_crop_mg_m2)) %>%
  left_join(site_info, by = "Site") %>%
  select(Date, Site, Distance_km, Compartment, Biomass_g_m2, 
         Cd_mg_g, Cd_standing_crop_mg_m2, Width_m) %>%
  arrange(Site, Compartment, Date)

# Calculate mean standing crop per site-compartment-date 
# (in case of multiple field replicates)
standing_crop_means <- algal_clean %>%
  group_by(Date, Site, Distance_km, Compartment, Width_m) %>%
  summarise(
    n_samples = n(),
    Mean_biomass_g_m2 = mean(Biomass_g_m2, na.rm = TRUE),
    Mean_Cd_mg_g = mean(Cd_mg_g, na.rm = TRUE),
    Mean_Cd_standing_crop_mg_m2 = mean(Cd_standing_crop_mg_m2, na.rm = TRUE),
    SD_Cd_standing_crop_mg_m2 = sd(Cd_standing_crop_mg_m2, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Site, Compartment, Date)

# Calculate delta standing crop between consecutive sampling dates
delta_standing_crop <- standing_crop_means %>%
  group_by(Site, Compartment) %>%
  arrange(Date) %>%
  mutate(
    # Previous date and standing crop
    Prev_date = lag(Date),
    Prev_standing_crop = lag(Mean_Cd_standing_crop_mg_m2),
    
    # Time difference (days)
    Delta_days = as.numeric(Date - Prev_date),
    
    # Change in standing crop (mg/m²)
    Delta_standing_crop_mg_m2 = Mean_Cd_standing_crop_mg_m2 - Prev_standing_crop,
    
    # Biological flux rate (mg/m²/day)
    Bio_flux_rate_mg_m2_day = ifelse(Delta_days > 0, 
                                     Delta_standing_crop_mg_m2 / Delta_days, 
                                     NA),
    
    # Scale to river reach (mg/km/day)
    Bio_flux_mg_km_day = Bio_flux_rate_mg_m2_day * Width_m * 1000,
    
    # Classify as uptake (positive) or release (negative)
    Flux_direction = case_when(
      Bio_flux_mg_km_day > 0 ~ "Net Uptake",
      Bio_flux_mg_km_day < 0 ~ "Net Release", 
      Bio_flux_mg_km_day == 0 ~ "No Change",
      TRUE ~ "Unknown"
    )
  ) %>%
  filter(!is.na(Prev_date)) %>%  # Remove first date (no previous to compare)
  ungroup() %>%
  arrange(Date, Distance_km, Compartment)

# Summary table by time period and location
flux_summary_by_period <- delta_standing_crop %>%
  group_by(Date, Site, Distance_km) %>%
  summarise(
    Date_range = paste(min(Prev_date), "to", min(Date)),
    Delta_days = first(Delta_days),
    n_compartments = n(),
    Compartments = paste(sort(Compartment), collapse = ", "),
    
    # Total biological flux across compartments
    Total_bio_flux_mg_km_day = sum(Bio_flux_mg_km_day, na.rm = TRUE),
    Mean_bio_flux_mg_km_day = mean(Bio_flux_mg_km_day, na.rm = TRUE),
    
    # Count direction of fluxes
    n_uptake = sum(Flux_direction == "Net Uptake", na.rm = TRUE),
    n_release = sum(Flux_direction == "Net Release", na.rm = TRUE),
    
    # Range of fluxes
    Min_flux_mg_km_day = min(Bio_flux_mg_km_day, na.rm = TRUE),
    Max_flux_mg_km_day = max(Bio_flux_mg_km_day, na.rm = TRUE),
    
    .groups = 'drop'
  ) %>%
  arrange(Date, Distance_km)

# Summary by site and compartment (across all time periods)
flux_summary_by_compartment <- delta_standing_crop %>%
  group_by(Site, Distance_km, Compartment) %>%
  summarise(
    n_time_periods = n(),
    Date_range = paste(min(Date), "to", max(Date)),
    
    # Mean biological flux over time
    Mean_bio_flux_mg_km_day = mean(Bio_flux_mg_km_day, na.rm = TRUE),
    SD_bio_flux_mg_km_day = sd(Bio_flux_mg_km_day, na.rm = TRUE),
    Min_bio_flux_mg_km_day = min(Bio_flux_mg_km_day, na.rm = TRUE),
    Max_bio_flux_mg_km_day = max(Bio_flux_mg_km_day, na.rm = TRUE),
    
    # Directional summary
    pct_uptake = sum(Flux_direction == "Net Uptake", na.rm = TRUE) / n() * 100,
    pct_release = sum(Flux_direction == "Net Release", na.rm = TRUE) / n() * 100,
    
    # Magnitude of changes
    Mean_abs_flux_mg_km_day = mean(abs(Bio_flux_mg_km_day), na.rm = TRUE),
    
    .groups = 'drop'
  ) %>%
  arrange(Distance_km, Compartment)

# Wide format table for temporal patterns
flux_wide_temporal <- delta_standing_crop %>%
  select(Date, Site, Compartment, Bio_flux_mg_km_day) %>%
  pivot_wider(
    names_from = Compartment,
    values_from = Bio_flux_mg_km_day,
    names_prefix = "Flux_"
  ) %>%
  arrange(Date, Site)

# Reach-scale biological balance (uptake vs release)
reach_balance <- delta_standing_crop %>%
  group_by(Site, Distance_km) %>%
  summarise(
    n_observations = n(),
    
    # Net biological flux (sum across compartments and time)
    Total_net_flux_mg_km_day = sum(Bio_flux_mg_km_day, na.rm = TRUE),
    Mean_net_flux_mg_km_day = mean(Bio_flux_mg_km_day, na.rm = TRUE),
    
    # Uptake vs release balance
    Total_uptake = sum(ifelse(Bio_flux_mg_km_day > 0, Bio_flux_mg_km_day, 0), na.rm = TRUE),
    Total_release = sum(ifelse(Bio_flux_mg_km_day < 0, abs(Bio_flux_mg_km_day), 0), na.rm = TRUE),
    Net_balance = Total_uptake - Total_release,
    
    # Classification
    Site_classification = case_when(
      Net_balance > 10 ~ "Net Biological Sink",
      Net_balance < -10 ~ "Net Biological Source", 
      TRUE ~ "Biological Equilibrium"
    ),
    
    .groups = 'drop'
  ) %>%
  arrange(Distance_km)

# Print results
cat("=== DELTA STANDING CROP BIOLOGICAL FLUX ANALYSIS ===\n")
cat("Method: Flux = Δ(Cd standing crop) / Δt\n")
cat("Positive flux = Net uptake, Negative flux = Net release\n\n")

cat("=== BIOLOGICAL FLUX BY TIME PERIOD AND SITE ===\n")
print(flux_summary_by_period)

cat("\n=== BIOLOGICAL FLUX BY COMPARTMENT (TEMPORAL AVERAGE) ===\n")
print(flux_summary_by_compartment)

cat("\n=== REACH-SCALE BIOLOGICAL BALANCE ===\n")
print(reach_balance)

cat("\n=== TEMPORAL PATTERNS (WIDE FORMAT) ===\n")
print(flux_wide_temporal)

cat("\n=== DETAILED DELTA CALCULATIONS ===\n")
detailed_output <- delta_standing_crop %>%
  select(Date, Site, Compartment, Prev_date, Delta_days, 
         Mean_Cd_standing_crop_mg_m2, Prev_standing_crop,
         Delta_standing_crop_mg_m2, Bio_flux_mg_km_day, Flux_direction)
print(detailed_output)

# Export results
write.csv(flux_summary_by_period, "biological_flux_by_period.csv", row.names = FALSE)
write.csv(flux_summary_by_compartment, "biological_flux_by_compartment_temporal.csv", row.names = FALSE)
write.csv(reach_balance, "reach_scale_biological_balance.csv", row.names = FALSE)
write.csv(flux_wide_temporal, "biological_flux_temporal_wide.csv", row.names = FALSE)
write.csv(detailed_output, "delta_standing_crop_detailed.csv", row.names = FALSE)

cat("\n=== FILES EXPORTED ===\n")
cat("- biological_flux_by_period.csv\n")
cat("- biological_flux_by_compartment_temporal.csv\n")
cat("- reach_scale_biological_balance.csv\n")
cat("- biological_flux_temporal_wide.csv\n")
cat("- delta_standing_crop_detailed.csv\n")

cat("\n=== INTERPRETATION GUIDE ===\n")
cat("Positive flux = Algae are net SINKS (removing Cd from water)\n")
cat("Negative flux = Algae are net SOURCES (releasing Cd to water)\n")
cat("Magnitude indicates strength of biological Cd demand\n")