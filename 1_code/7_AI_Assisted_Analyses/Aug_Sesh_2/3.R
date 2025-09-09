# Biological Cadmium Flux vs Water Column Analysis
# Test relationships between algae growth/decline and Cd fractions

library(dplyr)
library(tidyr)
library(lubridate)
library(broom)

# Read datasets
algal_data <- read.csv("2_incremental/AI_Subsets/Cd_redux.csv")
water_data <- read.csv("2_incremental/AI_Subsets/integrated_parsed.csv")

# Site characteristics  
site_info <- data.frame(
  Site = c("WS", "DL", "GR", "GC", "BG", "BN"),
  Width_m = c(11, 24, 19, 35, 27, 31),
  Distance_km = c(0, 44, 64.8, 89.2, 144.3, 167.8)
)

# First, let's check what columns exist in water_data
cat("Checking water data columns:\n")
print(colnames(water_data))

# Process water chemistry data
water_clean <- water_data %>%
  mutate(
    Date = as.Date(Date),
    Site = Site,
    
    # Cd concentrations (mg/L)
    UF_Cd_mg_L = UF_Cd,           # Dissolved (ultrafiltered)
    F_Cd_mg_L = F_Cd,             # Dissolved + colloidal (filtered)
    W_Cd_mg_L = W_Cd,             # Total (whole water)
    
    # Calculate fractions
    Dissolved_Cd_mg_L = UF_Cd_mg_L,
    Colloidal_Cd_mg_L = F_Cd_mg_L - UF_Cd_mg_L,
    Suspended_Cd_mg_L = W_Cd_mg_L - F_Cd_mg_L,
    
    # Water column Cd flux (mg/km/day) - rough estimate
    Discharge_m3_s = discharge_m3s,
    Total_Cd_flux_mg_km_day = ifelse(!is.na(W_Cd_mg_L) & !is.na(Discharge_m3_s),
                                     W_Cd_mg_L * Discharge_m3_s * 86400 * 1000 / 1000,
                                     NA),
    Dissolved_Cd_flux_mg_km_day = ifelse(!is.na(Dissolved_Cd_mg_L) & !is.na(Discharge_m3_s),
                                         Dissolved_Cd_mg_L * Discharge_m3_s * 86400 * 1000 / 1000,
                                         NA)
  ) %>%
  # Keep original Distance_km and add Width_m from site_info
  left_join(site_info %>% select(Site, Width_m), by = "Site") %>%
  select(Date, Site, Distance_km, Width_m, 
         Dissolved_Cd_mg_L, Colloidal_Cd_mg_L, Suspended_Cd_mg_L,
         UF_Cd_mg_L, F_Cd_mg_L, W_Cd_mg_L,
         Dissolved_pct, Colloidal_pct, Suspended_pct,
         Discharge_m3_s, Total_Cd_flux_mg_km_day, Dissolved_Cd_flux_mg_km_day) %>%
  arrange(Site, Date)

# Process algal data and calculate delta standing crop
cat("Checking algal data columns after join:\n")
algal_clean <- algal_data %>%
  mutate(
    Date = as.Date(SAMPLING_DATE, format = "%m/%d/%Y"),
    Site = SITE,
    Compartment = SAMPLE_DESCRIPTOR,
    Biomass_g_m2 = OM.AREA.g.m2,
    Cd_mg_g = Cd,
    Cd_standing_crop_mg_m2 = Biomass_g_m2 * Cd_mg_g
  ) %>%
  filter(!is.na(Date), !is.na(Site), !is.na(Compartment), 
         !is.na(Cd_standing_crop_mg_m2)) %>%
  left_join(site_info, by = "Site")

print(colnames(algal_clean))

# Calculate mean standing crop per site-compartment-date
standing_crop_means <- algal_clean %>%
  group_by(Date, Site, Distance_km, Compartment, Width_m) %>%
  summarise(
    Mean_biomass_g_m2 = mean(Biomass_g_m2, na.rm = TRUE),
    Mean_Cd_standing_crop_mg_m2 = mean(Cd_standing_crop_mg_m2, na.rm = TRUE),
    .groups = 'drop'
  )

# Calculate delta standing crop (biological flux)
delta_standing_crop <- standing_crop_means %>%
  group_by(Site, Compartment) %>%
  arrange(Date) %>%
  mutate(
    Prev_date = lag(Date),
    Prev_standing_crop = lag(Mean_Cd_standing_crop_mg_m2),
    Prev_biomass = lag(Mean_biomass_g_m2),
    
    Delta_days = as.numeric(Date - Prev_date),
    Delta_standing_crop_mg_m2 = Mean_Cd_standing_crop_mg_m2 - Prev_standing_crop,
    Delta_biomass_g_m2 = Mean_biomass_g_m2 - Prev_biomass,
    
    # Biological flux (mg/km/day)
    Bio_flux_mg_km_day = ifelse(Delta_days > 0, 
                                Delta_standing_crop_mg_m2 / Delta_days * Width_m * 1000, 
                                NA),
    
    # Biomass change rate (g/m²/day)
    Biomass_change_rate = ifelse(Delta_days > 0,
                                 Delta_biomass_g_m2 / Delta_days,
                                 NA),
    
    # Classify biological activity
    Bio_activity = case_when(
      Bio_flux_mg_km_day > 5 ~ "Strong Uptake",
      Bio_flux_mg_km_day > 0 ~ "Weak Uptake", 
      Bio_flux_mg_km_day < -5 ~ "Strong Release",
      Bio_flux_mg_km_day < 0 ~ "Weak Release",
      TRUE ~ "No Change"
    ),
    
    # Classify biomass change
    Biomass_trend = case_when(
      Biomass_change_rate > 1 ~ "Growing",
      Biomass_change_rate > 0 ~ "Slow Growth",
      Biomass_change_rate < -1 ~ "Declining", 
      Biomass_change_rate < 0 ~ "Slow Decline",
      TRUE ~ "Stable"
    )
  ) %>%
  filter(!is.na(Prev_date)) %>%
  ungroup()

# Aggregate biological fluxes by site and date (across all compartments)
bio_flux_total <- delta_standing_crop %>%
  group_by(Date, Site, Distance_km) %>%
  summarise(
    n_compartments = n(),
    Total_bio_flux_mg_km_day = sum(Bio_flux_mg_km_day, na.rm = TRUE),
    Mean_bio_flux_mg_km_day = mean(Bio_flux_mg_km_day, na.rm = TRUE),
    Total_biomass_change_rate = sum(Biomass_change_rate, na.rm = TRUE),
    
    # Dominant biological activity
    Dominant_activity = names(sort(table(Bio_activity), decreasing = TRUE))[1],
    Dominant_biomass_trend = names(sort(table(Biomass_trend), decreasing = TRUE))[1],
    
    .groups = 'drop'
  )

# Merge biological and water chemistry data
bio_water_merged <- bio_flux_total %>%
  left_join(water_clean, by = c("Date", "Site", "Distance_km")) %>%
  filter(!is.na(Total_bio_flux_mg_km_day), 
         !is.na(Dissolved_Cd_mg_L) | !is.na(Colloidal_Cd_mg_L) | !is.na(Suspended_Cd_mg_L))

# Calculate flux magnitude comparisons
flux_comparison <- bio_water_merged %>%
  mutate(
    # Biological flux as % of water column flux
    Bio_pct_of_total_flux = ifelse(!is.na(Total_Cd_flux_mg_km_day) & Total_Cd_flux_mg_km_day > 0,
                                   abs(Total_bio_flux_mg_km_day) / Total_Cd_flux_mg_km_day * 100,
                                   NA),
    
    # Flux magnitude categories
    Water_flux_category = case_when(
      Total_Cd_flux_mg_km_day > 500000 ~ "Very High",
      Total_Cd_flux_mg_km_day > 100000 ~ "High", 
      Total_Cd_flux_mg_km_day > 10000 ~ "Moderate",
      Total_Cd_flux_mg_km_day > 1000 ~ "Low",
      TRUE ~ "Very Low"
    ),
    
    Bio_flux_category = case_when(
      abs(Total_bio_flux_mg_km_day) > 100 ~ "High Bio Activity",
      abs(Total_bio_flux_mg_km_day) > 10 ~ "Moderate Bio Activity",
      abs(Total_bio_flux_mg_km_day) > 1 ~ "Low Bio Activity", 
      TRUE ~ "Negligible Bio Activity"
    )
  )

# Test relationships between biological activity and Cd fractions
correlation_analysis <- bio_water_merged %>%
  filter(!is.na(Total_bio_flux_mg_km_day)) %>%
  group_by(Site, Distance_km) %>%
  summarise(
    n_observations = n(),
    
    # Correlations with Cd concentrations
    cor_bio_dissolved = ifelse(n() > 3 & var(Total_bio_flux_mg_km_day, na.rm=T) > 0,
                               cor(Total_bio_flux_mg_km_day, Dissolved_Cd_mg_L, use="complete.obs"),
                               NA),
    cor_bio_colloidal = ifelse(n() > 3 & var(Total_bio_flux_mg_km_day, na.rm=T) > 0,
                               cor(Total_bio_flux_mg_km_day, Colloidal_Cd_mg_L, use="complete.obs"),
                               NA),
    cor_bio_suspended = ifelse(n() > 3 & var(Total_bio_flux_mg_km_day, na.rm=T) > 0,
                               cor(Total_bio_flux_mg_km_day, Suspended_Cd_mg_L, use="complete.obs"),
                               NA),
    
    # Correlations with biomass change
    cor_biomass_dissolved = ifelse(n() > 3 & var(Total_biomass_change_rate, na.rm=T) > 0,
                                   cor(Total_biomass_change_rate, Dissolved_Cd_mg_L, use="complete.obs"),
                                   NA),
    cor_biomass_colloidal = ifelse(n() > 3 & var(Total_biomass_change_rate, na.rm=T) > 0,
                                   cor(Total_biomass_change_rate, Colloidal_Cd_mg_L, use="complete.obs"),
                                   NA),
    cor_biomass_suspended = ifelse(n() > 3 & var(Total_biomass_change_rate, na.rm=T) > 0,
                                   cor(Total_biomass_change_rate, Suspended_Cd_mg_L, use="complete.obs"),
                                   NA),
    .groups = 'drop'
  ) %>%
  arrange(Distance_km)

# Detailed compartment-level analysis
compartment_analysis <- delta_standing_crop %>%
  left_join(water_clean, by = c("Date", "Site", "Distance_km")) %>%
  filter(!is.na(Bio_flux_mg_km_day), !is.na(Dissolved_Cd_mg_L)) %>%
  group_by(Site, Distance_km, Compartment) %>%
  summarise(
    n_obs = n(),
    Mean_bio_flux = mean(Bio_flux_mg_km_day, na.rm = TRUE),
    
    # Test your hypothesis: growing algae vs dissolved Cd
    cor_growth_dissolved = ifelse(n() > 3,
                                  cor(Biomass_change_rate, Dissolved_Cd_mg_L, use="complete.obs"),
                                  NA),
    
    # Dying algae vs suspended Cd
    cor_decline_suspended = ifelse(n() > 3,
                                   cor(-Biomass_change_rate, Suspended_Cd_mg_L, use="complete.obs"),
                                   NA),
    
    .groups = 'drop'
  ) %>%
  arrange(Distance_km, Compartment)

# Summary by contamination level (based on site)
contamination_summary <- flux_comparison %>%
  mutate(
    Contamination_level = case_when(
      Site %in% c("DL", "GR") ~ "High (Contaminated Reach)",
      Site %in% c("WS") ~ "Low (Upstream)",
      TRUE ~ "Moderate (Downstream)"
    )
  ) %>%
  group_by(Contamination_level) %>%
  summarise(
    n_observations = n(),
    Mean_dissolved_Cd = mean(Dissolved_Cd_mg_L, na.rm = TRUE),
    Mean_bio_flux = mean(Total_bio_flux_mg_km_day, na.rm = TRUE),
    Mean_bio_pct_water = mean(Bio_pct_of_total_flux, na.rm = TRUE),
    Typical_bio_activity = names(sort(table(Dominant_activity), decreasing = TRUE))[1],
    .groups = 'drop'
  )

# Print comprehensive results
cat("=== BIOLOGICAL vs WATER COLUMN Cd ANALYSIS ===\n")
cat("Testing hypothesis: Growing algae ↔ Dissolved Cd, Dying algae ↔ Suspended Cd\n\n")

cat("=== FLUX MAGNITUDE COMPARISON ===\n")
print(flux_comparison %>% 
        select(Date, Site, Total_bio_flux_mg_km_day, Total_Cd_flux_mg_km_day, 
               Bio_pct_of_total_flux, Dominant_activity) %>%
        arrange(Date, Site))

cat("\n=== CORRELATION ANALYSIS BY SITE ===\n")
print(correlation_analysis)

cat("\n=== COMPARTMENT-LEVEL ANALYSIS ===\n")
print(compartment_analysis)

cat("\n=== CONTAMINATION LEVEL SUMMARY ===\n")
print(contamination_summary)

cat("\n=== HYPOTHESIS TESTING RESULTS ===\n")
hypothesis_test <- bio_water_merged %>%
  filter(!is.na(Total_bio_flux_mg_km_day), !is.na(Dissolved_Cd_mg_L)) %>%
  summarise(
    # Overall correlations
    Overall_cor_bio_dissolved = cor(Total_bio_flux_mg_km_day, Dissolved_Cd_mg_L, use="complete.obs"),
    Overall_cor_bio_suspended = cor(Total_bio_flux_mg_km_day, Suspended_Cd_mg_L, use="complete.obs"),
    Overall_cor_biomass_dissolved = cor(Total_biomass_change_rate, Dissolved_Cd_mg_L, use="complete.obs"),
    
    # Your hypothesis predictions
    Hypothesis_1 = "Growing algae should correlate negatively with dissolved Cd",
    Hypothesis_2 = "Dying algae should correlate positively with suspended Cd",
    
    # Test results
    H1_result = ifelse(Overall_cor_biomass_dissolved < -0.2, "SUPPORTED", 
                       ifelse(Overall_cor_biomass_dissolved > 0.2, "OPPOSITE", "NO SIGNAL")),
    H2_result = ifelse(-Overall_cor_biomass_dissolved > 0.2 & 
                         cor(Total_biomass_change_rate, Suspended_Cd_mg_L, use="complete.obs") < -0.2, 
                       "SUPPORTED", "UNCLEAR")
  )
print(hypothesis_test)

# Export comprehensive results
write.csv(flux_comparison, "bio_vs_water_flux_comparison.csv", row.names = FALSE)
write.csv(correlation_analysis, "correlation_analysis_by_site.csv", row.names = FALSE)
write.csv(compartment_analysis, "compartment_level_analysis.csv", row.names = FALSE)
write.csv(bio_water_merged, "merged_bio_water_data.csv", row.names = FALSE)

cat("\n=== INTERPRETATION ===\n")
cat("Look for:\n")
cat("- Negative correlation: biomass growth vs dissolved Cd (uptake)\n") 
cat("- Positive correlation: biomass decline vs suspended Cd (release)\n")
cat("- Compartment differences: EPIL vs EPIP vs FILA\n")
cat("- Site differences: contaminated vs clean reaches\n")
cat("- Temporal patterns: seasonal algae-Cd relationships\n")