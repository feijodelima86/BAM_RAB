# Biological Cadmium Flux Analysis - Upper Clark Fork River
# Parse biological contributions to uptake by date and site

library(dplyr)
library(tidyr)
library(lubridate)

# Read the datasets
# Assuming water_data.csv and algal_data.csv are your input files
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
  # Convert date formats if needed
  mutate(
    Date = as.Date(SAMPLING_DATE, format = "%m/%d/%Y"),
    Site = SITE,
    Compartment = SAMPLE_DESCRIPTOR,
    Biomass_g_m2 = OM.AREA.g.m2,
    Cd_mg_g = Cd,
    SEP_days = SEP
  ) %>%
  # Remove rows with missing critical data
  filter(!is.na(Biomass_g_m2), !is.na(Cd_mg_g), !is.na(SEP_days)) %>%
  # Add site width information
  left_join(site_info, by = "Site") %>%
  select(Date, Site, Compartment, Biomass_g_m2, Cd_mg_g, SEP_days, Width_m, Distance_km)

# Calculate biological flux for each sample
biological_flux <- algal_clean %>%
  mutate(
    # Calculate Cd content per unit area (mg/m²)
    Cd_content_mg_m2 = Biomass_g_m2 * Cd_mg_g,
    
    # Calculate daily turnover rate (1/days)
    Turnover_rate_per_day = 1 / SEP_days,
    
    # Calculate daily Cd turnover (mg/m²/day)
    Daily_Cd_turnover_mg_m2_day = Cd_content_mg_m2 * Turnover_rate_per_day,
    
    # Calculate biological flux (mg/km/day)
    # Convert from m² to km: multiply by width (m) × 1000 (m/km)
    Bio_flux_mg_km_day = Daily_Cd_turnover_mg_m2_day * Width_m * 1000
  )

# Summarize by Date, Site, and Compartment
bio_flux_by_compartment <- biological_flux %>%
  group_by(Date, Site, Compartment, Distance_km) %>%
  summarise(
    n_samples = n(),
    Mean_biomass_g_m2 = mean(Biomass_g_m2, na.rm = TRUE),
    Mean_Cd_mg_g = mean(Cd_mg_g, na.rm = TRUE),
    Mean_SEP_days = mean(SEP_days, na.rm = TRUE),
    Mean_bio_flux_mg_km_day = mean(Bio_flux_mg_km_day, na.rm = TRUE),
    SD_bio_flux_mg_km_day = sd(Bio_flux_mg_km_day, na.rm = TRUE),
    Min_bio_flux_mg_km_day = min(Bio_flux_mg_km_day, na.rm = TRUE),
    Max_bio_flux_mg_km_day = max(Bio_flux_mg_km_day, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Date, Distance_km, Compartment)

# Summarize total biological flux by Date and Site (across all compartments)
bio_flux_by_site_date <- biological_flux %>%
  group_by(Date, Site, Distance_km) %>%
  summarise(
    n_samples = n(),
    n_compartments = n_distinct(Compartment),
    Compartments = paste(sort(unique(Compartment)), collapse = ", "),
    Total_bio_flux_mg_km_day = sum(Bio_flux_mg_km_day, na.rm = TRUE),
    Mean_bio_flux_mg_km_day = mean(Bio_flux_mg_km_day, na.rm = TRUE),
    SD_bio_flux_mg_km_day = sd(Bio_flux_mg_km_day, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Date, Distance_km)

# Summary statistics across the entire dataset
overall_summary <- biological_flux %>%
  group_by(Site, Compartment) %>%
  summarise(
    n_dates = n_distinct(Date),
    n_samples = n(),
    Mean_bio_flux_mg_km_day = mean(Bio_flux_mg_km_day, na.rm = TRUE),
    SD_bio_flux_mg_km_day = sd(Bio_flux_mg_km_day, na.rm = TRUE),
    Min_bio_flux_mg_km_day = min(Bio_flux_mg_km_day, na.rm = TRUE),
    Max_bio_flux_mg_km_day = max(Bio_flux_mg_km_day, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(Site, Compartment)

# Create a wide format table for easier comparison
bio_flux_wide <- bio_flux_by_compartment %>%
  select(Date, Site, Distance_km, Compartment, Mean_bio_flux_mg_km_day) %>%
  pivot_wider(
    names_from = Compartment,
    values_from = Mean_bio_flux_mg_km_day,
    names_prefix = "BioFlux_"
  ) %>%
  mutate(
    Total_BioFlux_mg_km_day = rowSums(select(., starts_with("BioFlux_")), na.rm = TRUE)
  ) %>%
  arrange(Date, Distance_km)

# Add water chemistry data for comparison (if needed)
if("UF_Cd" %in% colnames(water_data)) {
  water_clean <- water_data %>%
    mutate(
      Date = as.Date(Date),
      Site = Site,
      UF_Cd_mg_L = UF_Cd,
      Discharge_m3_s = discharge_m3s
    ) %>%
    select(Date, Site, UF_Cd_mg_L, Discharge_m3_s) %>%
    filter(!is.na(Date), !is.na(Site))
  
  # Merge with biological flux data
  bio_vs_water <- bio_flux_by_site_date %>%
    left_join(water_clean, by = c("Date", "Site")) %>%
    mutate(
      # Calculate water column Cd flux (very rough estimate)
      Water_Cd_flux_mg_km_day = ifelse(!is.na(UF_Cd_mg_L) & !is.na(Discharge_m3_s),
                                       UF_Cd_mg_L * Discharge_m3_s * 86400 * 1000 / 1000, # Convert to mg/km/day
                                       NA),
      # Calculate biological flux as percentage of water flux
      Bio_pct_of_water = ifelse(!is.na(Water_Cd_flux_mg_km_day) & Water_Cd_flux_mg_km_day > 0,
                                (Total_bio_flux_mg_km_day / Water_Cd_flux_mg_km_day) * 100,
                                NA)
    )
}

# Print results
cat("=== BIOLOGICAL FLUX BY COMPARTMENT ===\n")
print(bio_flux_by_compartment)

cat("\n=== BIOLOGICAL FLUX BY SITE AND DATE ===\n")
print(bio_flux_by_site_date)

cat("\n=== BIOLOGICAL FLUX - WIDE FORMAT ===\n")
print(bio_flux_wide)

cat("\n=== OVERALL SUMMARY BY SITE AND COMPARTMENT ===\n")
print(overall_summary)

# Export to CSV files
write.csv(bio_flux_by_compartment, "biological_flux_by_compartment.csv", row.names = FALSE)
write.csv(bio_flux_by_site_date, "biological_flux_by_site_date.csv", row.names = FALSE)
write.csv(bio_flux_wide, "biological_flux_wide_format.csv", row.names = FALSE)
write.csv(overall_summary, "biological_flux_summary.csv", row.names = FALSE)

# Optional: biological vs water comparison
if(exists("bio_vs_water")) {
  cat("\n=== BIOLOGICAL vs WATER COLUMN FLUX ===\n")
  print(bio_vs_water)
  write.csv(bio_vs_water, "biological_vs_water_flux.csv", row.names = FALSE)
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files exported:\n")
cat("- biological_flux_by_compartment.csv\n")
cat("- biological_flux_by_site_date.csv\n") 
cat("- biological_flux_wide_format.csv\n")
cat("- biological_flux_summary.csv\n")
if(exists("bio_vs_water")) {
  cat("- biological_vs_water_flux.csv\n")
}