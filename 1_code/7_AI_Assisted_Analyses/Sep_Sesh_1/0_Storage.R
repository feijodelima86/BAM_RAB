# Load required libraries
library(dplyr)
library(tidyr)

# Setup - assuming your data is in 'df2'
raw_data <- df1

# Remove rows with all NA values
raw_data <- raw_data[!is.na(raw_data$SITE), ]

# Step 1: Aggregate algae data by site, date, and algae type
algae_summary <- raw_data %>%
  filter(!is.na(SAMPLE_DESCRIPTOR), !is.na(Cd), !is.na(OM.AREA.g.m2)) %>%
  group_by(SITE, SAMPLING_DATE, SAMPLE_DESCRIPTOR) %>%
  summarise(
    Cd_mean = mean(Cd, na.rm = TRUE),
    OM_AREA_mean = mean(OM.AREA.g.m2, na.rm = TRUE),
    .groups = 'drop'
  )

# Step 2: Reshape to wide format
algae_wide <- algae_summary %>%
  pivot_wider(
    names_from = SAMPLE_DESCRIPTOR,
    values_from = c(Cd_mean, OM_AREA_mean),
    names_sep = "_"
  )

# Step 3: Get water chemistry data (one row per site/date)
water_data <- raw_data %>%
  filter(!is.na(SITE), !is.na(SAMPLING_DATE)) %>%
  group_by(SITE, SAMPLING_DATE) %>%
  summarise(
    UF_Cd = first(UF_Cd114),
    F_Cd = first(F_Cd114), 
    W_Cd = first(W_Cd114),
    discharge_m3s = first(discharge_m3s),
    Distance_km = first(Distance_km),
    Position = first(Position),
    .groups = 'drop'
  ) %>%
  filter(!is.na(UF_Cd), !is.na(discharge_m3s))

# Step 4: Combine algae and water data
data <- water_data %>%
  left_join(algae_wide, by = c("SITE", "SAMPLING_DATE")) %>%
  rename(
    Site = SITE,
    Date = SAMPLING_DATE,
    Cd_EPIL = Cd_mean_EPIL,
    OM_AREA_g_m2_EPIL = OM_AREA_mean_EPIL,
    Cd_EPIP = Cd_mean_EPIP, 
    OM_AREA_g_m2_EPIP = OM_AREA_mean_EPIP,
    Cd_FILA = Cd_mean_FILA,
    OM_AREA_g_m2_FILA = OM_AREA_mean_FILA
  )

# Setup parameters
site_order <- c("WS", "DL", "GR", "GC", "BG", "BN")
data$Site <- factor(data$Site, levels = site_order)
dates <- unique(data$Date)
n_dates <- length(dates)
reach_names <- c("WS-DL", "DL-GR", "GR-GC", "GC-BG", "BG-BN")
river_width <- c(11, 24, 19, 35, 27, 31)
reach_lengths <- c(44.0, 20.8, 24.4, 55.1, 23.5)

# 1. UF Cd flux per reach matrix
daily_uf_cd_flux <- matrix(NA, n_dates, length(site_order))
colnames(daily_uf_cd_flux) <- site_order
rownames(daily_uf_cd_flux) <- dates

for(d in 1:n_dates) {
  date_data <- data[data$Date == dates[d], ]
  for(s in 1:length(site_order)) {
    site_data <- date_data[date_data$Site == site_order[s], ]
    if(nrow(site_data) > 0) {
      uf_cd <- site_data$UF_Cd[1]
      discharge <- site_data$discharge_m3s[1]
      if(!is.na(uf_cd) & !is.na(discharge)) {
        daily_uf_cd_flux[d, s] <- uf_cd * discharge * 1000 * 86400
      }
    }
  }
}

# Calculate reach flux (average of upstream and downstream sites)
reach_uf_flux <- matrix(NA, n_dates, 5)
colnames(reach_uf_flux) <- reach_names
rownames(reach_uf_flux) <- dates

for(d in 1:n_dates) {
  reach_uf_flux[d, 1] <- mean(c(daily_uf_cd_flux[d, 1], daily_uf_cd_flux[d, 2]), na.rm = TRUE)
  reach_uf_flux[d, 2] <- mean(c(daily_uf_cd_flux[d, 2], daily_uf_cd_flux[d, 3]), na.rm = TRUE)
  reach_uf_flux[d, 3] <- mean(c(daily_uf_cd_flux[d, 3], daily_uf_cd_flux[d, 4]), na.rm = TRUE)
  reach_uf_flux[d, 4] <- mean(c(daily_uf_cd_flux[d, 4], daily_uf_cd_flux[d, 5]), na.rm = TRUE)
  reach_uf_flux[d, 5] <- mean(c(daily_uf_cd_flux[d, 5], daily_uf_cd_flux[d, 6]), na.rm = TRUE)
}

# 2. Stock per reach matrix
total_cd_stock <- matrix(NA, n_dates, length(site_order))
colnames(total_cd_stock) <- site_order
rownames(total_cd_stock) <- dates

for(d in 1:n_dates) {
  date_data <- data[data$Date == dates[d], ]
  for(s in 1:length(site_order)) {
    site_data <- date_data[date_data$Site == site_order[s], ]
    if(nrow(site_data) > 0) {
      # Calculate stock for each algae type (concentration * biomass)
      epil_stock <- ifelse(!is.na(site_data$Cd_EPIL[1]) & !is.na(site_data$OM_AREA_g_m2_EPIL[1]), 
                           site_data$Cd_EPIL[1] * site_data$OM_AREA_g_m2_EPIL[1], 0)
      epip_stock <- ifelse(!is.na(site_data$Cd_EPIP[1]) & !is.na(site_data$OM_AREA_g_m2_EPIP[1]), 
                           site_data$Cd_EPIP[1] * site_data$OM_AREA_g_m2_EPIP[1], 0)
      fil_stock <- ifelse(!is.na(site_data$Cd_FILA[1]) & !is.na(site_data$OM_AREA_g_m2_FILA[1]), 
                          site_data$Cd_FILA[1] * site_data$OM_AREA_g_m2_FILA[1], 0)
      
      total_cd_stock[d, s] <- epil_stock + epip_stock + fil_stock
    }
  }
}

# Calculate reach areas
reach_areas <- numeric(5)
for(r in 1:5) {
  avg_width <- (river_width[r] + river_width[r+1]) / 2
  reach_areas[r] <- reach_lengths[r] * 1000 * avg_width  # km * 1000 * m = m²
}

# Calculate reach stock (average stock per m² * reach area)
reach_stock <- matrix(NA, n_dates, 5)
colnames(reach_stock) <- reach_names
rownames(reach_stock) <- dates

for(d in 1:n_dates) {
  for(r in 1:5) {
    avg_stock <- mean(c(total_cd_stock[d, r], total_cd_stock[d, r+1]), na.rm = TRUE)
    reach_stock[d, r] <- avg_stock * reach_areas[r]
  }
}

# 3. % Cd storage relative to flux
percent_storage <- (reach_stock / reach_uf_flux) * 100

# Results
cat("Data Summary:\n")
cat("Sites processed:", paste(unique(data$Site), collapse = ", "), "\n")
cat("Dates processed:", paste(dates, collapse = ", "), "\n")
cat("Number of algae samples aggregated:", nrow(algae_summary), "\n\n")

cat("UF Cd Flux Matrix (g/day):\n")
print(round(reach_uf_flux, 2))

cat("\nStock Matrix (g):\n") 
print(round(reach_stock, 2))

cat("\n% Storage Matrix:\n")
print(round(percent_storage, 4))

# Optional: Create summary data frame for easier analysis
results_summary <- data.frame(
  Date = rep(dates, 5),
  Reach = rep(reach_names, each = n_dates),
  UF_Flux_g_day = as.vector(reach_uf_flux),
  Stock_g = as.vector(reach_stock), 
  Percent_Storage = as.vector(percent_storage)
)

cat("\nSummary statistics by reach:\n")
print(
  results_summary %>%
    group_by(Reach) %>%
    summarise(
      Mean_Flux = round(mean(UF_Flux_g_day, na.rm = TRUE), 2),
      Mean_Stock = round(mean(Stock_g, na.rm = TRUE), 2),
      Mean_Percent_Storage = round(mean(Percent_Storage, na.rm = TRUE), 4),
      .groups = 'drop'
    )
)