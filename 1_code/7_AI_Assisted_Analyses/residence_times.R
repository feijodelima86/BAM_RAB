# Simple R code for algal Cd flux analysis
# Works directly on objects without complex functions

# Read data (adjust path as needed)
data <- df2

# Define site order and distances
site_order <- c("WS", "DL", "GR", "GC", "BG", "BN")
site_distances <- c(0, 44, 64.8, 89.2, 144.3, 167.8)

# Create ordered factors for sites
data$Site <- factor(data$Site, levels = site_order)

# Sort data by date and site
data <- data[order(data$Date, data$Site), ]

# Get unique dates
dates <- unique(data$Date)

# Initialize storage matrices
n_dates <- length(dates)
n_reaches <- length(site_order) - 1

# Calculate daily Cd flux through each site (mg Cd/day)
# Simple: Cd concentration × daily discharge

# Initialize storage for daily fluxes
daily_cd_flux <- matrix(NA, n_dates, length(site_order))
colnames(daily_cd_flux) <- site_order
rownames(daily_cd_flux) <- dates

# Calculate for each date and site
for(d in 1:n_dates) {
  date_data <- data[data$Date == dates[d], ]
  
  for(s in 1:length(site_order)) {
    site_data <- date_data[date_data$Site == site_order[s], ]
    
    if(nrow(site_data) > 0) {
      cd_conc <- site_data$W_Cd[1]        # mg/L
      discharge <- site_data$discharge_m3s[1]  # m3/s
      
      if(!is.na(cd_conc) & !is.na(discharge)) {
        # Convert to mg Cd/day
        daily_cd_flux[d, s] <- cd_conc * discharge * 1000 * 86400
        # cd_conc (mg/L) × discharge (m3/s) × 1000 (L/m3) × 86400 (s/day)
      }
    }
  }
}

# Print results
cat("\n=== DAILY CD FLUX THROUGH SITES (mg Cd/day) ===\n")
print(round(daily_cd_flux, 2))

# Calculate averages
avg_daily_flux <- colMeans(daily_cd_flux, na.rm=TRUE)
cat("\n=== AVERAGE DAILY CD FLUX BY SITE ===\n")
print(round(avg_daily_flux, 2))

# Calculate Cd standing stock in algae (mg Cd/m2)
# Simple: Cd concentration × algal biomass for each compartment

# Initialize storage matrices
epil_cd_stock <- matrix(NA, n_dates, length(site_order))
epip_cd_stock <- matrix(NA, n_dates, length(site_order))
fil_cd_stock <- matrix(NA, n_dates, length(site_order))
total_cd_stock <- matrix(NA, n_dates, length(site_order))

colnames(epil_cd_stock) <- site_order
colnames(epip_cd_stock) <- site_order
colnames(fil_cd_stock) <- site_order
colnames(total_cd_stock) <- site_order
rownames(epil_cd_stock) <- dates
rownames(epip_cd_stock) <- dates
rownames(fil_cd_stock) <- dates
rownames(total_cd_stock) <- dates

# Calculate for each date and site
for(d in 1:n_dates) {
  date_data <- data[data$Date == dates[d], ]
  
  for(s in 1:length(site_order)) {
    site_data <- date_data[date_data$Site == site_order[s], ]
    
    if(nrow(site_data) > 0) {
      # Epilithon standing stock
      epil_cd <- site_data$Cd_EPIL[1]           # mg Cd/g biomass
      epil_biomass <- site_data$OM_AREA_g_m2_EPIL[1]  # g/m2
      
      if(!is.na(epil_cd) & !is.na(epil_biomass)) {
        epil_cd_stock[d, s] <- epil_cd * epil_biomass  # mg Cd/m2
      }
      
      # Epiphyton standing stock
      epip_cd <- site_data$Cd_EPIP[1]           # mg Cd/g biomass
      epip_biomass <- site_data$OM_AREA_g_m2_EPIP[1]  # g/m2
      
      if(!is.na(epip_cd) & !is.na(epip_biomass)) {
        epip_cd_stock[d, s] <- epip_cd * epip_biomass  # mg Cd/m2
      }
      
      # Filamentous algae standing stock
      fil_cd <- site_data$Cd_FILA[1]           # mg Cd/g biomass
      fil_biomass <- site_data$OM_AREA_g_m2_FILA[1]  # g/m2
      
      if(!is.na(fil_cd) & !is.na(fil_biomass)) {
        fil_cd_stock[d, s] <- fil_cd * fil_biomass  # mg Cd/m2
      }
      
      # Total standing stock
      total_cd_stock[d, s] <- sum(c(epil_cd_stock[d, s], epip_cd_stock[d, s], fil_cd_stock[d, s]), na.rm=TRUE)
    }
  }
}

# Print results
cat("\n=== EPILITHON CD STANDING STOCK (mg Cd/m2) ===\n")
print(round(epil_cd_stock, 4))

cat("\n=== EPIPHYTON CD STANDING STOCK (mg Cd/m2) ===\n")
print(round(epip_cd_stock, 4))

cat("\n=== FILAMENTOUS ALGAE CD STANDING STOCK (mg Cd/m2) ===\n")
print(round(fil_cd_stock, 4))

cat("\n=== TOTAL ALGAL CD STANDING STOCK (mg Cd/m2) ===\n")
print(round(total_cd_stock, 4))

# Calculate averages
avg_epil_stock <- colMeans(epil_cd_stock, na.rm=TRUE)
avg_epip_stock <- colMeans(epip_cd_stock, na.rm=TRUE)
avg_fil_stock <- colMeans(fil_cd_stock, na.rm=TRUE)
avg_total_stock <- colMeans(total_cd_stock, na.rm=TRUE)

cat("\n=== AVERAGE CD STANDING STOCK BY SITE ===\n")
stock_summary <- data.frame(
  Site = site_order,
  Epilithon = round(avg_epil_stock, 4),
  Epiphyton = round(avg_epip_stock, 4),
  Filamentous = round(avg_fil_stock, 4),
  Total = round(avg_total_stock, 4)
)
print(stock_summary)

# Fill missing values with site-specific means
total_cd_stock_filled <- total_cd_stock

for(s in 1:ncol(total_cd_stock_filled)) {
  site_mean <- mean(total_cd_stock_filled[, s], na.rm = TRUE)
  total_cd_stock_filled[is.na(total_cd_stock_filled[, s]), s] <- site_mean
}

# Define reach parameters
river_width <- c(11, 24, 19, 35, 27, 31) # WS, DL, GR, GC, BG, BN (meters)
reach_lengths <- c(44.0, 20.8, 24.4, 55.1, 23.5) # km
reach_names <- c("WS-DL", "DL-GR", "GR-GC", "GC-BG", "BG-BN")

total_cd_stock <- total_cd_stock_filled

# Calculate reach standing stock (total mg Cd per reach)
n_reaches <- length(reach_names)
reach_standing_stock <- matrix(NA, nrow = nrow(total_cd_stock), ncol = n_reaches)
colnames(reach_standing_stock) <- reach_names
rownames(reach_standing_stock) <- rownames(total_cd_stock)

# Calculate for each reach
for(r in 1:n_reaches) {
  # Get upstream and downstream site indices
  upstream_idx <- r
  downstream_idx <- r + 1
  
  # Calculate average width for this reach
  avg_width <- (river_width[upstream_idx] + river_width[downstream_idx]) / 2
  
  # Calculate reach area (km × m = km×m, convert to m2)
  reach_area <- reach_lengths[r] * 1000 * avg_width  # m2
  
  # For each date, calculate reach standing stock
  for(d in 1:nrow(total_cd_stock)) {
    # Average standing stock between upstream and downstream sites (mg Cd/m2)
    upstream_stock <- total_cd_stock[d, upstream_idx]
    downstream_stock <- total_cd_stock[d, downstream_idx]
    
    if(!is.na(upstream_stock) & !is.na(downstream_stock)) {
      avg_stock_per_m2 <- (upstream_stock + downstream_stock) / 2
      
      # Total standing stock in reach (mg Cd)
      reach_standing_stock[d, r] <- avg_stock_per_m2 * reach_area
    }
  }
}

# Print results
cat("\n=== REACH STANDING STOCK (total mg Cd per reach) ===\n")
print(round(reach_standing_stock, 2))

# Calculate reach areas for reference
reach_areas <- numeric(n_reaches)
for(r in 1:n_reaches) {
  avg_width <- (river_width[r] + river_width[r+1]) / 2
  reach_areas[r] <- reach_lengths[r] * 1000 * avg_width
}

cat("\n=== REACH AREAS (m2) ===\n")
reach_area_summary <- data.frame(
  Reach = reach_names,
  Length_km = reach_lengths,
  Avg_Width_m = round((river_width[1:5] + river_width[2:6]) / 2, 1),
  Area_m2 = round(reach_areas, 0)
)
print(reach_area_summary)

# Calculate average standing stock per reach
avg_reach_stock <- colMeans(reach_standing_stock, na.rm = TRUE)
cat("\n=== AVERAGE STANDING STOCK BY REACH ===\n")
reach_stock_summary <- data.frame(
  Reach = reach_names,
  Avg_Stock_mg_Cd = round(avg_reach_stock, 2),
  Avg_Stock_per_m2 = round(avg_reach_stock / reach_areas, 6)
)
print(reach_stock_summary)

# Create reach standing stock matrix (mg Cd/m2)
n_reaches <- 5
reach_names <- c("WS-DL", "DL-GR", "GR-GC", "GC-BG", "BG-BN")

reach_stock_matrix <- matrix(NA, nrow = nrow(total_cd_stock), ncol = n_reaches)
colnames(reach_stock_matrix) <- reach_names
rownames(reach_stock_matrix) <- rownames(total_cd_stock)

# Calculate for each reach by averaging upstream and downstream sites
for(d in 1:nrow(total_cd_stock)) {
  # WS-DL reach (average WS and DL)
  reach_stock_matrix[d, "WS-DL"] <- mean(c(total_cd_stock[d, "WS"], total_cd_stock[d, "DL"]), na.rm = TRUE)
  
  # DL-GR reach (average DL and GR)  
  reach_stock_matrix[d, "DL-GR"] <- mean(c(total_cd_stock[d, "DL"], total_cd_stock[d, "GR"]), na.rm = TRUE)
  
  # GR-GC reach (average GR and GC)
  reach_stock_matrix[d, "GR-GC"] <- mean(c(total_cd_stock[d, "GR"], total_cd_stock[d, "GC"]), na.rm = TRUE)
  
  # GC-BG reach (average GC and BG)
  reach_stock_matrix[d, "GC-BG"] <- mean(c(total_cd_stock[d, "GC"], total_cd_stock[d, "BG"]), na.rm = TRUE)
  
  # BG-BN reach (average BG and BN)
  reach_stock_matrix[d, "BG-BN"] <- mean(c(total_cd_stock[d, "BG"], total_cd_stock[d, "BN"]), na.rm = TRUE)
}

# Print the matrix
cat("\n=== REACH STANDING STOCK MATRIX (mg Cd/m2) ===\n")
print(round(reach_stock_matrix, 4))


# Define reach parameters
river_width <- c(11, 24, 19, 35, 27, 31) # WS, DL, GR, GC, BG, BN (meters)
reach_lengths <- c(44.0, 20.8, 24.4, 55.1, 23.5) # km
reach_names <- c("WS-DL", "DL-GR", "GR-GC", "GC-BG", "BG-BN")

# Calculate reach areas (m2)
reach_areas <- numeric(5)
for(r in 1:5) {
  avg_width <- (river_width[r] + river_width[r+1]) / 2
  reach_areas[r] <- reach_lengths[r] * 1000 * avg_width  # Convert km to m, then × width
}

# Create reach standing stock matrix (mg Cd/m2)
reach_stock_per_m2 <- matrix(NA, nrow = nrow(total_cd_stock), ncol = 5)
colnames(reach_stock_per_m2) <- reach_names
rownames(reach_stock_per_m2) <- rownames(total_cd_stock)

# Calculate standing stock per m2 for each reach
for(d in 1:nrow(total_cd_stock)) {
  # WS-DL reach (average WS and DL)
  reach_stock_per_m2[d, "WS-DL"] <- mean(c(total_cd_stock[d, "WS"], total_cd_stock[d, "DL"]), na.rm = TRUE)
  
  # DL-GR reach (average DL and GR)  
  reach_stock_per_m2[d, "DL-GR"] <- mean(c(total_cd_stock[d, "DL"], total_cd_stock[d, "GR"]), na.rm = TRUE)
  
  # GR-GC reach (average GR and GC)
  reach_stock_per_m2[d, "GR-GC"] <- mean(c(total_cd_stock[d, "GR"], total_cd_stock[d, "GC"]), na.rm = TRUE)
  
  # GC-BG reach (average GC and BG)
  reach_stock_per_m2[d, "GC-BG"] <- mean(c(total_cd_stock[d, "GC"], total_cd_stock[d, "BG"]), na.rm = TRUE)
  
  # BG-BN reach (average BG and BN)
  reach_stock_per_m2[d, "BG-BN"] <- mean(c(total_cd_stock[d, "BG"], total_cd_stock[d, "BN"]), na.rm = TRUE)
}

# Create total reach standing stock matrix (total mg Cd per reach)
total_reach_stock <- matrix(NA, nrow = nrow(total_cd_stock), ncol = 5)
colnames(total_reach_stock) <- reach_names
rownames(total_reach_stock) <- rownames(total_cd_stock)

# Calculate total standing stock for each reach
for(r in 1:5) {
  total_reach_stock[, r] <- reach_stock_per_m2[, r] * reach_areas[r]
}

# Print results
cat("\n=== REACH AREAS (m2) ===\n")
reach_area_df <- data.frame(
  Reach = reach_names,
  Length_km = reach_lengths,
  Avg_Width_m = round((river_width[1:5] + river_width[2:6]) / 2, 1),
  Area_m2 = round(reach_areas, 0)
)
print(reach_area_df)

cat("\n=== TOTAL REACH STANDING STOCK (mg Cd per reach) ===\n")
print(round(total_reach_stock, 2))

cat("\n=== AVERAGE TOTAL STANDING STOCK BY REACH ===\n")
avg_total_reach_stock <- colMeans(total_reach_stock, na.rm = TRUE)
print(round(avg_total_reach_stock, 2))


# Create reach daily flux matrix (mg Cd/day)
reach_daily_flux <- matrix(NA, nrow = nrow(daily_cd_flux), ncol = 5)
colnames(reach_daily_flux) <- reach_names
rownames(reach_daily_flux) <- rownames(daily_cd_flux)

# Calculate daily flux for each reach by averaging upstream and downstream sites
for(d in 1:nrow(daily_cd_flux)) {
  # WS-DL reach (average WS and DL fluxes)
  reach_daily_flux[d, "WS-DL"] <- mean(c(daily_cd_flux[d, "WS"], daily_cd_flux[d, "DL"]), na.rm = TRUE)
  
  # DL-GR reach (average DL and GR fluxes)  
  reach_daily_flux[d, "DL-GR"] <- mean(c(daily_cd_flux[d, "DL"], daily_cd_flux[d, "GR"]), na.rm = TRUE)
  
  # GR-GC reach (average GR and GC fluxes)
  reach_daily_flux[d, "GR-GC"] <- mean(c(daily_cd_flux[d, "GR"], daily_cd_flux[d, "GC"]), na.rm = TRUE)
  
  # GC-BG reach (average GC and BG fluxes)
  reach_daily_flux[d, "GC-BG"] <- mean(c(daily_cd_flux[d, "GC"], daily_cd_flux[d, "BG"]), na.rm = TRUE)
  
  # BG-BN reach (average BG and BN fluxes)
  reach_daily_flux[d, "BG-BN"] <- mean(c(daily_cd_flux[d, "BG"], daily_cd_flux[d, "BN"]), na.rm = TRUE)
}

# Print results
cat("\n=== REACH DAILY CD FLUX (mg Cd/day) ===\n")
print(round(reach_daily_flux, 2))

cat("\n=== AVERAGE DAILY CD FLUX BY REACH ===\n")
avg_reach_flux <- colMeans(reach_daily_flux, na.rm = TRUE)
print(round(avg_reach_flux, 2))

# Calculate what proportion of daily flux is stored in algae
storage_proportion <- (total_reach_stock / reach_daily_flux) * 100

cat("\n=== ALGAL CD STORAGE AS % OF DAILY FLUX ===\n")
print(round(storage_proportion, 4))

# Calculate average storage proportion by reach
avg_storage_proportion <- colMeans(storage_proportion, na.rm = TRUE)

cat("\n=== AVERAGE STORAGE PROPORTION BY REACH ===\n")
storage_summary <- data.frame(
  Reach = reach_names,
  Avg_Daily_Flux_mg = round(colMeans(reach_daily_flux, na.rm = TRUE), 0),
  Avg_Standing_Stock_mg = round(colMeans(total_reach_stock, na.rm = TRUE), 0),
  Storage_as_Pct_of_Daily_Flux = round(avg_storage_proportion, 4)
)
print(storage_summary)

# Calculate residence time (how many days of flux equals the standing stock)
residence_time_days <- total_reach_stock / reach_daily_flux

cat("\n=== CD RESIDENCE TIME IN ALGAE (days) ===\n")
print(round(residence_time_days, 3))

cat("\n=== AVERAGE RESIDENCE TIME BY REACH ===\n")
avg_residence_time <- colMeans(residence_time_days, na.rm = TRUE)
residence_summary <- data.frame(
  Reach = reach_names,
  Residence_Time_Days = round(avg_residence_time, 3)
)
print(residence_summary)

