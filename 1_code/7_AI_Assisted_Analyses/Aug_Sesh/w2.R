# Cd flux analysis properly normalized by river distance
# Assumes previous matrices already calculated

# River dimensions
river_width <- c(11, 24, 19, 35, 27, 31)  # WS, DL, GR, GC, BG, BN
site_km <- c(0, 44, 64.8, 89.2, 144.3, 167.8)
reach_length_km <- diff(site_km)
names(reach_length_km) <- paste(site_order[1:n_reaches], site_order[2:n_sites], sep="-")

# Average width for each reach
reach_width <- numeric(n_reaches)
for(r in 1:n_reaches) {
  reach_width[r] <- mean(c(river_width[r], river_width[r+1]))
}

# BIOLOGICAL FLUX NORMALIZED BY RIVER DISTANCE
# Convert to different spatial scales

# 1. Per m² (original)
# biological_flux_total is already in mg/m²/day

# 2. Per meter width × km length (total reach flux)
bio_flux_total_reach <- matrix(NA, n_intervals, n_reaches)
for(t in 1:n_intervals) {
  for(r in 1:n_reaches) {
    if(!is.na(biological_flux_total[t,r])) {
      # mg/m²/day × width(m) × length(km) × 1000(m/km) = mg/day for whole reach
      bio_flux_total_reach[t,r] <- biological_flux_total[t,r] * reach_width[r] * reach_length_km[r] * 1000
    }
  }
}

# 3. Per river km (intensity per unit river distance)
bio_flux_per_km <- matrix(NA, n_intervals, n_reaches)
for(t in 1:n_intervals) {
  for(r in 1:n_reaches) {
    if(!is.na(biological_flux_total[t,r])) {
      # mg/m²/day × width(m) = mg/day per meter of river
      # Then × 1000 = mg/day per km of river
      bio_flux_per_km[t,r] <- biological_flux_total[t,r] * reach_width[r] * 1000
    }
  }
}

rownames(bio_flux_per_km) <- rownames(biological_flux_total)
colnames(bio_flux_per_km) <- colnames(biological_flux_total)

# WATER COLUMN LOAD CHANGES - properly scaled
water_load_change <- matrix(NA, n_intervals, n_reaches)
water_load_per_km <- matrix(NA, n_intervals, n_reaches)

for(t in 1:n_intervals) {
  for(r in 1:n_reaches) {
    # Get average discharge and Cd for the reach
    q_avg_t1 <- mean(c(discharge[t, r], discharge[t, r+1]), na.rm=TRUE)
    q_avg_t2 <- mean(c(discharge[t+1, r], discharge[t+1, r+1]), na.rm=TRUE)
    
    cd_avg_t1 <- mean(c(w_cd[t, r], w_cd[t, r+1]), na.rm=TRUE)
    cd_avg_t2 <- mean(c(w_cd[t+1, r], w_cd[t+1, r+1]), na.rm=TRUE)
    
    if(!is.na(q_avg_t1) & !is.na(q_avg_t2) & !is.na(cd_avg_t1) & !is.na(cd_avg_t2)) {
      # Load change for whole reach (mg/day)
      load_t1 <- cd_avg_t1 * q_avg_t1 * 1000 * 86400  # mg/day
      load_t2 <- cd_avg_t2 * q_avg_t2 * 1000 * 86400  # mg/day
      water_load_change[t,r] <- load_t2 - load_t1
      
      # Load change per km of river
      water_load_per_km[t,r] <- water_load_change[t,r] / reach_length_km[r]
    }
  }
}

# PRINT NORMALIZED COMPARISONS
cat("\n=== REACH CHARACTERISTICS ===\n")
reach_info <- data.frame(
  Reach = names(reach_length_km),
  Length_km = reach_length_km,
  Avg_width_m = reach_width
)
print(reach_info)

cat("\n=== BIOLOGICAL FLUX AT DIFFERENT SCALES ===\n")
cat("\n1. Per unit area (mg/m²/day) - ORIGINAL:\n")
avg_per_area <- colMeans(biological_flux_total, na.rm=TRUE)
print(round(avg_per_area, 6))

cat("\n2. Per kilometer of river (mg/km/day) - NORMALIZED BY DISTANCE:\n")
avg_per_km <- colMeans(bio_flux_per_km, na.rm=TRUE)
print(round(avg_per_km, 2))

cat("\n3. Total for each reach (mg/day) - ACCOUNTS FOR REACH SIZE:\n")
avg_total <- colMeans(bio_flux_total_reach, na.rm=TRUE)
print(round(avg_total, 2))

# COMPARE BIOLOGICAL VS WATER ON PER-KM BASIS
cat("\n=== BIOLOGICAL VS WATER FLUX PER RIVER KM ===\n")
comparison_per_km <- matrix(NA, n_reaches, 5)
rownames(comparison_per_km) <- colnames(biological_flux_total)
colnames(comparison_per_km) <- c("Bio_mg/km/day", "Water_mg/km/day", "Bio_%", "Reach_km", "Width_m")

for(r in 1:n_reaches) {
  bio_mean <- mean(bio_flux_per_km[,r], na.rm=TRUE)
  water_mean <- mean(water_load_per_km[,r], na.rm=TRUE)
  
  comparison_per_km[r, 1] <- bio_mean
  comparison_per_km[r, 2] <- water_mean
  comparison_per_km[r, 3] <- ifelse(water_mean != 0, abs(bio_mean/water_mean) * 100, NA)
  comparison_per_km[r, 4] <- reach_length_km[r]
  comparison_per_km[r, 5] <- reach_width[r]
}

print(round(comparison_per_km, 2))

# IDENTIFY BIOLOGICAL HOTSPOTS (per km basis)
cat("\n=== BIOLOGICAL ACTIVITY INTENSITY (per river km) ===\n")
bio_intensity <- abs(avg_per_km)
ranked <- sort(bio_intensity, decreasing=TRUE)
cat("Ranked by biological activity intensity:\n")
for(i in 1:length(ranked)) {
  reach_name <- names(ranked)[i]
  direction <- ifelse(avg_per_km[reach_name] > 0, "RELEASE", "UPTAKE")
  cat("  ", i, ". ", reach_name, ": ", round(ranked[i], 2), " mg/km/day (", direction, ")\n", sep="")
}

# CALCULATE CUMULATIVE EFFECTS DOWNSTREAM
cat("\n=== CUMULATIVE BIOLOGICAL EFFECTS DOWNSTREAM ===\n")
cumulative_bio <- numeric(n_sites)
names(cumulative_bio) <- site_order

for(s in 2:n_sites) {
  # Sum all biological inputs/outputs upstream of this site
  total_upstream <- 0
  for(r in 1:(s-1)) {
    total_upstream <- total_upstream + mean(bio_flux_total_reach[,r], na.rm=TRUE)
  }
  cumulative_bio[s] <- total_upstream
  cat(site_order[s], " (", site_km[s], " km): ", round(total_upstream, 2), " mg/day cumulative\n", sep="")
}

# TIME SERIES FOR EACH REACH
cat("\n=== TEMPORAL PATTERNS BY REACH (mg/km/day) ===\n")
for(r in 1:n_reaches) {
  cat("\n", colnames(bio_flux_per_km)[r], " (", reach_length_km[r], " km):\n", sep="")
  for(t in 1:n_intervals) {
    if(!is.na(bio_flux_per_km[t,r])) {
      direction <- ifelse(bio_flux_per_km[t,r] > 0, "release", "uptake")
      cat("  ", rownames(bio_flux_per_km)[t], ": ", 
          sprintf("%+.2f", bio_flux_per_km[t,r]), " (", direction, ")\n", sep="")
    }
  }
}

# SPECIAL FOCUS ON DL-GR
cat("\n=== DL-GR DETAILED ANALYSIS ===\n")
dlgr_idx <- which(names(reach_length_km) == "DL-GR")
cat("Reach length: ", reach_length_km[dlgr_idx], " km\n")
cat("Average width: ", reach_width[dlgr_idx], " m\n")
cat("Total area: ", reach_length_km[dlgr_idx] * reach_width[dlgr_idx] / 1000, " km²\n")

cat("\nBiological flux normalized by distance:\n")
for(t in 1:n_intervals) {
  if(!is.na(bio_flux_per_km[t, dlgr_idx])) {
    cat("  ", rownames(bio_flux_per_km)[t], ": ", 
        round(bio_flux_per_km[t, dlgr_idx], 2), " mg/km/day\n", sep="")
  }
}

cat("\nWater load change normalized by distance:\n")
for(t in 1:n_intervals) {
  if(!is.na(water_load_per_km[t, dlgr_idx])) {
    cat("  ", rownames(bio_flux_per_km)[t], ": ", 
        round(water_load_per_km[t, dlgr_idx], 2), " mg/km/day\n", sep="")
  }
}

# CALCULATE BIOLOGICAL DEMAND VS SUPPLY
cat("\n=== BIOLOGICAL DEMAND VS WATER COLUMN SUPPLY ===\n")
for(r in 1:n_reaches) {
  reach_name <- colnames(biological_flux_total)[r]
  
  # For uptake periods (negative flux), what's the demand?
  uptake_periods <- which(biological_flux_total[,r] < 0)
  if(length(uptake_periods) > 0) {
    avg_demand <- mean(abs(biological_flux_total[uptake_periods, r])) * reach_width[r] * 1000  # mg/km/day
    
    # What's the average Cd supply in water?
    avg_cd <- mean(c(w_cd[uptake_periods, r], w_cd[uptake_periods, r+1]), na.rm=TRUE)
    avg_q <- mean(c(discharge[uptake_periods, r], discharge[uptake_periods, r+1]), na.rm=TRUE)
    water_supply <- avg_cd * avg_q * 1000 * 86400 / reach_length_km[r]  # mg/km/day
    
    cat("\n", reach_name, ":\n")
    cat("  Biological demand: ", round(avg_demand, 2), " mg/km/day\n")
    cat("  Water supply: ", round(water_supply, 2), " mg/km/day\n")
    cat("  Demand as % of supply: ", round(avg_demand/water_supply * 100, 2), "%\n")
  }
}

cat("\n=== KEY INSIGHTS ===\n")
cat("1. Biological activity normalized per km shows intensity independent of reach length\n")
cat("2. Longer reaches (WS-DL, BG-BN) have more total flux but may have lower intensity\n")
cat("3. DL-GR is short (", reach_length_km[dlgr_idx], " km) so contamination is concentrated\n")
cat("4. Biological demand is compared to water column supply to assess feasibility\n")