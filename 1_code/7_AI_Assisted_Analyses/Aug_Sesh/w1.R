# Simple R code for biological Cd flux analysis
# Using the integrated dataset with replicates

# Read the integrated dataset
data <- df1

# Convert date to Date format
data$SAMPLING_DATE <- as.Date(data$SAMPLING_DATE)

# Define site order and distances
site_order <- c("WS", "DL", "GR", "GC", "BG", "BN")
site_km <- c(0, 44, 64.8, 89.2, 144.3, 167.8)

# Get unique dates and sites
dates <- sort(unique(data$SAMPLING_DATE))
n_dates <- length(dates)
n_sites <- length(site_order)
n_reaches <- n_sites - 1

# Initialize storage matrices
epil_biomass <- matrix(NA, n_dates, n_sites)
epip_biomass <- matrix(NA, n_dates, n_sites)
fila_biomass <- matrix(NA, n_dates, n_sites)

epil_cd <- matrix(NA, n_dates, n_sites)
epip_cd <- matrix(NA, n_dates, n_sites)
fila_cd <- matrix(NA, n_dates, n_sites)

epil_sep <- matrix(NA, n_dates, n_sites)
epip_sep <- matrix(NA, n_dates, n_sites)
fila_sep <- matrix(NA, n_dates, n_sites)

# Water chemistry matrices
uf_cd <- matrix(NA, n_dates, n_sites)
f_cd <- matrix(NA, n_dates, n_sites)
w_cd <- matrix(NA, n_dates, n_sites)
discharge <- matrix(NA, n_dates, n_sites)

# Fill matrices - average replicates for each site/date/compartment
for(d in 1:n_dates) {
  for(s in 1:n_sites) {
    # Subset for this date and site
    subset <- data[data$SAMPLING_DATE == dates[d] & data$SITE == site_order[s], ]
    
    if(nrow(subset) > 0) {
      # Water chemistry (same for all rows at this site/date)
      uf_cd[d,s] <- subset$UF_Cd[1]
      f_cd[d,s] <- subset$F_Cd[1]
      w_cd[d,s] <- subset$W_Cd[1]
      discharge[d,s] <- subset$discharge_m3s[1]
      
      # EPIL
      epil_rows <- subset[subset$SAMPLE_DESCRIPTOR == "EPIL", ]
      if(nrow(epil_rows) > 0) {
        epil_biomass[d,s] <- mean(epil_rows$OM.AREA.g.m2, na.rm=TRUE)
        epil_cd[d,s] <- mean(epil_rows$Cd, na.rm=TRUE)
        epil_sep[d,s] <- mean(epil_rows$SEP, na.rm=TRUE)
      }
      
      # EPIP (exclude outlier SEP > 100)
      epip_rows <- subset[subset$SAMPLE_DESCRIPTOR == "EPIP", ]
      epip_rows <- epip_rows[!is.na(epip_rows$SEP) & epip_rows$SEP < 100, ]
      if(nrow(epip_rows) > 0) {
        epip_biomass[d,s] <- mean(epip_rows$OM.AREA.g.m2, na.rm=TRUE)
        epip_cd[d,s] <- mean(epip_rows$Cd, na.rm=TRUE)
        epip_sep[d,s] <- mean(epip_rows$SEP, na.rm=TRUE)
      }
      
      # FILA
      fila_rows <- subset[subset$SAMPLE_DESCRIPTOR == "FILA", ]
      if(nrow(fila_rows) > 0) {
        fila_biomass[d,s] <- mean(fila_rows$OM.AREA.g.m2, na.rm=TRUE)
        fila_cd[d,s] <- mean(fila_rows$Cd, na.rm=TRUE)
        fila_sep[d,s] <- mean(fila_rows$SEP, na.rm=TRUE)
      }
    }
  }
}

# Calculate biological flux for each reach and time interval
n_intervals <- n_dates - 1
biological_flux_total <- matrix(NA, n_intervals, n_reaches)
biological_flux_epil <- matrix(NA, n_intervals, n_reaches)
biological_flux_epip <- matrix(NA, n_intervals, n_reaches)
biological_flux_fila <- matrix(NA, n_intervals, n_reaches)

rownames(biological_flux_total) <- paste(dates[1:n_intervals], "to", dates[2:n_dates])
colnames(biological_flux_total) <- paste(site_order[1:n_reaches], site_order[2:n_sites], sep="-")

# Process each time interval
for(t in 1:n_intervals) {
  days_between <- as.numeric(dates[t+1] - dates[t])
  
  for(r in 1:n_reaches) {
    site_up <- r
    site_down <- r + 1
    
    # EPILITHON
    # Average biomass across reach at each time
    epil_bio_t1 <- mean(c(epil_biomass[t, site_up], epil_biomass[t, site_down]), na.rm=TRUE)
    epil_bio_t2 <- mean(c(epil_biomass[t+1, site_up], epil_biomass[t+1, site_down]), na.rm=TRUE)
    
    if(!is.na(epil_bio_t1) & !is.na(epil_bio_t2)) {
      epil_delta <- epil_bio_t2 - epil_bio_t1
      # Average tissue Cd concentration
      avg_cd <- mean(c(epil_cd[t, site_up], epil_cd[t, site_down],
                       epil_cd[t+1, site_up], epil_cd[t+1, site_down]), na.rm=TRUE)
      
      if(!is.na(avg_cd)) {
        # Flux: negative = uptake from water, positive = release to water
        biological_flux_epil[t,r] <- -epil_delta * avg_cd / days_between
      }
    }
    
    # EPIPHYTON
    epip_bio_t1 <- mean(c(epip_biomass[t, site_up], epip_biomass[t, site_down]), na.rm=TRUE)
    epip_bio_t2 <- mean(c(epip_biomass[t+1, site_up], epip_biomass[t+1, site_down]), na.rm=TRUE)
    
    if(!is.na(epip_bio_t1) & !is.na(epip_bio_t2)) {
      epip_delta <- epip_bio_t2 - epip_bio_t1
      avg_cd <- mean(c(epip_cd[t, site_up], epip_cd[t, site_down],
                       epip_cd[t+1, site_up], epip_cd[t+1, site_down]), na.rm=TRUE)
      
      if(!is.na(avg_cd)) {
        biological_flux_epip[t,r] <- -epip_delta * avg_cd / days_between
      }
    }
    
    # FILAMENTOUS
    fila_bio_t1 <- mean(c(fila_biomass[t, site_up], fila_biomass[t, site_down]), na.rm=TRUE)
    fila_bio_t2 <- mean(c(fila_biomass[t+1, site_up], fila_biomass[t+1, site_down]), na.rm=TRUE)
    
    if(!is.na(fila_bio_t1) & !is.na(fila_bio_t2)) {
      fila_delta <- fila_bio_t2 - fila_bio_t1
      avg_cd <- mean(c(fila_cd[t, site_up], fila_cd[t, site_down],
                       fila_cd[t+1, site_up], fila_cd[t+1, site_down]), na.rm=TRUE)
      
      if(!is.na(avg_cd)) {
        biological_flux_fila[t,r] <- -fila_delta * avg_cd / days_between
      }
    }
    
    # Total biological flux
    biological_flux_total[t,r] <- sum(c(biological_flux_epil[t,r], 
                                        biological_flux_epip[t,r], 
                                        biological_flux_fila[t,r]), na.rm=TRUE)
  }
}

# Calculate observed water column changes
water_cd_change <- matrix(NA, n_intervals, n_reaches)
dissolved_change <- matrix(NA, n_intervals, n_reaches)
colloidal_change <- matrix(NA, n_intervals, n_reaches)
suspended_change <- matrix(NA, n_intervals, n_reaches)
load_change <- matrix(NA, n_intervals, n_reaches)

for(t in 1:n_intervals) {
  for(r in 1:n_reaches) {
    # Average across reach at each time
    w_cd_t1 <- mean(c(w_cd[t, r], w_cd[t, r+1]), na.rm=TRUE)
    w_cd_t2 <- mean(c(w_cd[t+1, r], w_cd[t+1, r+1]), na.rm=TRUE)
    water_cd_change[t,r] <- w_cd_t2 - w_cd_t1
    
    # Dissolved
    uf_t1 <- mean(c(uf_cd[t, r], uf_cd[t, r+1]), na.rm=TRUE)
    uf_t2 <- mean(c(uf_cd[t+1, r], uf_cd[t+1, r+1]), na.rm=TRUE)
    dissolved_change[t,r] <- uf_t2 - uf_t1
    
    # Colloidal (F_Cd - UF_Cd)
    coll_t1 <- mean(c(f_cd[t, r] - uf_cd[t, r], 
                      f_cd[t, r+1] - uf_cd[t, r+1]), na.rm=TRUE)
    coll_t2 <- mean(c(f_cd[t+1, r] - uf_cd[t+1, r], 
                      f_cd[t+1, r+1] - uf_cd[t+1, r+1]), na.rm=TRUE)
    colloidal_change[t,r] <- coll_t2 - coll_t1
    
    # Suspended (W_Cd - F_Cd)
    susp_t1 <- mean(c(w_cd[t, r] - f_cd[t, r],
                      w_cd[t, r+1] - f_cd[t, r+1]), na.rm=TRUE)
    susp_t2 <- mean(c(w_cd[t+1, r] - f_cd[t+1, r],
                      w_cd[t+1, r+1] - f_cd[t+1, r+1]), na.rm=TRUE)
    suspended_change[t,r] <- susp_t2 - susp_t1
    
    # Load change (mg/s)
    q_avg_t1 <- mean(c(discharge[t, r], discharge[t, r+1]), na.rm=TRUE)
    q_avg_t2 <- mean(c(discharge[t+1, r], discharge[t+1, r+1]), na.rm=TRUE)
    load_t1 <- w_cd_t1 * q_avg_t1 * 1000  # Convert to mg/s
    load_t2 <- w_cd_t2 * q_avg_t2 * 1000
    load_change[t,r] <- load_t2 - load_t1
  }
}

# Print results
cat("\n=== BIOLOGICAL CD FLUX BY COMPARTMENT (mg/m2/day) ===\n")
cat("Negative = uptake from water, Positive = release to water\n\n")

cat("EPILITHON:\n")
print(round(biological_flux_epil, 6))

cat("\nEPIPHYTON:\n")
print(round(biological_flux_epip, 6))

cat("\nFILAMENTOUS:\n")
print(round(biological_flux_fila, 6))

cat("\nTOTAL BIOLOGICAL FLUX:\n")
print(round(biological_flux_total, 6))

cat("\n=== OBSERVED WATER COLUMN CHANGES (mg/L) ===\n")
cat("\nTotal Cd:\n")
print(round(water_cd_change, 6))

cat("\nDissolved (UF_Cd):\n")
print(round(dissolved_change, 6))

cat("\nColloidal (F_Cd - UF_Cd):\n")
print(round(colloidal_change, 6))

cat("\nSuspended (W_Cd - F_Cd):\n")
print(round(suspended_change, 6))

# Analyze biological release events
cat("\n=== BIOLOGICAL RELEASE ANALYSIS ===\n")
release_events <- which(biological_flux_total > 0, arr.ind=TRUE)
if(nrow(release_events) > 0) {
  for(i in 1:nrow(release_events)) {
    t_idx <- release_events[i,1]
    r_idx <- release_events[i,2]
    
    cat("\n", rownames(biological_flux_total)[t_idx], " at ", 
        colnames(biological_flux_total)[r_idx], ":\n", sep="")
    cat("  Total release: ", round(biological_flux_total[t_idx, r_idx], 6), " mg/m2/day\n")
    cat("  Water column response:\n")
    cat("    Dissolved: ", round(dissolved_change[t_idx, r_idx], 6), " mg/L\n")
    cat("    Colloidal: ", round(colloidal_change[t_idx, r_idx], 6), " mg/L\n")
    cat("    Suspended: ", round(suspended_change[t_idx, r_idx], 6), " mg/L\n")
  }
}

# Check compartment responsiveness to dissolved Cd
cat("\n=== COMPARTMENT RESPONSIVENESS TO DISSOLVED CD ===\n")
epil_cor <- cor(as.vector(epil_cd), as.vector(uf_cd), use="complete.obs")
epip_cor <- cor(as.vector(epip_cd), as.vector(uf_cd), use="complete.obs")
fila_cor <- cor(as.vector(fila_cd), as.vector(uf_cd), use="complete.obs")

cat("Correlation between tissue Cd and dissolved Cd (UF_Cd):\n")
cat("  EPIL: r = ", round(epil_cor, 3), "\n")
cat("  EPIP: r = ", round(epip_cor, 3), "\n")
cat("  FILA: r = ", round(fila_cor, 3), "\n")

# Focus on DL-GR reach
cat("\n=== DL-GR CONTAMINATION ANALYSIS ===\n")
dlgr_idx <- which(colnames(biological_flux_total) == "DL-GR")

cat("Biological flux in DL-GR across all periods:\n")
dlgr_bio <- biological_flux_total[, dlgr_idx]
dlgr_bio <- dlgr_bio[!is.na(dlgr_bio)]
cat("  Mean: ", round(mean(dlgr_bio), 6), " mg/m2/day\n")
cat("  Range: ", round(min(dlgr_bio), 6), " to ", round(max(dlgr_bio), 6), "\n")

cat("\nWater Cd change in DL-GR:\n")
dlgr_water <- water_cd_change[, dlgr_idx]
dlgr_water <- dlgr_water[!is.na(dlgr_water)]
cat("  Mean: ", round(mean(dlgr_water), 6), " mg/L\n")
cat("  Range: ", round(min(dlgr_water), 6), " to ", round(max(dlgr_water), 6), "\n")

cat("\nLoad change in DL-GR:\n")
dlgr_load <- load_change[, dlgr_idx]
dlgr_load <- dlgr_load[!is.na(dlgr_load)]
cat("  Mean: ", round(mean(dlgr_load), 3), " mg/s\n")

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Mean biological flux by reach (mg/m2/day):\n")
for(r in 1:n_reaches) {
  mean_flux <- mean(biological_flux_total[,r], na.rm=TRUE)
  cat("  ", colnames(biological_flux_total)[r], ": ", round(mean_flux, 6))
  if(mean_flux > 0) cat(" (net release)") else cat(" (net uptake)")
  cat("\n")
}

cat("\nContribution by compartment (% of total flux magnitude):\n")
total_abs <- sum(abs(biological_flux_total), na.rm=TRUE)
epil_pct <- sum(abs(biological_flux_epil), na.rm=TRUE) / total_abs * 100
epip_pct <- sum(abs(biological_flux_epip), na.rm=TRUE) / total_abs * 100
fila_pct <- sum(abs(biological_flux_fila), na.rm=TRUE) / total_abs * 100
cat("  EPIL: ", round(epil_pct, 1), "%\n")
cat("  EPIP: ", round(epip_pct, 1), "%\n")
cat("  FILA: ", round(fila_pct, 1), "%\n")

# Calculate unexplained flux (observed - biological)
cat("\n=== UNEXPLAINED FLUX ===\n")
cat("(Positive = additional sources, Negative = additional sinks)\n\n")
# Note: This is simplified - proper units conversion needed
unexplained <- water_cd_change * 1000  # Rough conversion
for(r in 1:n_reaches) {
  cat(colnames(biological_flux_total)[r], ":\n")
  for(t in 1:n_intervals) {
    if(!is.na(water_cd_change[t,r]) & !is.na(biological_flux_total[t,r])) {
      unexp <- water_cd_change[t,r] - (biological_flux_total[t,r] / 1000)  # Simplified
      if(abs(unexp) > 0.001) {
        cat("  ", rownames(biological_flux_total)[t], ": ", 
            round(unexp, 6), " mg/L\n", sep="")
      }
    }
  }
}

