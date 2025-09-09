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

# Storage for results
algal_uptake <- matrix(NA, n_dates, n_reaches)
observed_change <- matrix(NA, n_dates, n_reaches)
unexplained_flux <- matrix(NA, n_dates, n_reaches)
reach_names <- paste(site_order[1:n_reaches], site_order[2:(n_reaches+1)], sep="-")

# Process each date
for(d in 1:n_dates) {
  date_data <- data[data$Date == dates[d], ]
  
  # Process each reach
  for(r in 1:n_reaches) {
    upstream_site <- site_order[r]
    downstream_site <- site_order[r+1]
    
    # Get data for this reach
    up_data <- date_data[date_data$Site == upstream_site, ]
    down_data <- date_data[date_data$Site == downstream_site, ]
    
    # Skip if missing data
    if(nrow(up_data) == 0 | nrow(down_data) == 0) next
    
    # Calculate reach length (km)
    reach_length <- site_distances[r+1] - site_distances[r]
    
    # ALGAL CD UPTAKE CALCULATION
    # Average algal biomass and Cd content between sites
    # Using both epilithon (EPIL) and epiphyton (EPIP)
    
    # Epilithon
    epil_biomass_up <- up_data$OM_AREA_g_m2_EPIL[1]
    epil_biomass_down <- down_data$OM_AREA_g_m2_EPIL[1]
    epil_cd_up <- up_data$Cd_EPIL[1]
    epil_cd_down <- down_data$Cd_EPIL[1]
    epil_turnover_up <- up_data$SEP_EPIL[1]
    epil_turnover_down <- down_data$SEP_EPIL[1]
    
    # Epiphyton
    epip_biomass_up <- up_data$OM_AREA_g_m2_EPIP[1]
    epip_biomass_down <- down_data$OM_AREA_g_m2_EPIP[1]
    epip_cd_up <- up_data$Cd_EPIP[1]
    epip_cd_down <- down_data$Cd_EPIP[1]
    epip_turnover_up <- up_data$SEP_EPIP[1]
    epip_turnover_down <- down_data$SEP_EPIP[1]
    
    # Calculate uptake flux (mg Cd/m2/day) for each compartment
    # Using average of upstream and downstream values
    epil_flux <- 0
    if(!is.na(epil_biomass_up) & !is.na(epil_biomass_down)) {
      avg_epil_biomass <- (epil_biomass_up + epil_biomass_down) / 2
      avg_epil_cd <- mean(c(epil_cd_up, epil_cd_down), na.rm=TRUE)
      avg_epil_turnover <- mean(c(epil_turnover_up, epil_turnover_down), na.rm=TRUE)
      if(!is.na(avg_epil_cd) & !is.na(avg_epil_turnover)) {
        epil_flux <- avg_epil_biomass * avg_epil_cd * avg_epil_turnover
      }
    }
    
    epip_flux <- 0
    if(!is.na(epip_biomass_up) & !is.na(epip_biomass_down)) {
      avg_epip_biomass <- (epip_biomass_up + epip_biomass_down) / 2
      avg_epip_cd <- mean(c(epip_cd_up, epip_cd_down), na.rm=TRUE)
      avg_epip_turnover <- mean(c(epip_turnover_up, epip_turnover_down), na.rm=TRUE)
      if(!is.na(avg_epip_cd) & !is.na(avg_epip_turnover)) {
        epip_flux <- avg_epip_biomass * avg_epip_cd * avg_epip_turnover
      }
    }
    
    # Total algal uptake (mg Cd/m2/day)
    total_algal_flux <- epil_flux + epip_flux
    algal_uptake[d, r] <- total_algal_flux
    
    # OBSERVED CD CHANGE CALCULATION
    # Using total Cd (W_Cd) and discharge
    q_up <- up_data$discharge_m3s[1]
    q_down <- down_data$discharge_m3s[1]
    cd_up <- up_data$W_Cd[1]  # Total Cd (mg/L)
    cd_down <- down_data$W_Cd[1]
    
    if(!is.na(q_up) & !is.na(q_down) & !is.na(cd_up) & !is.na(cd_down)) {
      # Calculate loads (mg/s)
      load_up <- cd_up * q_up * 1000  # Convert m3/s to L/s
      load_down <- cd_down * q_down * 1000
      
      # Load change (mg/s)
      load_change <- load_down - load_up
      
      # Convert to per area basis
      # Assume average river width of 20m (adjust as needed)
      river_width <- 20  # meters
      reach_area <- reach_length * 1000 * river_width  # m2
      
      # Convert to mg/m2/day
      observed_change[d, r] <- (load_change * 86400) / reach_area  # 86400 s/day
      
      # Calculate unexplained flux
      unexplained_flux[d, r] <- observed_change[d, r] + total_algal_flux  # + because uptake removes Cd
    }
  }
}

# Create summary table
rownames(algal_uptake) <- dates
colnames(algal_uptake) <- reach_names
rownames(observed_change) <- dates
colnames(observed_change) <- reach_names
rownames(unexplained_flux) <- dates
colnames(unexplained_flux) <- reach_names

# Print results
cat("\n=== ALGAL CD UPTAKE (mg Cd/m2/day) ===\n")
print(round(algal_uptake, 6))

cat("\n=== OBSERVED CD CHANGE (mg Cd/m2/day) ===\n")
print(round(observed_change, 6))

cat("\n=== UNEXPLAINED FLUX (mg Cd/m2/day) ===\n")
cat("Positive = source, Negative = additional sink beyond algae\n")
print(round(unexplained_flux, 6))

# Calculate reach averages
cat("\n=== REACH AVERAGES ===\n")
avg_algal <- colMeans(algal_uptake, na.rm=TRUE)
avg_observed <- colMeans(observed_change, na.rm=TRUE)
avg_unexplained <- colMeans(unexplained_flux, na.rm=TRUE)

reach_summary <- data.frame(
  Reach = reach_names,
  Algal_Uptake = round(avg_algal, 6),
  Observed_Change = round(avg_observed, 6),
  Unexplained_Flux = round(avg_unexplained, 6)
)
print(reach_summary)

# Highlight DL-GR contamination
cat("\n=== DL-GR CONTAMINATION SIGNAL ===\n")
dlgr_index <- which(reach_names == "DL-GR")
cat("Average unexplained flux in DL-GR:", round(avg_unexplained[dlgr_index], 6), "mg Cd/m2/day\n")
cat("This represents contamination input since tributaries assumed clean\n")

# Simple base R plots
par(mfrow=c(2,2))

# Plot 1: Algal uptake by reach
barplot(avg_algal, names.arg=reach_names, 
        main="Average Algal Cd Uptake", 
        ylab="mg Cd/m2/day", xlab="Reach",
        col="lightgreen")

# Plot 2: Observed change by reach
barplot(avg_observed, names.arg=reach_names,
        main="Observed Cd Change", 
        ylab="mg Cd/m2/day", xlab="Reach",
        col=ifelse(avg_observed > 0, "red", "blue"))
abline(h=0)

# Plot 3: Unexplained flux by reach
barplot(avg_unexplained, names.arg=reach_names,
        main="Unexplained Flux (after algal uptake)", 
        ylab="mg Cd/m2/day", xlab="Reach",
        col=ifelse(reach_names == "DL-GR", "orange", "gray"))
abline(h=0)
text(dlgr_index*1.2-0.5, avg_unexplained[dlgr_index], 
     "Contamination", pos=3, cex=0.8)

# Plot 4: Time series of unexplained flux for DL-GR
plot(1:n_dates, unexplained_flux[, dlgr_index], 
     type="b", pch=19,
     main="DL-GR Unexplained Flux Over Time",
     xlab="Sampling Date Index", ylab="mg Cd/m2/day",
     xaxt="n")
axis(1, at=1:n_dates, labels=1:n_dates)
abline(h=0, lty=2)
abline(h=mean(unexplained_flux[, dlgr_index], na.rm=TRUE), col="red", lty=2)
legend("topright", c("Observed", "Mean"), 
       col=c("black", "red"), lty=c(1,2), cex=0.8)

# Reset plot parameters
par(mfrow=c(1,1))

# Additional diagnostics
cat("\n=== QUALITY CHECK ===\n")
cat("Reaches with positive unexplained flux (potential sources):\n")
for(i in 1:length(reach_names)) {
  if(avg_unexplained[i] > 0) {
    cat("  ", reach_names[i], ":", round(avg_unexplained[i], 6), "mg Cd/m2/day\n")
  }
}

cat("\nReaches with negative unexplained flux (additional sinks):\n")
for(i in 1:length(reach_names)) {
  if(avg_unexplained[i] < 0) {
    cat("  ", reach_names[i], ":", round(avg_unexplained[i], 6), "mg Cd/m2/day\n")
  }
}

# Calculate fraction of Cd change explained by algae
cat("\n=== ALGAL CONTRIBUTION ===\n")
for(i in 1:length(reach_names)) {
  if(avg_observed[i] != 0) {
    pct_explained <- abs(avg_algal[i] / avg_observed[i]) * 100
    cat(reach_names[i], ": Algae explain", round(pct_explained, 1), "% of observed change\n")
  }
}