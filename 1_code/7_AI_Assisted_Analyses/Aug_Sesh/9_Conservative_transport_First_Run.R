# Stage 1: Conservative Transport Model
# Cadmium transport in Upper Clark Fork River

# Load data
data <- df2

# Step 1: Data Preparation
# Sort data by Site and Date
data <- data[order(data$Distance_km, data$Date), ]

# Calculate three Cd fractions explicitly
data$dissolved_cd <- data$UF_Cd
data$colloidal_cd <- data$F_Cd - data$UF_Cd
data$suspended_cd <- data$W_Cd - data$F_Cd

# Create site order for downstream progression
site_order <- c("WS", "DL", "GR", "GC", "BG", "BN")
data$site_num <- match(data$Site, site_order)

# Display basic data structure
print("Data structure:")
print(unique(data[,c("Site", "Distance_km", "site_num")]))
print(paste("Date range:", min(data$Date), "to", max(data$Date)))
print(paste("Number of dates:", length(unique(data$Date))))

# Step 2: Calculate Reach-Scale Changes
# Get unique dates for analysis
dates <- unique(data$Date)
print(paste("Analyzing", length(dates), "dates"))

# Create empty dataframe to store reach-scale calculations
reach_data <- data.frame()

for(date in dates) {
  print(paste("Processing date:", date))
  
  # Subset data for this date
  daily_data <- data[data$Date == date, ]
  daily_data <- daily_data[order(daily_data$Distance_km), ]
  
  print(paste("  Sites available:", paste(daily_data$Site, collapse=", ")))
  print(paste("  Number of sites:", nrow(daily_data)))
  
  # Check if we have enough sites for reach analysis
  if(nrow(daily_data) < 2) {
    print(paste("  Skipping date", date, "- insufficient sites"))
    next
  }
  
  # Calculate reach-by-reach changes
  for(i in 2:nrow(daily_data)) {
    upstream_site <- daily_data[i-1, ]
    downstream_site <- daily_data[i, ]
    
    # Check for missing data
    if(nrow(upstream_site) == 0 || nrow(downstream_site) == 0) {
      print(paste("  Warning: Missing site data for reach", i))
      next
    }
    
    # Check for missing key variables
    if(is.na(upstream_site$discharge_m3s) || is.na(downstream_site$discharge_m3s) ||
       is.na(upstream_site$UF_Cd) || is.na(downstream_site$UF_Cd)) {
      print(paste("  Skipping reach due to missing discharge or Cd data"))
      next
    }
    
    # Calculate reach characteristics
    reach_length <- downstream_site$Distance_km - upstream_site$Distance_km
    q_upstream <- upstream_site$discharge_m3s
    q_downstream <- downstream_site$discharge_m3s
    delta_q <- q_downstream - q_upstream
    
    # Store reach data
    reach_row <- data.frame(
      date = date,
      upstream_site = as.character(upstream_site$Site),
      downstream_site = as.character(downstream_site$Site),
      upstream_km = upstream_site$Distance_km,
      downstream_km = downstream_site$Distance_km,
      reach_length_km = reach_length,
      q_upstream = q_upstream,
      q_downstream = q_downstream,
      delta_q = delta_q,
      # Upstream Cd concentrations
      dissolved_up = upstream_site$dissolved_cd,
      colloidal_up = upstream_site$colloidal_cd,
      suspended_up = upstream_site$suspended_cd,
      # Downstream observed Cd concentrations
      dissolved_down_obs = downstream_site$dissolved_cd,
      colloidal_down_obs = downstream_site$colloidal_cd,
      suspended_down_obs = downstream_site$suspended_cd
    )
    
    reach_data <- rbind(reach_data, reach_row)
  }
}

print("Reach-scale data created")
print(paste("Number of reaches analyzed:", nrow(reach_data)))

# Convert Excel serial dates to proper dates if needed
if(all(reach_data$date > 10000)) {
  print("Converting Excel serial dates...")
  reach_data$date_proper <- as.Date(reach_data$date, origin="1900-01-01") - 2
  reach_data$date_char <- format(reach_data$date_proper, "%Y-%m-%d")
  print("Date range:")
  print(paste(min(reach_data$date_char), "to", max(reach_data$date_char)))
}

# Summary of flow changes
print("")
print("=== FLOW ANALYSIS ===")
print(paste("Reaches with flow gains (tributaries):", sum(reach_data$delta_q > 0.1, na.rm=T)))
print(paste("Reaches with flow losses (evap/irrigation):", sum(reach_data$delta_q < -0.1, na.rm=T)))
print(paste("Reaches with minimal flow change:", sum(abs(reach_data$delta_q) <= 0.1, na.rm=T)))

if(sum(reach_data$delta_q > 0.1, na.rm=T) > 0) {
  print(paste("Mean flow gain where delta_q > 0.1:", round(mean(reach_data$delta_q[reach_data$delta_q > 0.1], na.rm=T), 2), "m3/s"))
}
if(sum(reach_data$delta_q < -0.1, na.rm=T) > 0) {
  print(paste("Mean flow loss where delta_q < -0.1:", round(mean(reach_data$delta_q[reach_data$delta_q < -0.1], na.rm=T), 2), "m3/s"))
}

# Step 3: Apply Conservative Mixing Model
# Calculate what downstream concentrations should be under conservative transport

for(i in 1:nrow(reach_data)) {
  # Get reach parameters
  q_up <- reach_data$q_upstream[i]
  q_down <- reach_data$q_downstream[i]
  delta_q <- reach_data$delta_q[i]
  
  # Conservative transport scenarios
  if(abs(delta_q) < 0.1) {
    # No significant flow change - concentrations should remain constant
    reach_data$dissolved_pred[i] <- reach_data$dissolved_up[i]
    reach_data$colloidal_pred[i] <- reach_data$colloidal_up[i]
    reach_data$suspended_pred[i] <- reach_data$suspended_up[i]
    reach_data$tributary_dissolved[i] <- NA
    reach_data$tributary_colloidal[i] <- NA
    reach_data$tributary_suspended[i] <- NA
    reach_data$transport_type[i] <- "No change"
    
  } else if(delta_q > 0) {
    # Tributary input - calculate implied tributary concentrations
    # Working backwards: if observed concentration = conservative mixing, what was Cd_tributary?
    # Cd_down_obs = (Cd_up * Q_up + Cd_trib * delta_Q) / Q_down
    # Therefore: Cd_trib = (Cd_down_obs * Q_down - Cd_up * Q_up) / delta_Q
    
    reach_data$tributary_dissolved[i] <- (reach_data$dissolved_down_obs[i] * q_down - reach_data$dissolved_up[i] * q_up) / delta_q
    reach_data$tributary_colloidal[i] <- (reach_data$colloidal_down_obs[i] * q_down - reach_data$colloidal_up[i] * q_up) / delta_q
    reach_data$tributary_suspended[i] <- (reach_data$suspended_down_obs[i] * q_down - reach_data$suspended_up[i] * q_up) / delta_q
    
    # Calculate predicted concentrations assuming these tributary concentrations
    reach_data$dissolved_pred[i] <- (reach_data$dissolved_up[i] * q_up + reach_data$tributary_dissolved[i] * delta_q) / q_down
    reach_data$colloidal_pred[i] <- (reach_data$colloidal_up[i] * q_up + reach_data$tributary_colloidal[i] * delta_q) / q_down
    reach_data$suspended_pred[i] <- (reach_data$suspended_up[i] * q_up + reach_data$tributary_suspended[i] * delta_q) / q_down
    reach_data$transport_type[i] <- "Tributary input"
    
  } else {
    # Water loss from main stem (evaporation, irrigation, etc.)
    # Assume proportional loss of solutes with water - concentrations remain constant
    reach_data$dissolved_pred[i] <- reach_data$dissolved_up[i]
    reach_data$colloidal_pred[i] <- reach_data$colloidal_up[i]
    reach_data$suspended_pred[i] <- reach_data$suspended_up[i]
    reach_data$tributary_dissolved[i] <- NA  # No tributary, just water loss
    reach_data$tributary_colloidal[i] <- NA
    reach_data$tributary_suspended[i] <- NA
    reach_data$transport_type[i] <- "Water loss"
  }
}

# Step 4: Calculate Deviations from Conservative Transport
reach_data$dissolved_residual <- reach_data$dissolved_down_obs - reach_data$dissolved_pred
reach_data$colloidal_residual <- reach_data$colloidal_down_obs - reach_data$colloidal_pred  
reach_data$suspended_residual <- reach_data$suspended_down_obs - reach_data$suspended_pred

# Calculate relative deviations (%)
reach_data$dissolved_rel_dev <- (reach_data$dissolved_residual / reach_data$dissolved_pred) * 100
reach_data$colloidal_rel_dev <- (reach_data$colloidal_residual / reach_data$colloidal_pred) * 100
reach_data$suspended_rel_dev <- (reach_data$suspended_residual / reach_data$suspended_pred) * 100

# Step 5: Summary Statistics
print("=== CONSERVATIVE TRANSPORT MODEL RESULTS ===")
print("")

# Categorize reaches by transport type
tributary_reaches <- reach_data[reach_data$delta_q > 0.1 & !is.na(reach_data$delta_q), ]
water_loss_reaches <- reach_data[reach_data$delta_q < -0.1 & !is.na(reach_data$delta_q), ]
stable_reaches <- reach_data[abs(reach_data$delta_q) <= 0.1 & !is.na(reach_data$delta_q), ]

print(paste("Reaches with tributary inputs:", nrow(tributary_reaches)))
print(paste("Reaches with water losses:", nrow(water_loss_reaches)))
print(paste("Reaches with stable flow:", nrow(stable_reaches)))

if(nrow(tributary_reaches) > 0) {
  print("")
  print("Mean Tributary Cd Concentrations (mg/L) - where tributaries exist:")
  print(paste("Dissolved:", round(mean(tributary_reaches$tributary_dissolved, na.rm=T), 6)))
  print(paste("Colloidal:", round(mean(tributary_reaches$tributary_colloidal, na.rm=T), 6))) 
  print(paste("Suspended:", round(mean(tributary_reaches$tributary_suspended, na.rm=T), 6)))
  
  # Check for negative tributary concentrations (indicating net removal)
  neg_trib <- sum(tributary_reaches$tributary_dissolved < 0, na.rm=T)
  if(neg_trib > 0) {
    print(paste("NOTE:", neg_trib, "tributary reaches show negative implied Cd concentrations"))
    print("This suggests net removal processes (biological uptake, precipitation, etc.)")
  }
}

print("")
print("Mean Absolute Deviations from Conservative Transport:")
print(paste("All reaches - Dissolved:", round(mean(abs(reach_data$dissolved_residual), na.rm=T), 6), "mg/L"))
print(paste("All reaches - Colloidal:", round(mean(abs(reach_data$colloidal_residual), na.rm=T), 6), "mg/L"))
print(paste("All reaches - Suspended:", round(mean(abs(reach_data$suspended_residual), na.rm=T), 6), "mg/L"))

if(nrow(tributary_reaches) > 0) {
  print("")
  print("Tributary reaches only:")
  print(paste("Dissolved deviation:", round(mean(abs(tributary_reaches$dissolved_residual), na.rm=T), 6), "mg/L"))
}

if(nrow(water_loss_reaches) > 0) {
  print("")
  print("Water loss reaches only:")
  print(paste("Dissolved deviation:", round(mean(abs(water_loss_reaches$dissolved_residual), na.rm=T), 6), "mg/L"))
  print("(Large deviations here indicate non-proportional solute loss or other processes)")
}

print("")
print("Mean Relative Deviations (%) - excluding NA values:")
dissolved_rel_clean <- reach_data$dissolved_rel_dev[!is.na(reach_data$dissolved_rel_dev) & is.finite(reach_data$dissolved_rel_dev)]
colloidal_rel_clean <- reach_data$colloidal_rel_dev[!is.na(reach_data$colloidal_rel_dev) & is.finite(reach_data$colloidal_rel_dev)]
suspended_rel_clean <- reach_data$suspended_rel_dev[!is.na(reach_data$suspended_rel_dev) & is.finite(reach_data$suspended_rel_dev)]

if(length(dissolved_rel_clean) > 0) print(paste("Dissolved:", round(mean(abs(dissolved_rel_clean)), 1), "%"))
if(length(colloidal_rel_clean) > 0) print(paste("Colloidal:", round(mean(abs(colloidal_rel_clean)), 1), "%"))
if(length(suspended_rel_clean) > 0) print(paste("Suspended:", round(mean(abs(suspended_rel_clean)), 1), "%"))

# Step 6: Visualize Patterns
par(mfrow=c(2,2))

# Plot 1: Tributary Cd concentrations by reach (only where tributaries exist)
trib_reaches <- reach_data[reach_data$delta_q > 0.1 & !is.na(reach_data$tributary_dissolved), ]
if(nrow(trib_reaches) > 0) {
  y_max <- max(c(trib_reaches$tributary_dissolved, trib_reaches$tributary_colloidal, trib_reaches$tributary_suspended), na.rm=T)
  plot(trib_reaches$upstream_km, trib_reaches$tributary_dissolved, 
       pch=16, col="blue", ylim=c(0, y_max),
       xlab="Upstream Distance (km)", ylab="Tributary Cd (mg/L)", 
       main="Estimated Tributary Cd Concentrations")
  points(trib_reaches$upstream_km, trib_reaches$tributary_colloidal, pch=16, col="green")
  points(trib_reaches$upstream_km, trib_reaches$tributary_suspended, pch=16, col="red")
  legend("topright", c("Dissolved", "Colloidal", "Suspended"), col=c("blue", "green", "red"), pch=16)
} else {
  plot(1, type="n", xlab="", ylab="", main="No significant tributary inputs detected")
}

# Plot 2: Observed vs Predicted (Dissolved) - exclude NA values
valid_data <- reach_data[!is.na(reach_data$dissolved_pred) & !is.na(reach_data$dissolved_down_obs), ]
if(nrow(valid_data) > 0) {
  plot(valid_data$dissolved_pred, valid_data$dissolved_down_obs,
       pch=16, col="blue", xlab="Predicted Dissolved Cd", ylab="Observed Dissolved Cd",
       main="Conservative Transport: Dissolved Cd")
  abline(0,1, lty=2)
  
  # Add R-squared
  r_sq <- cor(valid_data$dissolved_pred, valid_data$dissolved_down_obs)^2
  text(min(valid_data$dissolved_pred), max(valid_data$dissolved_down_obs), 
       paste("R² =", round(r_sq, 3)), pos=4)
}

# Plot 3: Residuals by reach - exclude NA values
valid_resid <- reach_data[!is.na(reach_data$dissolved_residual), ]
if(nrow(valid_resid) > 0) {
  plot(valid_resid$upstream_km, valid_resid$dissolved_residual,
       pch=16, col="blue", xlab="Upstream Distance (km)", ylab="Residual (Obs - Pred)",
       main="Deviations from Conservative Transport")
  
  # Add other fractions if available
  if(any(!is.na(reach_data$colloidal_residual))) {
    points(valid_resid$upstream_km, valid_resid$colloidal_residual, pch=16, col="green")
  }
  if(any(!is.na(reach_data$suspended_residual))) {
    points(valid_resid$upstream_km, valid_resid$suspended_residual, pch=16, col="red")
  }
  
  abline(h=0, lty=2)
  legend("topright", c("Dissolved", "Colloidal", "Suspended"), col=c("blue", "green", "red"), pch=16)
}

# Plot 4: Temporal patterns in deviations
if(exists("date_char", where=reach_data)) {
  # Use converted dates
  boxplot(reach_data$dissolved_rel_dev ~ reach_data$date_char, 
          main="Dissolved Cd Deviations by Date", xlab="Date", ylab="Relative Deviation (%)",
          las=2)
} else {
  # Use original date format
  boxplot(reach_data$dissolved_rel_dev ~ reach_data$date, 
          main="Dissolved Cd Deviations by Date", xlab="Date", ylab="Relative Deviation (%)",
          las=2)
}

# Display key results table
print("")
print("=== REACH-BY-REACH ANALYSIS ===")
summary_table <- reach_data[,c("date", "upstream_site", "downstream_site", "delta_q", "transport_type",
                               "dissolved_down_obs", "dissolved_pred", "dissolved_residual",
                               "tributary_dissolved")]
print(head(summary_table, 15))

# Save results
write.csv(reach_data, "conservative_transport_results.csv", row.names=FALSE)
print("")
print("Results saved to 'conservative_transport_results.csv'")