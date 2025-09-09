# Fraction-Specific Cadmium Analysis
# Analyzing dissolved, colloidal, and suspended fractions separately

# Assumes df is loaded and transport_results exists from previous analysis

# === CALCULATE FRACTION CONCENTRATIONS ===
# These are already in the data but let's verify
df$dissolved_Cd <- df$W_Cd * df$Dissolved_pct
df$colloidal_Cd <- df$W_Cd * df$Colloidal_pct  
df$suspended_Cd <- df$W_Cd * df$Suspended_pct

# Calculate loads for each fraction
df$dissolved_load <- df$dissolved_Cd * df$discharge_m3s
df$colloidal_load <- df$colloidal_Cd * df$discharge_m3s
df$suspended_load <- df$suspended_Cd * df$discharge_m3s

# Get minimum concentrations for each fraction (needed for mass balance)
min_dissolved <- min(df$dissolved_Cd)
min_colloidal <- min(df$colloidal_Cd)
min_suspended <- min(df$suspended_Cd)

cat("=== MINIMUM CONCENTRATIONS BY FRACTION ===\n")
cat("Min dissolved:", round(min_dissolved, 6), "\n")
cat("Min colloidal:", round(min_colloidal, 6), "\n")
cat("Min suspended:", round(min_suspended, 6), "\n\n")

# === THREE-PANEL DISCHARGE-CONCENTRATION PLOT ===
par(mfrow = c(1,3))

# Panel 1: Dissolved
plot(df$discharge_m3s, df$dissolved_Cd,
     log = "xy",
     main = "Dissolved Cd",
     xlab = "Discharge (m³/s)",
     ylab = "Dissolved Cd (µg/L)",
     pch = 16, col = as.numeric(factor(df$Site)))
legend("topright", legend = unique(df$Site), 
       col = 1:length(unique(df$Site)), pch = 16, cex = 0.6)

# Panel 2: Colloidal  
plot(df$discharge_m3s, df$colloidal_Cd,
     log = "xy",
     main = "Colloidal Cd",
     xlab = "Discharge (m³/s)",
     ylab = "Colloidal Cd (µg/L)",
     pch = 16, col = as.numeric(factor(df$Site)))

# Panel 3: Suspended
plot(df$discharge_m3s, df$suspended_Cd,
     log = "xy", 
     main = "Suspended Cd",
     xlab = "Discharge (m³/s)",
     ylab = "Suspended Cd (µg/L)",
     pch = 16, col = as.numeric(factor(df$Site)))

# === FRACTION-SPECIFIC DILUTION FROM MASS BALANCE ===
cat("\n=== FRACTION-SPECIFIC DILUTION EFFECTS ===\n")

# For each fraction, calculate how much concentration change is due to dilution
fraction_dilution <- list()
counter <- 1

for(date in dates) {
  date_data <- df[df$Date == date,]
  
  for(i in 1:nrow(transport_results$actual_segments)) {
    upstream <- date_data[date_data$Site == transport_results$actual_segments$upstream_site[i],]
    downstream <- date_data[date_data$Site == transport_results$actual_segments$downstream_site[i],]
    
    if(nrow(upstream) == 0 | nrow(downstream) == 0) next
    
    Q_change <- downstream$discharge_m3s - upstream$discharge_m3s
    
    for(frac in c("dissolved", "colloidal", "suspended")) {
      conc_col <- paste0(frac, "_Cd")
      
      # Calculate expected concentration from dilution alone
      if(Q_change > 0) {
        # Gaining reach - dilution with clean water
        min_conc <- get(paste0("min_", frac))
        tributary_conc <- min_conc
        expected_conc <- (upstream[[conc_col]] * upstream$discharge_m3s + 
                            tributary_conc * Q_change) / downstream$discharge_m3s
      } else {
        # Losing reach - concentration increases
        expected_conc <- upstream[[conc_col]]  # No dilution in losing reach
      }
      
      actual_conc <- downstream[[conc_col]]
      
      # Dilution vs process effects
      dilution_change <- expected_conc - upstream[[conc_col]]
      process_change <- actual_conc - expected_conc
      
      fraction_dilution[[counter]] <- data.frame(
        Date = date,
        reach = transport_results$actual_segments$reach_name[i],
        fraction = frac,
        upstream_conc = upstream[[conc_col]],
        expected_conc = expected_conc,
        actual_conc = actual_conc,
        dilution_change = dilution_change,
        process_change = process_change,
        dilution_pct = (dilution_change / upstream[[conc_col]]) * 100,
        process_pct = (process_change / expected_conc) * 100
      )
      counter <- counter + 1
    }
  }
}

fraction_dilution <- do.call(rbind, fraction_dilution)

# Summary by fraction
dilution_by_fraction <- fraction_dilution %>%
  group_by(fraction) %>%
  summarise(
    mean_dilution_pct = mean(dilution_pct, na.rm = TRUE),
    mean_process_pct = mean(process_pct, na.rm = TRUE),
    sd_process_pct = sd(process_pct, na.rm = TRUE)
  )

cat("Dilution % = concentration change from physical mixing\n")
cat("Process % = concentration change from biological/geochemical processes\n\n")
print(dilution_by_fraction)

cat("\nInterpretation:\n")
most_affected <- dilution_by_fraction$fraction[which.max(abs(dilution_by_fraction$mean_process_pct))]
cat("- Fraction most affected by non-dilution processes:", most_affected, "\n")

if(any(dilution_by_fraction$mean_process_pct > 0)) {
  sources <- dilution_by_fraction$fraction[dilution_by_fraction$mean_process_pct > 0]
  cat("- Fractions showing net sources:", paste(sources, collapse = ", "), "\n")
}

if(any(dilution_by_fraction$mean_process_pct < 0)) {
  sinks <- dilution_by_fraction$fraction[dilution_by_fraction$mean_process_pct < 0]
  cat("- Fractions showing net sinks:", paste(sinks, collapse = ", "), "\n")
}

# === FRACTION MASS BALANCE BY REACH ===
# Get unique dates for analysis
dates <- unique(df$Date)

# Quick analysis of which fraction dominates residuals
fraction_residuals <- list()
counter <- 1

for(date in dates) {
  date_data <- df[df$Date == date,]
  
  for(i in 1:nrow(transport_results$actual_segments)) {
    upstream <- date_data[date_data$Site == transport_results$actual_segments$upstream_site[i],]
    downstream <- date_data[date_data$Site == transport_results$actual_segments$downstream_site[i],]
    
    if(nrow(upstream) == 0 | nrow(downstream) == 0) next
    
    Q_change <- downstream$discharge_m3s - upstream$discharge_m3s
    
    # Calculate for each fraction
    for(frac in c("dissolved", "colloidal", "suspended")) {
      load_col <- paste0(frac, "_load")
      conc_col <- paste0(frac, "_Cd")
      min_conc <- get(paste0("min_", frac))
      
      if(Q_change > 0) {
        expected <- upstream[[load_col]] + (Q_change * min_conc)
      } else {
        expected <- upstream[[load_col]] * (downstream$discharge_m3s / upstream$discharge_m3s)
      }
      
      actual <- downstream[[load_col]]
      residual_pct <- ((actual - expected) / expected) * 100
      
      fraction_residuals[[counter]] <- data.frame(
        Date = date,
        reach = transport_results$actual_segments$reach_name[i],
        fraction = frac,
        residual_pct = residual_pct
      )
      counter <- counter + 1
    }
  }
}

fraction_residuals <- do.call(rbind, fraction_residuals)

# Summary by fraction
fraction_summary <- fraction_residuals %>%
  group_by(fraction) %>%
  summarise(
    mean_residual = mean(residual_pct, na.rm = TRUE),
    sd_residual = sd(residual_pct, na.rm = TRUE),
    max_residual = max(residual_pct, na.rm = TRUE),
    min_residual = min(residual_pct, na.rm = TRUE)
  )

cat("\n=== RESIDUALS BY FRACTION ===\n")
print(fraction_summary)

# Which fraction shows largest deviations?
cat("\n=== FRACTION WITH LARGEST DEVIATIONS ===\n")
most_variable <- fraction_summary$fraction[which.max(fraction_summary$sd_residual)]
cat("Most variable fraction:", most_variable, "\n")
cat("This suggests", most_variable, "Cd is most affected by non-transport processes\n")

# === HOTSPOTS BY FRACTION ===
fraction_hotspots <- fraction_residuals %>%
  group_by(reach, fraction) %>%
  summarise(mean_residual = mean(residual_pct, na.rm = TRUE)) %>%
  filter(abs(mean_residual) > 20) %>%
  arrange(desc(abs(mean_residual)))

if(nrow(fraction_hotspots) > 0) {
  cat("\n=== FRACTION-SPECIFIC HOTSPOTS (>20% deviation) ===\n")
  print(fraction_hotspots)
}

# === STORE FRACTION RESULTS ===
fraction_results <- list(
  fraction_dilution = fraction_dilution,
  dilution_by_fraction = dilution_by_fraction,
  fraction_residuals = fraction_residuals,
  fraction_summary = fraction_summary,
  fraction_hotspots = fraction_hotspots
)

cat("\n=== FRACTION ANALYSIS COMPLETE ===\n")
cat("Results stored in 'fraction_results'\n")
