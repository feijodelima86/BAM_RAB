# Cadmium Mass Balance Analysis - Simple Transport Model
# Direct implementation for THIS specific dataset only

library(tidyverse)

# Read data - assuming it's loaded as 'df'
# For this example, you'll need to load your data as 'df'

# === SETUP ===
# Calculate loads directly
df$W_Cd_load <- df$W_Cd * df$discharge_m3s

# Get minimum Cd for tributary assumption
min_Cd <- min(df$W_Cd)
cat("Minimum Cd concentration for tributary water:", min_Cd, "\n\n")

# Get site order from upstream to downstream
sites_ordered <- unique(df[order(df$Distance_km), c("Site", "Distance_km")])
print(sites_ordered)

# Define actual river segments (no skipping)
actual_segments <- data.frame(
  upstream_site = sites_ordered$Site[-nrow(sites_ordered)],
  downstream_site = sites_ordered$Site[-1],
  upstream_km = sites_ordered$Distance_km[-nrow(sites_ordered)],
  downstream_km = sites_ordered$Distance_km[-1]
)
actual_segments$reach_length <- actual_segments$downstream_km - actual_segments$upstream_km
actual_segments$reach_name <- paste(actual_segments$upstream_site, 
                                    actual_segments$downstream_site, sep = "-")

cat("\n=== RIVER SEGMENTS ===\n")
print(actual_segments)

# === MASS BALANCE CALCULATIONS ===
# Get unique dates for analysis
dates <- unique(df$Date)

# Store results
load_analysis <- list()
counter <- 1

cat("\n=== MASS BALANCE BY DATE AND REACH ===\n")

for(date in dates) {
  date_data <- df[df$Date == date,]
  
  for(i in 1:nrow(actual_segments)) {
    # Get upstream and downstream data
    upstream <- date_data[date_data$Site == actual_segments$upstream_site[i],]
    downstream <- date_data[date_data$Site == actual_segments$downstream_site[i],]
    
    # Skip if missing data
    if(nrow(upstream) == 0 | nrow(downstream) == 0) next
    
    # Calculate discharge change
    Q_change <- downstream$discharge_m3s - upstream$discharge_m3s
    
    # Calculate expected load based on transport type
    if(Q_change > 0) {
      # GAINING reach - add tributary water at min_Cd
      tributary_flow <- Q_change
      expected_load <- upstream$W_Cd_load + (tributary_flow * min_Cd)
      transport_type <- "Gaining"
    } else {
      # LOSING reach - remove water at upstream concentration  
      expected_load <- upstream$W_Cd_load * (downstream$discharge_m3s / upstream$discharge_m3s)
      transport_type <- "Losing"
    }
    
    # Calculate residual
    actual_load <- downstream$W_Cd_load
    residual <- actual_load - expected_load
    residual_pct <- (residual / expected_load) * 100
    
    # Store results
    load_analysis[[counter]] <- data.frame(
      Date = date,
      reach = actual_segments$reach_name[i],
      upstream_site = actual_segments$upstream_site[i],
      downstream_site = actual_segments$downstream_site[i],
      distance_km = actual_segments$downstream_km[i],
      Q_upstream = upstream$discharge_m3s,
      Q_downstream = downstream$discharge_m3s,
      Q_change = Q_change,
      transport_type = transport_type,
      Cd_upstream = upstream$W_Cd,
      Cd_downstream = downstream$W_Cd,
      upstream_load = upstream$W_Cd_load,
      actual_load = actual_load,
      expected_load = expected_load,
      residual = residual,
      residual_pct = residual_pct
    )
    counter <- counter + 1
  }
}

# Combine results
load_analysis <- do.call(rbind, load_analysis)

# === IDENTIFY HOTSPOTS ===
cat("\n=== HOTSPOT IDENTIFICATION (>20% deviation) ===\n")

# Calculate average residual by reach
residual_summary <- load_analysis %>%
  group_by(reach, upstream_site, downstream_site) %>%
  summarise(
    mean_residual_pct = mean(residual_pct),
    sd_residual_pct = sd(residual_pct),
    n_dates = n(),
    transport_type = unique(transport_type)[1]
  ) %>%
  arrange(desc(abs(mean_residual_pct)))

print(residual_summary)

# Identify consistent hotspots
hotspots <- residual_summary[abs(residual_summary$mean_residual_pct) > 20,]
cat("\n=== REACHES WITH >20% AVERAGE DEVIATION ===\n")
if(nrow(hotspots) > 0) {
  print(hotspots)
} else {
  cat("No reaches exceed 20% average deviation from transport\n")
}

# === SIMPLE VISUALIZATIONS ===
# Plot 1: Deviation from transport by reach
par(mfrow = c(2,1))

# Average residuals by reach
barplot(residual_summary$mean_residual_pct, 
        names.arg = residual_summary$reach,
        main = "Average % Deviation from Physical Transport",
        ylab = "Residual (%)",
        col = ifelse(abs(residual_summary$mean_residual_pct) > 20, "red", "gray"))
abline(h = c(-20, 0, 20), lty = c(2, 1, 2))

# Residuals over time for each reach
plot(load_analysis$distance_km, load_analysis$residual_pct,
     main = "Spatial Pattern of Residuals",
     xlab = "Distance (km)", ylab = "Residual (%)",
     pch = 16, col = as.numeric(factor(load_analysis$Date)))
abline(h = c(-20, 0, 20), lty = c(2, 1, 2))
legend("topright", legend = unique(load_analysis$Date), 
       col = 1:length(unique(load_analysis$Date)), pch = 16, cex = 0.7)

# === KEY FINDINGS ===
cat("\n=== KEY FINDINGS ===\n")
cat("1. Total reaches analyzed:", nrow(actual_segments), "\n")
cat("2. Dates analyzed:", length(unique(load_analysis$Date)), "\n")
cat("3. Reaches with >20% deviation:", nrow(hotspots), "\n")
cat("4. Range of residuals:", round(range(load_analysis$residual_pct), 1), "%\n")

# Identify sources vs sinks
sources <- load_analysis[load_analysis$residual_pct > 20,]
sinks <- load_analysis[load_analysis$residual_pct < -20,]

cat("\n5. Observations with >20% positive residual (SOURCES):", nrow(sources), "\n")
if(nrow(sources) > 0) {
  cat("   Reaches:", unique(sources$reach), "\n")
}

cat("\n6. Observations with >20% negative residual (SINKS):", nrow(sinks), "\n")
if(nrow(sinks) > 0) {
  cat("   Reaches:", unique(sinks$reach), "\n")
}

# === DILUTION EFFECTS FROM MASS BALANCE ===
cat("\n=== DILUTION EFFECTS BY REACH ===\n")

# Calculate dilution effect for each reach
dilution_analysis <- load_analysis %>%
  mutate(
    # Expected concentration if only dilution occurred
    expected_conc = expected_load / Q_downstream,
    # Actual dilution effect (negative = dilution, positive = concentration)
    dilution_effect = (Cd_downstream - Cd_upstream) - (expected_conc - Cd_upstream),
    # Percent change due to dilution vs other processes
    dilution_pct = ((expected_conc - Cd_upstream) / Cd_upstream) * 100,
    process_pct = ((Cd_downstream - expected_conc) / expected_conc) * 100
  )

# Summary of dilution vs process effects
dilution_summary <- dilution_analysis %>%
  group_by(reach) %>%
  summarise(
    mean_dilution_pct = mean(dilution_pct),
    mean_process_pct = mean(process_pct),
    transport_type = unique(transport_type)[1]
  )

cat("Dilution % = change in concentration due to physical mixing\n")
cat("Process % = change in concentration due to biological/geochemical processes\n\n")
print(dilution_summary)

# === STORE KEY OBJECTS FOR FURTHER ANALYSIS ===
transport_results <- list(
  min_Cd = min_Cd,
  actual_segments = actual_segments,
  load_analysis = load_analysis,
  residual_summary = residual_summary,
  hotspots = hotspots,
  dilution_analysis = dilution_analysis,
  dilution_summary = dilution_summary
)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Key objects stored in 'transport_results'\n")

