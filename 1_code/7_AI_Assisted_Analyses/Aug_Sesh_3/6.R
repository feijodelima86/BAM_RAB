# Three-panel plot showing process effects (dilution removed) by fraction
# This shows the biological/geochemical signal after accounting for dilution

# Add distance information to fraction_dilution
fraction_dilution_spatial <- merge(
  fraction_dilution,
  transport_results$actual_segments[,c("reach_name", "downstream_km")],
  by.x = "reach", by.y = "reach_name"
)

# Get ordered unique dates
unique_dates <- sort(unique(fraction_dilution_spatial$Date))
n_dates <- length(unique_dates)

# Create blue gradient from light to dark
blue_colors <- colorRampPalette(c("#E6F3FF", "#0033CC"))(n_dates)

# Create named vector for easy lookup
date_colors <- setNames(blue_colors, unique_dates)


# Set up three panels
par(mfrow = c(3,1))

# Panel 3: Suspended fraction process effects
suspended_process <- fraction_dilution_spatial[fraction_dilution_spatial$fraction == "suspended",]
plot(suspended_process$downstream_km, suspended_process$process_pct,
     main = "Suspended: Process Effects",
     xlab = "Distance (km)", ylab = "Process Effect (%)",
     pch = 16, col  = date_colors[as.character(dissolved_data$Date)],
     ylim = range(fraction_dilution_spatial$process_pct, na.rm = TRUE))
abline(h = c(-20, 0, 20), lty = c(2, 1, 2), col = "gray")
abline(h = 0, lty = 1, lwd = 2)

# Panel 2: Colloidal fraction process effects
colloidal_process <- fraction_dilution_spatial[fraction_dilution_spatial$fraction == "colloidal",]
plot(colloidal_process$downstream_km, colloidal_process$process_pct,
     main = "Colloidal: Process Effects",
     xlab = "Distance (km)", ylab = "Process Effect (%)",
     pch = 16, col  = date_colors[as.character(dissolved_data$Date)],
     ylim = range(fraction_dilution_spatial$process_pct, na.rm = TRUE))
abline(h = c(-20, 0, 20), lty = c(2, 1, 2), col = "gray")
abline(h = 0, lty = 1, lwd = 2)

# Panel 1: Dissolved fraction process effects
dissolved_process <- fraction_dilution_spatial[fraction_dilution_spatial$fraction == "dissolved",]
plot(dissolved_process$downstream_km, dissolved_process$process_pct,
     main = "Dissolved: Process Effects",
     xlab = "Distance (km)", ylab = "Process Effect (%)",
     pch = 16, col  = date_colors[as.character(dissolved_data$Date)],
     ylim = range(fraction_dilution_spatial$process_pct, na.rm = TRUE))
abline(h = c(-20, 0, 20), lty = c(2, 1, 2), col = "gray")
abline(h = 0, lty = 1, lwd = 2)
legend("topright", legend = unique(dissolved_process$Date), 
       col = date_colors[as.character(dissolved_data$Date)],, pch = 16, cex = 0.5)
text(min(dissolved_process$downstream_km), max(fraction_dilution_spatial$process_pct)*0.9,
     "Positive = Source", adj = 0, cex = 0.7)
text(min(dissolved_process$downstream_km), min(fraction_dilution_spatial$process_pct)*0.9,
     "Negative = Sink", adj = 0, cex = 0.7)

# Reset to single panel
par(mfrow = c(1,1))

# Compare process effects between fractions at each location
cat("\n=== PROCESS EFFECTS BY FRACTION AND LOCATION ===\n")
process_comparison <- fraction_dilution_spatial %>%
  group_by(fraction, downstream_km) %>%
  summarise(
    mean_process = mean(process_pct, na.rm = TRUE),
    consistency = ifelse(all(process_pct > 0, na.rm = TRUE), "Always Source",
                         ifelse(all(process_pct < 0, na.rm = TRUE), "Always Sink", 
                                "Variable"))
  ) %>%
  pivot_wider(names_from = fraction, values_from = c(mean_process, consistency))

print(process_comparison)

# Identify which fraction has strongest process effects at each location
cat("\n=== DOMINANT PROCESSES BY LOCATION ===\n")
for(km in unique(fraction_dilution_spatial$downstream_km)) {
  km_data <- fraction_dilution_spatial[fraction_dilution_spatial$downstream_km == km,]
  mean_by_frac <- tapply(km_data$process_pct, km_data$fraction, mean, na.rm = TRUE)
  dominant <- names(mean_by_frac)[which.max(abs(mean_by_frac))]
  cat("At", km, "km: Strongest process in", dominant, "fraction (", 
      round(mean_by_frac[dominant], 1), "%)\n")
}