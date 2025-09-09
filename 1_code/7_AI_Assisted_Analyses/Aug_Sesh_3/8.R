# Create progressive blue color scheme for dates
# Earlier dates = lighter blue, Later dates = darker blue

# Get ordered unique dates
unique_dates <- sort(unique(fraction_dilution_spatial$Date))
n_dates <- length(unique_dates)

# Create blue gradient from light to dark
blue_colors <- colorRampPalette(c("#E6F3FF", "#0033CC"))(n_dates)

# Create named vector for easy lookup
date_colors <- setNames(blue_colors, unique_dates)

# === EXAMPLE: Three-panel plot with blue gradient ===
par(mfrow = c(1,3))

# Panel 1: Dissolved fraction
dissolved_data <- fraction_dilution_spatial[fraction_dilution_spatial$fraction == "dissolved",]
plot(dissolved_data$downstream_km, dissolved_data$residual_pct,
     main = "Dissolved Cd Residuals",
     xlab = "Distance (km)", ylab = "Residual (%)",
     pch = 16, col = date_colors[as.character(dissolved_data$Date)],
     ylim = range(fraction_dilution_spatial$residual_pct, na.rm = TRUE))
abline(h = c(-20, 0, 20), lty = c(2, 1, 2))
legend("topright", legend = unique_dates, 
       col = blue_colors, pch = 16, cex = 0.5,
       title = "Date")

# Panel 2: Colloidal fraction
colloidal_data <- fraction_residuals_spatial[fraction_residuals_spatial$fraction == "colloidal",]
plot(colloidal_data$downstream_km, colloidal_data$residual_pct,
     main = "Colloidal Cd Residuals",
     xlab = "Distance (km)", ylab = "Residual (%)",
     pch = 16, col = date_colors[as.character(colloidal_data$Date)],
     ylim = range(fraction_residuals_spatial$residual_pct, na.rm = TRUE))
abline(h = c(-20, 0, 20), lty = c(2, 1, 2))

# Panel 3: Suspended fraction
suspended_data <- fraction_residuals_spatial[fraction_residuals_spatial$fraction == "suspended",]
plot(suspended_data$downstream_km, suspended_data$residual_pct,
     main = "Suspended Cd Residuals",
     xlab = "Distance (km)", ylab = "Residual (%)",
     pch = 16, col = date_colors[as.character(suspended_data$Date)],
     ylim = range(fraction_residuals_spatial$residual_pct, na.rm = TRUE))
abline(h = c(-20, 0, 20), lty = c(2, 1, 2))

par(mfrow = c(1,1))

# === ALTERNATIVE: Different blue palettes ===
# You can also try these other blue gradients:

# Option 2: Sky to Navy
# blue_colors <- colorRampPalette(c("#87CEEB", "#000080"))(n_dates)

# Option 3: Pastel to Royal Blue  
# blue_colors <- colorRampPalette(c("#B8D4E3", "#4169E1"))(n_dates)

# Option 4: Using built-in blues
# library(RColorBrewer)
# blue_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(n_dates)

# === FOR ANY EXISTING PLOT ===
# To apply to any plot, just replace the col argument:
# Instead of: col = as.numeric(factor(data$Date))
# Use: col = date_colors[as.character(data$Date)]