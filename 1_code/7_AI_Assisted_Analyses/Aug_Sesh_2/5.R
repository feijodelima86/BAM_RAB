# Simple Base R Plots for Biological-Chemical Relationships
# Using data from the full spectrum analysis

# =================== BASIC RELATIONSHIP PLOTS ===================

# Set up plotting parameters
par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))

# Color scheme for compartments
comp_colors <- c("EPIL" = "brown", "EPIP" = "green", "FILA" = "blue")

# Plot 1: Biological change vs Dissolved Cd change by compartment
plot(full_spectrum_data$Weighted_biomass_change, full_spectrum_data$Dissolved_change_interp,
     col = comp_colors[full_spectrum_data$Compartment], pch = 16, cex = 0.8,
     xlab = "Weighted Biomass Change", ylab = "Dissolved Cd Change Rate",
     main = "Biomass vs Dissolved Cd")
legend("topright", legend = names(comp_colors), col = comp_colors, pch = 16, cex = 0.8)
abline(lm(Dissolved_change_interp ~ Weighted_biomass_change, data = full_spectrum_data), lty = 2)

# Plot 2: Biological change vs Suspended Cd change by compartment  
plot(full_spectrum_data$Weighted_biomass_change, full_spectrum_data$Suspended_change_interp,
     col = comp_colors[full_spectrum_data$Compartment], pch = 16, cex = 0.8,
     xlab = "Weighted Biomass Change", ylab = "Suspended Cd Change Rate", 
     main = "Biomass vs Suspended Cd")
legend("topright", legend = names(comp_colors), col = comp_colors, pch = 16, cex = 0.8)
abline(lm(Suspended_change_interp ~ Weighted_biomass_change, data = full_spectrum_data), lty = 2)

# Plot 3: Biological change vs Colloidal Cd change by compartment
plot(full_spectrum_data$Weighted_biomass_change, full_spectrum_data$Colloidal_change_interp,
     col = comp_colors[full_spectrum_data$Compartment], pch = 16, cex = 0.8,
     xlab = "Weighted Biomass Change", ylab = "Colloidal Cd Change Rate",
     main = "Biomass vs Colloidal Cd") 
legend("topright", legend = names(comp_colors), col = comp_colors, pch = 16, cex = 0.8)
abline(lm(Colloidal_change_interp ~ Weighted_biomass_change, data = full_spectrum_data), lty = 2)

# Plot 4: Cd change vs Dissolved Cd change by compartment
plot(full_spectrum_data$Weighted_Cd_change, full_spectrum_data$Dissolved_change_interp,
     col = comp_colors[full_spectrum_data$Compartment], pch = 16, cex = 0.8,
     xlab = "Weighted Cd Standing Crop Change", ylab = "Dissolved Cd Change Rate",
     main = "Cd Crop vs Dissolved Cd")
legend("topright", legend = names(comp_colors), col = comp_colors, pch = 16, cex = 0.8)
abline(lm(Dissolved_change_interp ~ Weighted_Cd_change, data = full_spectrum_data), lty = 2)

# Plot 5: Cd change vs Suspended Cd change by compartment
plot(full_spectrum_data$Weighted_Cd_change, full_spectrum_data$Suspended_change_interp,
     col = comp_colors[full_spectrum_data$Compartment], pch = 16, cex = 0.8,
     xlab = "Weighted Cd Standing Crop Change", ylab = "Suspended Cd Change Rate",
     main = "Cd Crop vs Suspended Cd")
legend("topright", legend = names(comp_colors), col = comp_colors, pch = 16, cex = 0.8)
abline(lm(Suspended_change_interp ~ Weighted_Cd_change, data = full_spectrum_data), lty = 2)

# Plot 6: Three-way comparison
plot(full_spectrum_data$Dissolved_change_interp, full_spectrum_data$Suspended_change_interp,
     col = comp_colors[full_spectrum_data$Compartment], pch = 16, cex = 0.8,
     xlab = "Dissolved Cd Change Rate", ylab = "Suspended Cd Change Rate",
     main = "Dissolved vs Suspended Cd")
legend("topright", legend = names(comp_colors), col = comp_colors, pch = 16, cex = 0.8)
abline(lm(Suspended_change_interp ~ Dissolved_change_interp, data = full_spectrum_data), lty = 2)

# =================== TEMPORAL PATTERNS ===================

par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# Convert Date to numeric for plotting
full_spectrum_data$Date_numeric <- as.numeric(full_spectrum_data$Date)

# Plot 7: Biomass change over time by compartment
plot(full_spectrum_data$Date_numeric, full_spectrum_data$Weighted_biomass_change,
     col = comp_colors[full_spectrum_data$Compartment], pch = 16, cex = 0.8,
     xlab = "Date", ylab = "Weighted Biomass Change",
     main = "Biomass Change Over Time", xaxt = "n")
# Add date labels
date_breaks <- seq(min(full_spectrum_data$Date_numeric, na.rm = TRUE), 
                   max(full_spectrum_data$Date_numeric, na.rm = TRUE), 
                   length.out = 5)
axis(1, at = date_breaks, labels = as.Date(date_breaks, origin = "1970-01-01"), cex.axis = 0.8)
legend("topright", legend = names(comp_colors), col = comp_colors, pch = 16, cex = 0.8)

# Plot 8: Dissolved Cd change over time
plot(full_spectrum_data$Date_numeric, full_spectrum_data$Dissolved_change_interp,
     pch = 16, cex = 0.8, col = "red",
     xlab = "Date", ylab = "Dissolved Cd Change Rate", 
     main = "Dissolved Cd Change Over Time", xaxt = "n")
axis(1, at = date_breaks, labels = as.Date(date_breaks, origin = "1970-01-01"), cex.axis = 0.8)
abline(h = 0, lty = 2, col = "gray")

# Plot 9: Suspended Cd change over time
plot(full_spectrum_data$Date_numeric, full_spectrum_data$Suspended_change_interp,
     pch = 16, cex = 0.8, col = "blue",
     xlab = "Date", ylab = "Suspended Cd Change Rate",
     main = "Suspended Cd Change Over Time", xaxt = "n")
axis(1, at = date_breaks, labels = as.Date(date_breaks, origin = "1970-01-01"), cex.axis = 0.8)
abline(h = 0, lty = 2, col = "gray")

# Plot 10: All three Cd fractions over time
plot(full_spectrum_data$Date_numeric, full_spectrum_data$Dissolved_change_interp,
     pch = 16, cex = 0.6, col = "red", ylim = range(c(full_spectrum_data$Dissolved_change_interp,
                                                      full_spectrum_data$Colloidal_change_interp,
                                                      full_spectrum_data$Suspended_change_interp), na.rm = TRUE),
     xlab = "Date", ylab = "Cd Change Rate",
     main = "All Cd Fractions Over Time", xaxt = "n")
points(full_spectrum_data$Date_numeric, full_spectrum_data$Colloidal_change_interp, 
       pch = 17, cex = 0.6, col = "orange")
points(full_spectrum_data$Date_numeric, full_spectrum_data$Suspended_change_interp, 
       pch = 18, cex = 0.6, col = "blue")
axis(1, at = date_breaks, labels = as.Date(date_breaks, origin = "1970-01-01"), cex.axis = 0.8)
abline(h = 0, lty = 2, col = "gray")
legend("topright", legend = c("Dissolved", "Colloidal", "Suspended"), 
       col = c("red", "orange", "blue"), pch = c(16, 17, 18), cex = 0.8)

# =================== SPATIAL PATTERNS ===================

par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

# Plot 11: Biomass change vs distance by compartment
plot(full_spectrum_data$Distance_km, full_spectrum_data$Weighted_biomass_change,
     col = comp_colors[full_spectrum_data$Compartment], pch = 16, cex = 0.8,
     xlab = "Distance (km)", ylab = "Weighted Biomass Change",
     main = "Biomass Change vs Distance")
legend("topright", legend = names(comp_colors), col = comp_colors, pch = 16, cex = 0.8)

# Plot 12: Dissolved Cd change vs distance
plot(full_spectrum_data$Distance_km, full_spectrum_data$Dissolved_change_interp,
     pch = 16, cex = 0.8, col = "red",
     xlab = "Distance (km)", ylab = "Dissolved Cd Change Rate",
     main = "Dissolved Cd Change vs Distance")
abline(h = 0, lty = 2, col = "gray")

# Plot 13: Suspended Cd change vs distance  
plot(full_spectrum_data$Distance_km, full_spectrum_data$Suspended_change_interp,
     pch = 16, cex = 0.8, col = "blue",
     xlab = "Distance (km)", ylab = "Suspended Cd Change Rate",
     main = "Suspended Cd Change vs Distance")
abline(h = 0, lty = 2, col = "gray")

# Plot 14: Site-level averages
site_averages <- aggregate(cbind(Weighted_biomass_change, Dissolved_change_interp, Suspended_change_interp) ~ 
                             Site + Distance_km, data = full_spectrum_data, FUN = mean, na.rm = TRUE)

plot(site_averages$Distance_km, site_averages$Weighted_biomass_change,
     pch = 16, cex = 1.2, col = "darkgreen",
     xlab = "Distance (km)", ylab = "Mean Change Rates",
     main = "Site Averages vs Distance",
     ylim = range(c(site_averages$Weighted_biomass_change, 
                    site_averages$Dissolved_change_interp * 1000,
                    site_averages$Suspended_change_interp * 1000), na.rm = TRUE))
points(site_averages$Distance_km, site_averages$Dissolved_change_interp * 1000, 
       pch = 17, cex = 1.2, col = "red")
points(site_averages$Distance_km, site_averages$Suspended_change_interp * 1000, 
       pch = 18, cex = 1.2, col = "blue")
text(site_averages$Distance_km, site_averages$Weighted_biomass_change, 
     labels = site_averages$Site, pos = 3, cex = 0.8)
legend("topright", legend = c("Biomass Change", "Dissolved Cd × 1000", "Suspended Cd × 1000"), 
       col = c("darkgreen", "red", "blue"), pch = c(16, 17, 18), cex = 0.8)
abline(h = 0, lty = 2, col = "gray")

# =================== COMPARTMENT-SPECIFIC CORRELATION PLOTS ===================

par(mfrow = c(1, 3), mar = c(4, 4, 3, 2))

# Individual compartment plots with trend lines
compartments <- c("EPIL", "EPIP", "FILA")

for(i in 1:3) {
  comp_data <- subset(full_spectrum_data, Compartment == compartments[i])
  
  plot(comp_data$Weighted_biomass_change, comp_data$Suspended_change_interp,
       pch = 16, cex = 1, col = comp_colors[compartments[i]],
       xlab = "Weighted Biomass Change", ylab = "Suspended Cd Change Rate",
       main = paste(compartments[i], "- Biomass vs Suspended Cd"))
  
  # Add trend line if enough data
  if(nrow(comp_data) > 5 && var(comp_data$Weighted_biomass_change, na.rm = TRUE) > 0) {
    lm_fit <- lm(Suspended_change_interp ~ Weighted_biomass_change, data = comp_data)
    abline(lm_fit, col = "black", lwd = 2)
    
    # Add correlation coefficient
    cor_val <- cor(comp_data$Weighted_biomass_change, comp_data$Suspended_change_interp, use = "complete.obs")
    text(min(comp_data$Weighted_biomass_change, na.rm = TRUE), 
         max(comp_data$Suspended_change_interp, na.rm = TRUE),
         labels = paste("r =", round(cor_val, 3)), pos = 4, cex = 1.2)
  }
  
  abline(h = 0, v = 0, lty = 2, col = "gray")
}

# Reset plotting parameters
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))