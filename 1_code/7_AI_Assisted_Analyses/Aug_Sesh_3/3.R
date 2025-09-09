# Multi-panel scatterplot of deviations from simple transport by date

# Get unique dates and set up multi-panel layout
unique_dates <- sort(unique(load_analysis$Date))
n_dates <- length(unique_dates)

# Set up panel layout (adjust rows/cols as needed)
n_cols <- 3
n_rows <- ceiling(n_dates / n_cols)
par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))

# Create a plot for each date
for(date in unique_dates) {
  # Get data for this date
  date_data <- load_analysis[load_analysis$Date == date,]
  
  # Match to actual segments to get distances
  date_data <- merge(date_data, actual_segments, by = "reach")
  date_data$mid_distance <- (date_data$upstream_km + date_data$downstream_km) / 2
  
  # Create scatterplot
  plot(date_data$mid_distance, date_data$residual_pct,
       xlim = range(sites_ordered$Distance_km),
       ylim = range(load_analysis$residual_pct, na.rm = TRUE),
       xlab = "Distance (km)",
       ylab = "Residual (%)",
       main = as.character(date),
       pch = 19,
       col = ifelse(date_data$residual_pct > 0, "red", "blue"),
       cex = 1.5)
  
  # Add reference lines
  abline(h = 0, lty = 2, col = "gray40")
  abline(h = c(-20, 20), lty = 3, col = "gray70")
  
  # Add segment labels if not too crowded
  if(nrow(date_data) <= 6) {
    text(date_data$mid_distance, date_data$residual_pct,
         labels = date_data$reach,
         pos = 3, cex = 0.6, offset = 0.3)
  }
  
  # Add grid
  grid(col = "gray90")
}

# Reset layout
par(mfrow = c(1, 1))

# Alternative version: Single plot with dates as different colors
plot(NULL, 
     xlim = range(sites_ordered$Distance_km),
     ylim = range(load_analysis$residual_pct, na.rm = TRUE),
     xlab = "Distance (km)",
     ylab = "Residual (% deviation from expected)",
     main = "Deviations from Simple Transport - All Dates")

# Color palette for dates
colors <- rainbow(n_dates)

for(i in 1:n_dates) {
  date <- unique_dates[i]
  date_data <- load_analysis[load_analysis$Date == date,]
  date_data <- merge(date_data, actual_segments, by = "reach")
  date_data$mid_distance <- (date_data$upstream_km + date_data$downstream_km) / 2
  
  points(date_data$mid_distance, date_data$residual_pct,
         col = colors[i], pch = 19, cex = 1.5)
  lines(date_data$mid_distance, date_data$residual_pct,
        col = colors[i], lty = 2)
}

abline(h = 0, lty = 1, col = "black", lwd = 2)
abline(h = c(-20, 20), lty = 3, col = "gray50")
grid(col = "gray90")

legend("topright", 
       legend = as.character(unique_dates),
       col = colors, 
       pch = 19, 
       cex = 0.7,
       ncol = 2)