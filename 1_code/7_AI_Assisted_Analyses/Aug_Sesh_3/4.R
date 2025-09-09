# Step 1: Simple Transport Framework Analysis
# Only actual river segments and requested plots

library(tidyverse)

# ============================================================================
# BASIC SETUP
# ============================================================================

# Add calculated columns
df$W_Cd_load <- df$W_Cd * df$discharge_m3s

# Find minimum Cd for tributary assumption
min_Cd <- min(df$W_Cd)

cat("=== BASIC INFO ===\n")
cat("Minimum W_Cd for tributary assumption:", min_Cd, "\n")
cat("Sites in dataset:", unique(df$Site), "\n")
cat("Distance range:", range(df$Distance_km), "km\n\n")

# ============================================================================
# IDENTIFY ACTUAL RIVER SEGMENTS
# ============================================================================

# Get unique sites ordered by distance
sites_ordered <- df %>%
  select(Site, Distance_km) %>%
  distinct() %>%
  arrange(Distance_km)

cat("Sites in order (upstream to downstream):\n")
print(sites_ordered)
cat("\n")

# Create list of actual adjacent segments
actual_segments <- data.frame(
  upstream_site = sites_ordered$Site[-nrow(sites_ordered)],
  downstream_site = sites_ordered$Site[-1],
  upstream_km = sites_ordered$Distance_km[-nrow(sites_ordered)],
  downstream_km = sites_ordered$Distance_km[-1]
)
actual_segments$reach <- paste(actual_segments$upstream_site, "->", 
                               actual_segments$downstream_site)

cat("Actual river segments:\n")
print(actual_segments$reach)
cat("\n")

# ============================================================================
# MASS BALANCE FOR ACTUAL SEGMENTS ONLY
# ============================================================================

# Function to calculate mass balance for a specific date and segment
calculate_segment_balance <- function(date, upstream_site, downstream_site) {
  upstream <- df[df$Date == date & df$Site == upstream_site,]
  downstream <- df[df$Date == date & df$Site == downstream_site,]
  
  if(nrow(upstream) == 0 || nrow(downstream) == 0) return(NULL)
  
  # Calculate loads
  upstream_load <- upstream$W_Cd * upstream$discharge_m3s
  downstream_load <- downstream$W_Cd * downstream$discharge_m3s
  
  # Discharge change (tributary input or abstraction)
  Q_change <- downstream$discharge_m3s - upstream$discharge_m3s
  
  # Expected load depends on whether gaining or losing water
  if(Q_change > 0) {
    # Gaining reach: tributary adds water with minimum Cd
    expected_new_load <- Q_change * min_Cd
    expected_total_load <- upstream_load + expected_new_load
  } else {
    # Losing reach: abstraction removes water at upstream concentration
    # Expected downstream load is proportional to discharge ratio
    expected_total_load <- upstream_load * (downstream$discharge_m3s / upstream$discharge_m3s)
  }
  
  # Residual
  residual <- downstream_load - expected_total_load
  residual_pct <- 100 * residual / expected_total_load
  
  return(data.frame(
    Date = date,
    reach = paste(upstream_site, "->", downstream_site),
    residual = residual,
    residual_pct = residual_pct
  ))
}

# Calculate for all dates and actual segments
all_results <- list()
for(i in 1:nrow(actual_segments)) {
  for(date in unique(df$Date)) {
    result <- calculate_segment_balance(date, 
                                        actual_segments$upstream_site[i],
                                        actual_segments$downstream_site[i])
    if(!is.null(result)) {
      all_results[[length(all_results) + 1]] <- result
    }
  }
}

load_analysis <- do.call(rbind, all_results)

# Summary by reach
residual_summary <- load_analysis %>%
  group_by(reach) %>%
  summarise(
    n_obs = n(),
    mean_residual_pct = mean(residual_pct, na.rm = TRUE),
    sd_residual_pct = sd(residual_pct, na.rm = TRUE)
  ) %>%
  arrange(desc(abs(mean_residual_pct)))

cat("=== MASS BALANCE DEVIATIONS FROM SIMPLE TRANSPORT ===\n")
cat("Positive = more Cd than expected from dilution\n")
cat("Negative = less Cd than expected\n\n")
print(residual_summary)
cat("\n")

# ============================================================================
# PLOT 1: DEVIATIONS FROM SIMPLE TRANSPORT
# ============================================================================

# Order reaches by distance for plotting
load_analysis$reach <- factor(load_analysis$reach, 
                              levels = actual_segments$reach)

boxplot(residual_pct ~ reach, data = load_analysis,
        main = "Deviations from Simple Transport Model",
        ylab = "Residual (% difference from expected)",
        xlab = "River Segment",
        las = 2,
        col = ifelse(residual_summary$mean_residual_pct > 0, "coral", "lightblue"))
abline(h = 0, col = "red", lty = 2, lwd = 2)
abline(h = c(-20, 20), col = "gray", lty = 3)

# ============================================================================
# PLOT 2: THREE PANEL DISCHARGE-CONCENTRATION BY PARTICLE SIZE
# ============================================================================

# Calculate concentrations for each fraction
df$Dissolved_conc <- df$W_Cd * df$Dissolved_pct / 100
df$Colloidal_conc <- df$W_Cd * df$Colloidal_pct / 100
df$Suspended_conc <- df$W_Cd * df$Suspended_pct / 100

par(mfrow = c(1, 3))

# Dissolved
valid_dissolved <- df$Dissolved_conc > 0 & !is.na(df$Dissolved_conc)
plot(df$discharge_m3s[valid_dissolved], 
     df$Dissolved_conc[valid_dissolved],
     log = "xy",
     main = "Dissolved",
     xlab = "Discharge (m3/s)",
     ylab = "Dissolved Cd",
     col = "blue", pch = 16)
if(sum(valid_dissolved) > 2) {
  dissolved_model <- lm(log(Dissolved_conc) ~ log(discharge_m3s), 
                        data = df[valid_dissolved,])
  abline(dissolved_model, col = "blue", lwd = 2)
  legend("topright", legend = paste("slope =", round(coef(dissolved_model)[2], 2)),
         bty = "n", text.col = "blue")
}

# Colloidal  
valid_colloidal <- df$Colloidal_conc > 0 & !is.na(df$Colloidal_conc)
plot(df$discharge_m3s[valid_colloidal], 
     df$Colloidal_conc[valid_colloidal],
     log = "xy",
     main = "Colloidal",
     xlab = "Discharge (m3/s)",
     ylab = "Colloidal Cd",
     col = "green4", pch = 16)
if(sum(valid_colloidal) > 2) {
  colloidal_model <- lm(log(Colloidal_conc) ~ log(discharge_m3s), 
                        data = df[valid_colloidal,])
  abline(colloidal_model, col = "green4", lwd = 2)
  legend("topright", legend = paste("slope =", round(coef(colloidal_model)[2], 2)),
         bty = "n", text.col = "green4")
}

# Suspended
valid_suspended <- df$Suspended_conc > 0 & !is.na(df$Suspended_conc)
plot(df$discharge_m3s[valid_suspended], 
     df$Suspended_conc[valid_suspended],
     log = "xy",
     main = "Suspended",
     xlab = "Discharge (m3/s)",
     ylab = "Suspended Cd",
     col = "brown", pch = 16)
if(sum(valid_suspended) > 2) {
  suspended_model <- lm(log(Suspended_conc) ~ log(discharge_m3s), 
                        data = df[valid_suspended,])
  abline(suspended_model, col = "brown", lwd = 2)
  legend("topright", legend = paste("slope =", round(coef(suspended_model)[2], 2)),
         bty = "n", text.col = "brown")
}

par(mfrow = c(1, 1))

# ============================================================================
# SUMMARY
# ============================================================================

cat("=== KEY FINDINGS ===\n")
cat("Segments with >20% deviation from transport:\n")
hotspots <- residual_summary %>% filter(abs(mean_residual_pct) > 20)
if(nrow(hotspots) > 0) {
  print(hotspots)
} else {
  cat("None\n")
}

# Save results
transport_results <- list(
  actual_segments = actual_segments,
  residual_summary = residual_summary,
  load_analysis = load_analysis
)