# Step 1: Basic Transport Framework Analysis for Cadmium
# Simple, direct code for your specific dataset

library(tidyverse)
library(ggplot2)

# ============================================================================
# PART 1: BASIC DATA SETUP
# ============================================================================

# Assume df is already loaded
# Add calculated columns
df$W_Cd_load <- df$W_Cd * df$discharge_m3s
df$Date <- as.Date(df$Date)

# Order sites by distance
df <- df %>% arrange(Distance_km, Date)

# Quick summary
cat("=== DATA SUMMARY ===\n")
cat("Sites:", unique(df$Site), "\n")
cat("Distance range:", range(df$Distance_km), "km\n")
cat("Date range:", range(df$Date), "\n")
cat("Number of observations:", nrow(df), "\n\n")

# ============================================================================
# PART 2: DISCHARGE PATTERNS
# ============================================================================

# Calculate discharge changes between consecutive sites for each date
discharge_changes <- df %>%
  arrange(Date, Distance_km) %>%
  group_by(Date) %>%
  mutate(
    discharge_change = discharge_m3s - lag(discharge_m3s),
    upstream_site = lag(Site),
    reach = paste(upstream_site, "->", Site)
  ) %>%
  filter(!is.na(discharge_change))

# Summary by reach
discharge_summary <- discharge_changes %>%
  group_by(reach) %>%
  summarise(
    mean_Q_change = mean(discharge_change),
    n_gains = sum(discharge_change > 0),
    n_losses = sum(discharge_change < 0)
  )

cat("=== DISCHARGE CHANGES BY REACH ===\n")
print(discharge_summary)
cat("\n")

# Plot discharge along river
plot(df$Distance_km, df$discharge_m3s, 
     main = "Discharge vs Distance",
     xlab = "Distance (km)", 
     ylab = "Discharge (m3/s)",
     col = as.numeric(factor(df$Date)))
legend("topright", legend = unique(df$Date), 
       col = 1:length(unique(df$Date)), pch = 1, cex = 0.5)

# ============================================================================
# PART 3: IDENTIFY MINIMUM CD FOR TRIBUTARY ASSUMPTION
# ============================================================================

min_Cd <- min(df$W_Cd)
min_Cd_5pct <- quantile(df$W_Cd, 0.05)

cat("=== MINIMUM CADMIUM VALUES ===\n")
cat("Absolute minimum W_Cd:", min_Cd, "\n")
cat("5th percentile W_Cd:", min_Cd_5pct, "\n")

# Where does minimum occur?
min_location <- df[which.min(df$W_Cd), c("Site", "Date", "Distance_km", "W_Cd")]
cat("Minimum occurs at:\n")
print(min_location)
cat("\n")

# Cd by site
cd_by_site <- df %>%
  group_by(Site, Distance_km) %>%
  summarise(
    mean_W_Cd = mean(W_Cd),
    min_W_Cd = min(W_Cd),
    max_W_Cd = max(W_Cd)
  )

cat("=== CADMIUM BY SITE ===\n")
print(cd_by_site)
cat("\n")

# ============================================================================
# PART 4: BASIC SPATIAL PATTERNS
# ============================================================================

# Plot Cd vs distance
plot(df$Distance_km, df$W_Cd,
     main = "W_Cd Concentration vs Distance",
     xlab = "Distance (km)",
     ylab = "W_Cd",
     col = "blue")
lines(lowess(df$Distance_km, df$W_Cd), col = "red", lwd = 2)

# Plot discharge vs Cd concentration
plot(df$discharge_m3s, df$W_Cd,
     log = "xy",
     main = "Discharge-Concentration Relationship",
     xlab = "Discharge (m3/s)",
     ylab = "W_Cd")
abline(lm(log(W_Cd) ~ log(discharge_m3s), data = df), col = "red")

# Calculate dilution slope
dilution_model <- lm(log(W_Cd) ~ log(discharge_m3s), data = df)
slope <- coef(dilution_model)[2]
cat("=== DILUTION ANALYSIS ===\n")
cat("Power law slope:", round(slope, 3), "\n")
cat("(slope = -1 means perfect dilution)\n\n")

# ============================================================================
# PART 5: SIMPLE MASS BALANCE
# ============================================================================

# Calculate load changes and expected values
load_analysis <- df %>%
  arrange(Date, Distance_km) %>%
  group_by(Date) %>%
  mutate(
    # Basic load changes
    load_change = W_Cd_load - lag(W_Cd_load),
    Q_change = discharge_m3s - lag(discharge_m3s),
    
    # Expected load if tributary has minimum Cd
    expected_new_load = Q_change * min_Cd,
    expected_total_load = lag(W_Cd_load) + expected_new_load,
    
    # Residual (what can't be explained by transport)
    residual = W_Cd_load - expected_total_load,
    residual_pct = 100 * residual / expected_total_load,
    
    # Reach info
    upstream_site = lag(Site),
    reach = paste(upstream_site, "->", Site)
  ) %>%
  filter(!is.na(load_change))

# Summary of residuals by reach
residual_summary <- load_analysis %>%
  group_by(reach) %>%
  summarise(
    mean_residual_pct = mean(residual_pct, na.rm = TRUE),
    n_sources = sum(residual > 0, na.rm = TRUE),
    n_sinks = sum(residual < 0, na.rm = TRUE)
  ) %>%
  arrange(desc(abs(mean_residual_pct)))

cat("=== MASS BALANCE RESIDUALS ===\n")
cat("Positive = more Cd than expected from transport\n")
cat("Negative = less Cd than expected from transport\n\n")
print(residual_summary)
cat("\n")

# Plot residuals
boxplot(residual_pct ~ reach, data = load_analysis,
        main = "Deviations from Simple Transport",
        ylab = "Residual (%)",
        las = 2)
abline(h = 0, col = "red", lty = 2)

# ============================================================================
# PART 6: ENVIRONMENTAL CONTEXT
# ============================================================================

# Average environmental conditions by site
env_by_site <- df %>%
  group_by(Site, Distance_km) %>%
  summarise(
    mean_pH = mean(pH, na.rm = TRUE),
    mean_DOC = mean(DOC, na.rm = TRUE),
    mean_temp = mean(daily_temp, na.rm = TRUE),
    mean_alkalinity = mean(Alkalinity, na.rm = TRUE)
  )

cat("=== ENVIRONMENTAL CONDITIONS BY SITE ===\n")
print(env_by_site)
cat("\n")

# Plot environmental gradients
par(mfrow = c(2, 2))
plot(env_by_site$Distance_km, env_by_site$mean_pH, 
     type = "b", main = "pH", xlab = "Distance (km)", ylab = "pH")
plot(env_by_site$Distance_km, env_by_site$mean_DOC,
     type = "b", main = "DOC", xlab = "Distance (km)", ylab = "DOC")
plot(env_by_site$Distance_km, env_by_site$mean_temp,
     type = "b", main = "Temperature", xlab = "Distance (km)", ylab = "Temp (C)")
plot(env_by_site$Distance_km, env_by_site$mean_alkalinity,
     type = "b", main = "Alkalinity", xlab = "Distance (km)", ylab = "Alkalinity")
par(mfrow = c(1, 1))

# ============================================================================
# SUMMARY FOR STEP 2
# ============================================================================

cat("=== KEY FINDINGS FOR STEP 2 ===\n")
cat("Tributary Cd assumption:", min_Cd, "\n")
cat("Dilution slope:", round(slope, 3), "\n")

# Identify hotspots (>20% deviation from transport)
hotspots <- residual_summary %>%
  filter(abs(mean_residual_pct) > 20) %>%
  pull(reach)

cat("Hotspot reaches (>20% deviation):", 
    if(length(hotspots) > 0) paste(hotspots, collapse = ", ") else "None", "\n")

# Save key results
transport_results <- list(
  min_Cd = min_Cd,
  min_Cd_5pct = min_Cd_5pct,
  dilution_slope = slope,
  discharge_summary = discharge_summary,
  cd_by_site = cd_by_site,
  residual_summary = residual_summary,
  load_analysis = load_analysis
)

cat("\nResults saved in 'transport_results' object\n")

