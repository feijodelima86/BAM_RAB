# ===============================================================
# COMPREHENSIVE CADMIUM SPECIATION ANALYSIS
# ===============================================================


# Load required libraries
library(dplyr)
library(tidyr)
library(betareg)
library(car)

# Read subset data

df <- read.csv("2_incremental/AI_Subsets/cd.csv")

# Calculate cadmium fractions (using average of wavelengths 111 and 114)
df$UF_Cd <- (df$UF_Cd114)  # Truly dissolved
df$F_Cd <- (df$F_Cd114)      # Dissolved + colloidal
df$W_Cd <- (df$W_Cd114)      # Total (whole water)

# Calculate fraction percentages with negative handling
df$Dissolved_pct <- ifelse((df$UF_Cd / df$W_Cd) < 0, 0, (df$UF_Cd / df$W_Cd) * 100)
df$Colloidal_pct <- ifelse(((df$F_Cd - df$UF_Cd) / df$W_Cd) < 0, 0, 
                           ((df$F_Cd - df$UF_Cd) / df$W_Cd) * 100)
df$Suspended_pct <- ifelse(((df$W_Cd - df$F_Cd) / df$W_Cd) < 0, 0, 
                           ((df$W_Cd - df$F_Cd) / df$W_Cd) * 100)

# Add flags for fraction presence/absence
df$Dissolved_present <- df$UF_Cd >= 0 & (df$UF_Cd / df$W_Cd) >= 0
df$Colloidal_present <- (df$F_Cd - df$UF_Cd) >= 0
df$Suspended_present <- (df$W_Cd - df$F_Cd) >= 0

# Get unique water samples
speciation_unique <- df %>%
  filter(!is.na(Dissolved_pct) & !is.na(Colloidal_pct) & !is.na(Suspended_pct)) %>%
  distinct(SAMPLING_DATE, SITE, .keep_all = TRUE) %>%
  select(Date = SAMPLING_DATE, Site = SITE, 
         Dissolved_pct, Colloidal_pct, Suspended_pct,
         Dissolved_present, Colloidal_present, Suspended_present,
         UF_Cd, F_Cd, W_Cd) %>%
  arrange(Date, Site)


# Add correct spatial context
site_distances <- data.frame(
  Site = c("WS", "DL", "GR", "GC", "BG", "BN"),
  Distance_km = c(0, 44, 64.8, 89.2, 144.3, 167.8),
  Position = c("Upstream", "Upper-Middle", "Middle", "Middle-Lower", "Lower", "Downstream")
)

# Merge with speciation data
speciation_spatial <- merge(speciation_unique, site_distances, by = "Site") %>%
  arrange(Distance_km, Date)


# ===============================================================
# RECONSTRUCTED TOTAL CADMIUM CONCENTRATION ANALYSIS - GAMMA GLM (CURVES)
# ===============================================================

# Function to prepare data for Gamma GLM (reconstructed total from fractions)
prepare_gamma_data_total <- function(data) {
  # Calculate individual fractions
  # Dissolved = UF_Cd
  # Colloidal = F_Cd - UF_Cd
  # Suspended = W_Cd - F_Cd
  
  clean_data <- data[!is.na(data$UF_Cd) & !is.na(data$F_Cd) & !is.na(data$W_Cd), ]
  
  # Calculate each fraction, setting negative values to 0
  clean_data$dissolved <- ifelse(clean_data$UF_Cd < 0, 0, clean_data$UF_Cd)
  clean_data$colloidal <- ifelse((clean_data$F_Cd - clean_data$UF_Cd) < 0, 0, 
                                 clean_data$F_Cd - clean_data$UF_Cd)
  clean_data$suspended <- ifelse((clean_data$W_Cd - clean_data$F_Cd) < 0, 0, 
                                 clean_data$W_Cd - clean_data$F_Cd)
  
  # Reconstruct total as sum of positive fractions
  clean_data$fraction_conc <- clean_data$dissolved + clean_data$colloidal + clean_data$suspended
  
  # Filter for positive total values
  clean_data <- clean_data[clean_data$fraction_conc > 0, ]
  clean_data$fraction_name <- "Reconstructed Total"
  
  # Convert date variables
  clean_data$Date_parsed <- as.Date(clean_data$Date, format = "%m/%d/%Y")
  clean_data$Day_of_year <- as.numeric(format(clean_data$Date_parsed, "%j"))
  
  # Create colors
  min_day <- min(clean_data$Day_of_year)
  max_day <- max(clean_data$Day_of_year)
  color_func <- colorRampPalette(c("lightblue", "navy"))
  clean_data$point_color <- color_func(100)[round(((clean_data$Day_of_year - min_day) / (max_day - min_day)) * 99) + 1]
  
  return(clean_data)
}

# Function to fit Gamma GLM with log link (ensures positive predictions)
fit_gamma_curve_model <- function(data) {
  glm(fraction_conc ~ Distance_km + Day_of_year + Distance_km * Day_of_year, 
      data = data, family = Gamma(link = "log"))
}

# Function to create smooth curve predictions
create_curve_predictions <- function(model, data) {
  distance_seq <- seq(0, 170, by = 5)  # Finer resolution for smooth curves
  early_day <- min(data$Day_of_year)
  mid_day <- median(data$Day_of_year)
  late_day <- max(data$Day_of_year)
  
  # Predict on response scale (automatically back-transformed from log)
  pred_early <- predict(model, newdata = data.frame(Distance_km = distance_seq, Day_of_year = early_day), type = "response")
  pred_mid <- predict(model, newdata = data.frame(Distance_km = distance_seq, Day_of_year = mid_day), type = "response")
  pred_late <- predict(model, newdata = data.frame(Distance_km = distance_seq, Day_of_year = late_day), type = "response")
  
  color_func <- colorRampPalette(c("lightblue", "navy"))
  
  return(list(
    distance_seq = distance_seq,
    pred_early = pred_early,
    pred_mid = pred_mid,
    pred_late = pred_late,
    early_day = early_day,
    mid_day = mid_day,
    late_day = late_day,
    early_color = color_func(100)[1],
    mid_color = color_func(100)[50],
    late_color = color_func(100)[100]
  ))
}

# ===============================================================
# ANALYZE RECONSTRUCTED TOTAL CONCENTRATION WITH CURVES
# ===============================================================

# Prepare data for reconstructed total
curve_total <- prepare_gamma_data_total(speciation_spatial)

# Fit Gamma GLM model with log link (creates curves)
model_curve_total <- fit_gamma_curve_model(curve_total)

# Create smooth curve predictions
pred_curve_total <- create_curve_predictions(model_curve_total, curve_total)

# ===============================================================
# CREATE PLOT WITH SMOOTH CURVES
# ===============================================================

par(mar = c(4, 4, 3, 1))

# Plot: Reconstructed Total Concentration
plot(curve_total$Distance_km, curve_total$fraction_conc,
     xlab = "Distance (km)", ylab = "Total Cd (mg/L)",
     main = "Reconstructed Total Cd Concentration", pch = 16, col = curve_total$point_color, cex = 1.2,
     ylim = c(0, max(curve_total$fraction_conc) * 1.1))
lines(pred_curve_total$distance_seq, pred_curve_total$pred_early, col = pred_curve_total$early_color, lwd = 3)
lines(pred_curve_total$distance_seq, pred_curve_total$pred_mid, col = pred_curve_total$mid_color, lwd = 3)
lines(pred_curve_total$distance_seq, pred_curve_total$pred_late, col = pred_curve_total$late_color, lwd = 3)

# Add legend
legend("topleft", legend = c("Early Season", "Mid Season", "Late Season"),
       col = c(pred_curve_total$early_color, pred_curve_total$mid_color, pred_curve_total$late_color),
       lwd = 3, cex = 0.8)

# ===============================================================
# ANOVA TEST FOR CURVE MODEL
# ===============================================================

library(car)

print("=== RECONSTRUCTED TOTAL CONCENTRATION ANOVA (GAMMA CURVES) ===")
summary(model_curve_total)
Anova(model_curve_total, type = "II")

# ===============================================================
# INTERPRETATION OF CURVE COEFFICIENTS
# ===============================================================

print("\n===============================================================")
print("CURVE MODEL INTERPRETATION (Log Link):")
print("===============================================================")

# Extract and interpret coefficients
distance_coef <- coef(model_curve_total)["Distance_km"]
percent_change_per_km <- (exp(distance_coef) - 1) * 100

print("Reconstructed Total Cd:")
print(paste("  - Changes by", round(percent_change_per_km, 2), "% per km downstream"))
print(paste("  - Over 168 km:", round((exp(distance_coef * 168) - 1) * 100, 1), "% total change"))

print("===============================================================")


# Deviance-based pseudo R-squared
pseudo_r2_dev <- 1 - (model_curve_total$deviance/model_curve_total$null.deviance)
print(paste("Pseudo R-squared (Deviance):", round(pseudo_r2_dev, 3)))
