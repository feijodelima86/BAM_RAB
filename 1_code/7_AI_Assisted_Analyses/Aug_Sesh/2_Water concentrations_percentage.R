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

# Prepare data - filter out unrealistic values and ensure suspended fraction exists
clean_data_suspended <- speciation_spatial[speciation_spatial$Suspended_pct >= 0 & 
                                             speciation_spatial$Suspended_pct <= 100 &
                                             !is.na(speciation_spatial$Suspended_pct), ]

# Function to prepare data for each fraction
prepare_fraction_data <- function(data, fraction_name) {
  if (fraction_name == "dissolved") {
    clean_data <- data[data$Dissolved_pct <= 100 & data$Dissolved_pct > 5, ]
    clean_data$fraction_pct <- clean_data$Dissolved_pct
  } else if (fraction_name == "colloidal") {
    clean_data <- data[data$Colloidal_pct >= 0 & data$Colloidal_pct <= 100 & 
                         data$Colloidal_present == TRUE, ]
    clean_data$fraction_pct <- clean_data$Colloidal_pct
  } else if (fraction_name == "suspended") {
    clean_data <- data[data$Suspended_pct >= 0 & data$Suspended_pct <= 100, ]
    clean_data$fraction_pct <- clean_data$Suspended_pct
  }
  
  # Convert date variables
  clean_data$Date_parsed <- as.Date(clean_data$Date, format = "%m/%d/%Y")
  clean_data$Day_of_year <- as.numeric(format(clean_data$Date_parsed, "%j"))
  
  # Convert to proportions for Beta regression
  clean_data$fraction_prop <- clean_data$fraction_pct / 100
  epsilon <- 0.001
  clean_data$fraction_prop_adj <- pmax(pmin(clean_data$fraction_prop, 1-epsilon), epsilon)
  
  # Create colors
  min_day <- min(clean_data$Day_of_year)
  max_day <- max(clean_data$Day_of_year)
  color_func <- colorRampPalette(c("lightblue", "navy"))
  clean_data$point_color <- color_func(100)[round(((clean_data$Day_of_year - min_day) / (max_day - min_day)) * 99) + 1]
  
  return(clean_data)
}

# Function to fit Beta regression model
fit_beta_model <- function(data) {
  betareg(fraction_prop_adj ~ Distance_km + Day_of_year + Distance_km * Day_of_year, data = data)
}

# Function to create predictions
create_predictions <- function(model, data) {
  distance_seq <- seq(0, 170, by = 10)
  early_day <- min(data$Day_of_year)
  mid_day <- median(data$Day_of_year)
  late_day <- max(data$Day_of_year)
  
  pred_early <- predict(model, newdata = data.frame(Distance_km = distance_seq, Day_of_year = early_day), type = "response") * 100
  pred_mid <- predict(model, newdata = data.frame(Distance_km = distance_seq, Day_of_year = mid_day), type = "response") * 100
  pred_late <- predict(model, newdata = data.frame(Distance_km = distance_seq, Day_of_year = late_day), type = "response") * 100
  
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
# ANALYZE ALL THREE FRACTIONS
# ===============================================================

# Prepare data for each fraction
clean_dissolved <- prepare_fraction_data(speciation_spatial, "dissolved")
clean_colloidal <- prepare_fraction_data(speciation_spatial, "colloidal") 
clean_suspended <- prepare_fraction_data(speciation_spatial, "suspended")

# Fit models
model_dissolved <- fit_beta_model(clean_dissolved)
model_colloidal <- fit_beta_model(clean_colloidal)
model_suspended <- fit_beta_model(clean_suspended)

# Model Summaries
summary(model_dissolved)
Anova(model_dissolved, type = "II")

summary(model_colloidal)
Anova(model_colloidal, type = "II")

summary(model_suspended)
Anova(model_suspended, type = "II")

# Create predictions
pred_dissolved <- create_predictions(model_dissolved, clean_dissolved)
pred_colloidal <- create_predictions(model_colloidal, clean_colloidal)
pred_suspended <- create_predictions(model_suspended, clean_suspended)

# ===============================================================
# CREATE SIDE-BY-SIDE PLOTS
# ===============================================================

par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

# Plot 1: Dissolved Fraction
plot(clean_dissolved$Distance_km, clean_dissolved$fraction_pct,
     xlab = "Distance (km)", ylab = "Dissolved Fraction (%)",
     main = "Dissolved Cd", pch = 16, col = clean_dissolved$point_color, cex = 1.2,
     ylim = c(0, 100))
lines(pred_dissolved$distance_seq, pred_dissolved$pred_early, col = pred_dissolved$early_color, lwd = 3)
lines(pred_dissolved$distance_seq, pred_dissolved$pred_mid, col = pred_dissolved$mid_color, lwd = 3)
lines(pred_dissolved$distance_seq, pred_dissolved$pred_late, col = pred_dissolved$late_color, lwd = 3)

# Plot 2: Colloidal Fraction
plot(clean_colloidal$Distance_km, clean_colloidal$fraction_pct,
     xlab = "Distance (km)", ylab = "Colloidal Fraction (%)",
     main = "Colloidal Cd", pch = 16, col = clean_colloidal$point_color, cex = 1.2,
     ylim = c(0, 100))
lines(pred_colloidal$distance_seq, pred_colloidal$pred_early, col = pred_colloidal$early_color, lwd = 3)
lines(pred_colloidal$distance_seq, pred_colloidal$pred_mid, col = pred_colloidal$mid_color, lwd = 3)
lines(pred_colloidal$distance_seq, pred_colloidal$pred_late, col = pred_colloidal$late_color, lwd = 3)

# Plot 3: Suspended Fraction
plot(clean_suspended$Distance_km, clean_suspended$fraction_pct,
     xlab = "Distance (km)", ylab = "Suspended Fraction (%)",
     main = "Suspended Cd", pch = 16, col = clean_suspended$point_color, cex = 1.2,
     ylim = c(0, 100))
lines(pred_suspended$distance_seq, pred_suspended$pred_early, col = pred_suspended$early_color, lwd = 3)
lines(pred_suspended$distance_seq, pred_suspended$pred_mid, col = pred_suspended$mid_color, lwd = 3)
lines(pred_suspended$distance_seq, pred_suspended$pred_late, col = pred_suspended$late_color, lwd = 3)

# Add legend to the last plot
legend("topleft", legend = c("Early Season", "Mid Season", "Late Season"),
       col = c(pred_suspended$early_color, pred_suspended$mid_color, pred_suspended$late_color),
       lwd = 3, cex = 0.8)

par(mfrow = c(1, 1))







