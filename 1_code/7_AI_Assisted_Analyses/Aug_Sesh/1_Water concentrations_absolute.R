# ===============================================================
# ABSOLUTE CADMIUM CONCENTRATION ANALYSIS - GAMMA GLM (CURVES)
# ===============================================================

# Function to prepare data for Gamma GLM (positive concentrations only)
prepare_gamma_data <- function(data, fraction_name) {
  if (fraction_name == "dissolved") {
    clean_data <- data[!is.na(data$UF_Cd) & data$UF_Cd > 0, ]
    clean_data$fraction_conc <- clean_data$UF_Cd
    clean_data$fraction_name <- "Dissolved"
  } else if (fraction_name == "colloidal") {
    # Calculate colloidal concentration (F_Cd - UF_Cd)
    clean_data <- data[!is.na(data$F_Cd) & !is.na(data$UF_Cd) & 
                         (data$F_Cd - data$UF_Cd) > 0, ]
    clean_data$fraction_conc <- clean_data$F_Cd - clean_data$UF_Cd
    clean_data$fraction_name <- "Colloidal"
  } else if (fraction_name == "suspended") {
    # Calculate suspended concentration (W_Cd - F_Cd)
    clean_data <- data[!is.na(data$W_Cd) & !is.na(data$F_Cd) & 
                         (data$W_Cd - data$F_Cd) > 0, ]
    clean_data$fraction_conc <- clean_data$W_Cd - clean_data$F_Cd
    clean_data$fraction_name <- "Suspended"
  }
  
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
# ANALYZE ALL THREE CONCENTRATIONS WITH CURVES
# ===============================================================

# Prepare data for each fraction
curve_dissolved <- prepare_gamma_data(speciation_spatial, "dissolved")
curve_colloidal <- prepare_gamma_data(speciation_spatial, "colloidal") 
curve_suspended <- prepare_gamma_data(speciation_spatial, "suspended")

# Fit Gamma GLM models with log link (creates curves)
model_curve_dissolved <- fit_gamma_curve_model(curve_dissolved)
model_curve_colloidal <- fit_gamma_curve_model(curve_colloidal)
model_curve_suspended <- fit_gamma_curve_model(curve_suspended)

# Create smooth curve predictions
pred_curve_dissolved <- create_curve_predictions(model_curve_dissolved, curve_dissolved)
pred_curve_colloidal <- create_curve_predictions(model_curve_colloidal, curve_colloidal)
pred_curve_suspended <- create_curve_predictions(model_curve_suspended, curve_suspended)

# ===============================================================
# CREATE SIDE-BY-SIDE PLOTS WITH SMOOTH CURVES
# ===============================================================

par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

# Plot 1: Dissolved Concentration (Exponential Decay Curve)
plot(curve_dissolved$Distance_km, curve_dissolved$fraction_conc,
     xlab = "Distance (km)", ylab = "Dissolved Cd (mg/L)",
     main = "Dissolved Cd Concentration", pch = 16, col = curve_dissolved$point_color, cex = 1.2,
     ylim = c(0, max(curve_dissolved$fraction_conc) * 1.1))
lines(pred_curve_dissolved$distance_seq, pred_curve_dissolved$pred_early, col = pred_curve_dissolved$early_color, lwd = 3)
lines(pred_curve_dissolved$distance_seq, pred_curve_dissolved$pred_mid, col = pred_curve_dissolved$mid_color, lwd = 3)
lines(pred_curve_dissolved$distance_seq, pred_curve_dissolved$pred_late, col = pred_curve_dissolved$late_color, lwd = 3)

# Plot 2: Colloidal Concentration (Growth/Decay Curve)
plot(curve_colloidal$Distance_km, curve_colloidal$fraction_conc,
     xlab = "Distance (km)", ylab = "Colloidal Cd (mg/L)",
     main = "Colloidal Cd Concentration", pch = 16, col = curve_colloidal$point_color, cex = 1.2,
     ylim = c(0, max(curve_colloidal$fraction_conc) * 1.1))
lines(pred_curve_colloidal$distance_seq, pred_curve_colloidal$pred_early, col = pred_curve_colloidal$early_color, lwd = 3)
lines(pred_curve_colloidal$distance_seq, pred_curve_colloidal$pred_mid, col = pred_curve_colloidal$mid_color, lwd = 3)
lines(pred_curve_colloidal$distance_seq, pred_curve_colloidal$pred_late, col = pred_curve_colloidal$late_color, lwd = 3)

# Plot 3: Suspended Concentration (Exponential Growth Curve)
plot(curve_suspended$Distance_km, curve_suspended$fraction_conc,
     xlab = "Distance (km)", ylab = "Suspended Cd (mg/L)",
     main = "Suspended Cd Concentration", pch = 16, col = curve_suspended$point_color, cex = 1.2,
     ylim = c(0, max(curve_suspended$fraction_conc) * 1.1))
lines(pred_curve_suspended$distance_seq, pred_curve_suspended$pred_early, col = pred_curve_suspended$early_color, lwd = 3)
lines(pred_curve_suspended$distance_seq, pred_curve_suspended$pred_mid, col = pred_curve_suspended$mid_color, lwd = 3)
lines(pred_curve_suspended$distance_seq, pred_curve_suspended$pred_late, col = pred_curve_suspended$late_color, lwd = 3)

# Add legend to the last plot
legend("topleft", legend = c("Early Season", "Mid Season", "Late Season"),
       col = c(pred_curve_suspended$early_color, pred_curve_suspended$mid_color, pred_curve_suspended$late_color),
       lwd = 3, cex = 0.8)

par(mfrow = c(1, 1))

# ===============================================================
# ANOVA TESTS FOR CURVE MODELS
# ===============================================================

library(car)

print("=== DISSOLVED CONCENTRATION ANOVA (GAMMA CURVES) ===")
summary(model_curve_dissolved)
Anova(model_curve_dissolved, type = "II")

print("\n=== COLLOIDAL CONCENTRATION ANOVA (GAMMA CURVES) ===")
summary(model_curve_colloidal)
Anova(model_curve_colloidal, type = "II")

print("\n=== SUSPENDED CONCENTRATION ANOVA (GAMMA CURVES) ===")
summary(model_curve_suspended)
Anova(model_curve_suspended, type = "II")

# ===============================================================
# INTERPRETATION OF CURVE COEFFICIENTS
# ===============================================================

print("\n===============================================================")
print("CURVE MODEL INTERPRETATION (Log Link):")
print("===============================================================")

# Extract and interpret coefficients
for(model_name in c("Dissolved", "Colloidal", "Suspended")) {
  if(model_name == "Dissolved") model <- model_curve_dissolved
  if(model_name == "Colloidal") model <- model_curve_colloidal  
  if(model_name == "Suspended") model <- model_curve_suspended
  
  distance_coef <- coef(model)["Distance_km"]
  percent_change_per_km <- (exp(distance_coef) - 1) * 100
  
  print(paste(model_name, "Cd:"))
  print(paste("  - Changes by", round(percent_change_per_km, 2), "% per km downstream"))
  print(paste("  - Over 168 km:", round((exp(distance_coef * 168) - 1) * 100, 1), "% total change"))
}

print("===============================================================")