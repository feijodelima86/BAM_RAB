# Hypothesis 2: Compartment Interactions Analysis
# Testing whether compartment abundances affect metal concentrations in other compartments
# Based on corrected nested fraction understanding from Hypothesis 1

# Load required libraries
library(dplyr)
library(tidyr)

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

##### GLM #####

# Prepare data for GLM
clean_data <- speciation_spatial[speciation_spatial$Dissolved_pct <= 100, ]

# Convert date to proper format and extract temporal variables
clean_data$Date_parsed <- as.Date(clean_data$Date, format = "%m/%d/%Y")
clean_data$Month <- as.numeric(format(clean_data$Date_parsed, "%m"))
clean_data$Day_of_year <- as.numeric(format(clean_data$Date_parsed, "%j"))

# Ensure suspended_pct is positive for log link
clean_data$Suspended_pct_adj <- pmax(clean_data$Suspended_pct, 0.1)

# Create color mapping based on Day_of_year for points
min_day <- min(clean_data$Day_of_year)
max_day <- max(clean_data$Day_of_year)

# Create color function that maps day of year to blue gradient
color_func <- colorRampPalette(c("lightblue", "navy"))
clean_data$point_color <- color_func(100)[round(((clean_data$Day_of_year - min_day) / (max_day - min_day)) * 99) + 1]

# Fit GLM with interaction
glm_interaction <- glm(Suspended_pct_adj ~ Distance_km + Day_of_year + Distance_km * Day_of_year,
                       data = clean_data,
                       family = gaussian(link = "log"))

# Model summary
summary(glm_interaction)

# Create predictions for different times of year (matching the data range)
distance_seq <- seq(0, 170, by = 10)

# Use actual day range from data
early_day <- min(clean_data$Day_of_year)  # Earliest day in data
mid_day <- median(clean_data$Day_of_year)  # Median day
late_day <- max(clean_data$Day_of_year)   # Latest day in data

# Predictions for early, mid, late season
pred_early <- data.frame(Distance_km = distance_seq, Day_of_year = early_day)
pred_early_link <- predict(glm_interaction, newdata = pred_early, se.fit = TRUE)
pred_early_response <- exp(pred_early_link$fit)

pred_mid <- data.frame(Distance_km = distance_seq, Day_of_year = mid_day)
pred_mid_link <- predict(glm_interaction, newdata = pred_mid, se.fit = TRUE)
pred_mid_response <- exp(pred_mid_link$fit)

pred_late <- data.frame(Distance_km = distance_seq, Day_of_year = late_day)
pred_late_link <- predict(glm_interaction, newdata = pred_late, se.fit = TRUE)
pred_late_response <- exp(pred_late_link$fit)

# Get matching colors for lines
early_color <- color_func(100)[1]    # Lightest blue
mid_color <- color_func(100)[50]     # Medium blue  
late_color <- color_func(100)[100]   # Darkest blue

# Plot with matching colors
plot(clean_data$Distance_km, clean_data$Suspended_pct,
     xlab = "Distance Downstream (km)", 
     ylab = "Suspended Fraction (%)",
     main = "Seasonal Changes in Downstream Cd Precipitation",
     pch = 16, col = clean_data$point_color, cex = 1.5,
     ylim = c(0, max(clean_data$Suspended_pct) * 1.1))

# Add seasonal prediction lines with matching colors
lines(distance_seq, pred_early_response, col = early_color, lwd = 4, lty = 1)
lines(distance_seq, pred_mid_response, col = mid_color, lwd = 4, lty = 1)
lines(distance_seq, pred_late_response, col = late_color, lwd = 4, lty = 1)



# Print results
interaction_effect <- summary(glm_interaction)$coefficients["Distance_km:Day_of_year", ]
print("Interaction Effect (Distance × Day of Year):")
print(interaction_effect)


#### GLM2 ####

# Prepare data for dissolved fraction analysis
clean_data <- speciation_spatial[speciation_spatial$Dissolved_pct <= 100 & speciation_spatial$Dissolved_pct > 0, ]

# Convert date to proper format and extract temporal variables
clean_data$Date_parsed <- as.Date(clean_data$Date, format = "%m/%d/%Y")
clean_data$Month <- as.numeric(format(clean_data$Date_parsed, "%m"))
clean_data$Day_of_year <- as.numeric(format(clean_data$Date_parsed, "%j"))

# Ensure dissolved_pct is positive for log link
clean_data$Dissolved_pct_adj <- pmax(clean_data$Dissolved_pct, 0.1)

# Create color mapping based on Day_of_year for points
min_day <- min(clean_data$Day_of_year)
max_day <- max(clean_data$Day_of_year)

# Create color function that maps day of year to blue gradient
color_func <- colorRampPalette(c("lightblue", "navy"))
clean_data$point_color <- color_func(100)[round(((clean_data$Day_of_year - min_day) / (max_day - min_day)) * 99) + 1]

# Fit GLM with interaction for DISSOLVED fraction
glm_dissolved <- glm(Dissolved_pct_adj ~ Distance_km + Day_of_year + Distance_km * Day_of_year,
                     data = clean_data,
                     family = gaussian(link = "log"))

# Model summary
summary(glm_dissolved)

# Create predictions for different times of year
distance_seq <- seq(0, 170, by = 10)

# Use actual day range from data
early_day <- min(clean_data$Day_of_year)  # Earliest day in data
mid_day <- median(clean_data$Day_of_year)  # Median day
late_day <- max(clean_data$Day_of_year)   # Latest day in data

# Predictions for early, mid, late season
pred_early <- data.frame(Distance_km = distance_seq, Day_of_year = early_day)
pred_early_link <- predict(glm_dissolved, newdata = pred_early, se.fit = TRUE)
pred_early_response <- exp(pred_early_link$fit)

pred_mid <- data.frame(Distance_km = distance_seq, Day_of_year = mid_day)
pred_mid_link <- predict(glm_dissolved, newdata = pred_mid, se.fit = TRUE)
pred_mid_response <- exp(pred_mid_link$fit)

pred_late <- data.frame(Distance_km = distance_seq, Day_of_year = late_day)
pred_late_link <- predict(glm_dissolved, newdata = pred_late, se.fit = TRUE)
pred_late_response <- exp(pred_late_link$fit)

# Get matching colors for lines
early_color <- color_func(100)[1]    # Lightest blue
mid_color <- color_func(100)[50]     # Medium blue  
late_color <- color_func(100)[100]   # Darkest blue

# Plot dissolved fraction with matching colors
plot(clean_data$Distance_km, clean_data$Dissolved_pct,
     xlab = "Distance Downstream (km)", 
     ylab = "Dissolved Fraction (%)",
     main = "Seasonal Changes in Downstream Dissolved Cd",
     pch = 16, col = clean_data$point_color, cex = 1.5,
     ylim = c(0, max(clean_data$Dissolved_pct) * 1.1))

# Add seasonal prediction lines with matching colors
lines(distance_seq, pred_early_response, col = early_color, lwd = 4, lty = 1)
lines(distance_seq, pred_mid_response, col = mid_color, lwd = 4, lty = 1)
lines(distance_seq, pred_late_response, col = late_color, lwd = 4, lty = 1)

# Add legend
legend("topright", 
       legend = c(paste("Early Season (Day", early_day, ")"), 
                  paste("Mid Season (Day", mid_day, ")"),
                  paste("Late Season (Day", late_day, ")")),
       col = c(early_color, mid_color, late_color),
       lwd = 4,
       cex = 0.9)

# Print results
interaction_effect_dissolved <- summary(glm_dissolved)$coefficients["Distance_km:Day_of_year", ]
print("Interaction Effect for DISSOLVED Cd (Distance × Day of Year):")
print(interaction_effect_dissolved)

# Calculate slope changes through season for dissolved fraction
distance_coef <- summary(glm_dissolved)$coefficients["Distance_km", "Estimate"]
interaction_coef <- summary(glm_dissolved)$coefficients["Distance_km:Day_of_year", "Estimate"]

early_slope_dissolved <- distance_coef + interaction_coef * early_day
late_slope_dissolved <- distance_coef + interaction_coef * late_day

print(paste("Early season slope (dissolved):", round(early_slope_dissolved, 6)))
print(paste("Late season slope (dissolved):", round(late_slope_dissolved, 6)))
print(paste("Slope change through season (dissolved):", round(late_slope_dissolved - early_slope_dissolved, 6)))