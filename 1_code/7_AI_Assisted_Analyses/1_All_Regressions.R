# Linear Regression Plots: EPIL vs EPIP Metal Concentrations
# 3x4 panel layout for all 12 elements using ALL individual samples (no averaging)

# Load required libraries
library(dplyr)
library(tidyr)

# Read data
metal_data <- read.csv("2_incremental/7_AI_ANALYSIS/metal_analysis_subset.csv")

# Define metal configurations
metal_configs <- list(
  As = list(tissue = "As", name = "Arsenic"),
  Cd = list(tissue = "Cd", name = "Cadmium"),
  Co = list(tissue = "Co", name = "Cobalt"),
  Cr = list(tissue = "Cr", name = "Chromium"),
  Cu = list(tissue = "Cu", name = "Copper"),
  Fe = list(tissue = "Fer.1", name = "Iron"),
  Mn = list(tissue = "Mn", name = "Manganese"),
  Mo = list(tissue = "Mo", name = "Molybdenum"),
  Ni = list(tissue = "Ni", name = "Nickel"),
  Pb = list(tissue = "Pb", name = "Lead"),
  Se = list(tissue = "Se", name = "Selenium"),
  Zn = list(tissue = "Zn", name = "Zinc")
)

# Function to prepare data for EPIL vs EPIP comparison using ALL individual samples
prepare_epil_epip_data <- function(data, tissue_col) {
  
  # Clean and filter data - keep all individual samples
  clean_data <- data %>%
    filter(!is.na(.data[[tissue_col]]) & 
             SAMPLE_DESCRIPTOR %in% c("EPIL", "EPIP") &
             !is.na(SITE) & !is.na(SAMPLING_DATE) &
             is.finite(.data[[tissue_col]]))
  
  if(nrow(clean_data) == 0) return(NULL)
  
  # Create unique sample ID for tracking
  clean_data <- clean_data %>%
    mutate(sample_id = paste(SITE, SAMPLING_DATE, SAMPLE_DESCRIPTOR, row_number(), sep = "_")) %>%
    select(SITE, SAMPLING_DATE, SAMPLE_DESCRIPTOR, sample_id, metal_conc = all_of(tissue_col))
  
  # Separate EPIL and EPIP samples
  epil_data <- clean_data %>%
    filter(SAMPLE_DESCRIPTOR == "EPIL") %>%
    select(SITE, SAMPLING_DATE, sample_id, epil_conc = metal_conc)
  
  epip_data <- clean_data %>%
    filter(SAMPLE_DESCRIPTOR == "EPIP") %>%
    select(SITE, SAMPLING_DATE, sample_id, epip_conc = metal_conc)
  
  # For each site-date combination, create all possible EPIL-EPIP pairs
  site_dates <- unique(clean_data[, c("SITE", "SAMPLING_DATE")])
  
  all_pairs <- data.frame()
  
  for(i in 1:nrow(site_dates)) {
    site <- site_dates$SITE[i]
    date <- site_dates$SAMPLING_DATE[i]
    
    epil_subset <- epil_data %>% filter(SITE == site & SAMPLING_DATE == date)
    epip_subset <- epip_data %>% filter(SITE == site & SAMPLING_DATE == date)
    
    if(nrow(epil_subset) > 0 & nrow(epip_subset) > 0) {
      # Create all combinations of EPIL and EPIP samples for this site-date
      for(j in 1:nrow(epil_subset)) {
        for(k in 1:nrow(epip_subset)) {
          all_pairs <- rbind(all_pairs, data.frame(
            SITE = site,
            SAMPLING_DATE = date,
            epil_sample_id = epil_subset$sample_id[j],
            epip_sample_id = epip_subset$sample_id[k],
            conc_EPIL = epil_subset$epil_conc[j],
            conc_EPIP = epip_subset$epip_conc[k]
          ))
        }
      }
    }
  }
  
  return(all_pairs)
}

# Function to create individual regression plot with automatic outlier detection
create_regression_plot <- function(data, metal_symbol, metal_name) {
  
  if(is.null(data) || nrow(data) < 10) {  # Need more points for outlier detection
    # Create empty plot for insufficient data
    plot(1, 1, type = "n", xlab = "", ylab = "", main = paste(metal_symbol, "-", metal_name),
         axes = FALSE)
    text(1, 1, "Insufficient\nData", cex = 1.2)
    return()
  }
  
  # Calculate initial regression to identify outliers
  lm_fit_initial <- lm(conc_EPIP ~ conc_EPIL, data = data)
  
  # Calculate Cook's distance to identify influential outliers
  cooks_dist <- cooks.distance(lm_fit_initial)
  
  # Define outliers as points with Cook's distance > 4/n (common threshold)
  n_points <- nrow(data)
  cooks_threshold <- 4 / n_points
  outliers <- cooks_dist > cooks_threshold
  
  # Also identify extreme outliers using standardized residuals
  std_residuals <- rstandard(lm_fit_initial)
  extreme_outliers <- abs(std_residuals) > 3  # Points beyond 3 standard deviations
  
  # Combine both criteria
  all_outliers <- outliers | extreme_outliers
  
  # Calculate statistics for both full and cleaned datasets
  # Full dataset
  cor_test_full <- cor.test(data$conc_EPIL, data$conc_EPIP, method = "pearson")
  r_full <- round(cor_test_full$estimate, 3)
  r2_full <- round(summary(lm_fit_initial)$r.squared, 3)
  p_full <- cor_test_full$p.value
  
  # Dataset without outliers (if we have enough points left)
  if(sum(!all_outliers) >= 5) {
    data_clean <- data[!all_outliers, ]
    cor_test_clean <- cor.test(data_clean$conc_EPIL, data_clean$conc_EPIP, method = "pearson")
    lm_fit_clean <- lm(conc_EPIP ~ conc_EPIL, data = data_clean)
    
    r_clean <- round(cor_test_clean$estimate, 3)
    r2_clean <- round(summary(lm_fit_clean)$r.squared, 3)
    p_clean <- cor_test_clean$p.value
    n_clean <- nrow(data_clean)
    n_outliers <- sum(all_outliers)
  } else {
    # Too many outliers removed, use full dataset
    data_clean <- data
    lm_fit_clean <- lm_fit_initial
    cor_test_clean <- cor_test_full
    
    r_clean <- r_full
    r2_clean <- r2_full
    p_clean <- p_full
    n_clean <- n_points
    n_outliers <- 0
    all_outliers <- rep(FALSE, n_points)
  }
  
  # Format p-values
  p_formatted_full <- ifelse(p_full < 0.001, "< 0.001", paste("=", round(p_full, 3)))
  p_formatted_clean <- ifelse(p_clean < 0.001, "< 0.001", paste("=", round(p_clean, 3)))
  
  # Create the plot
  # First plot normal points in blue
  plot(data$conc_EPIL[!all_outliers], data$conc_EPIP[!all_outliers],
       xlab = "EPIL Concentration (mg/g)",
       ylab = "EPIP Concentration (mg/g)",
       main = paste(metal_symbol, "-", metal_name),
       pch = 16, col = "steelblue", cex = 0.8,
       cex.main = 1.1, cex.lab = 1.0, cex.axis = 0.9,
       xlim = range(data$conc_EPIL), ylim = range(data$conc_EPIP))
  
  # Add outliers in gray
  if(any(all_outliers)) {
    points(data$conc_EPIL[all_outliers], data$conc_EPIP[all_outliers],
           pch = 16, col = "gray60", cex = 0.8)
  }
  
  # Add regression line using CLEANED data (without outliers)
  abline(lm_fit_clean, col = "red", lwd = 2)
  
  # Add confidence interval for cleaned regression
  x_seq <- seq(min(data$conc_EPIL), max(data$conc_EPIL), length.out = 100)
  pred_interval <- predict(lm_fit_clean, newdata = data.frame(conc_EPIL = x_seq), 
                           interval = "confidence", level = 0.95)
  
  # Add confidence bands
  lines(x_seq, pred_interval[,"lwr"], col = "red", lty = 2, lwd = 1)
  lines(x_seq, pred_interval[,"upr"], col = "red", lty = 2, lwd = 1)
  
  # Create statistics text
  if(n_outliers > 0) {
    stats_text <- paste0(
      "All data (n=", n_points, "):\n",
      "r = ", r_full, ", R² = ", r2_full, "\n",
      "p ", p_formatted_full, "\n\n",
      "Clean data (n=", n_clean, "):\n",
      "r = ", r_clean, ", R² = ", r2_clean, "\n",
      "p ", p_formatted_clean, "\n",
      "Outliers removed: ", n_outliers
    )
  } else {
    stats_text <- paste0(
      "r = ", r_clean, "\n",
      "R² = ", r2_clean, "\n",
      "p ", p_formatted_clean, "\n",
      "n = ", n_clean, "\n",
      "No outliers detected"
    )
  }
  
  # Position text in upper left corner
  x_range <- range(data$conc_EPIL)
  y_range <- range(data$conc_EPIP)
  
  text(x_range[1] + 0.05 * diff(x_range), 
       y_range[2] - 0.05 * diff(y_range),
       stats_text, 
       adj = c(0, 1), 
       cex = 0.7, 
       col = "black",
       font = 1)
  
  # Add legend if outliers were detected
  if(any(all_outliers)) {
    legend("bottomright", 
           legend = c("Normal data", "Outliers", "Regression (clean)"),
           col = c("steelblue", "gray60", "red"),
           pch = c(16, 16, NA),
           lty = c(NA, NA, 1),
           lwd = c(NA, NA, 2),
           cex = 0.7,
           bg = "white")
  }
}

# Set up the plotting layout
cat("Generating EPIL vs EPIP regression plots using ALL individual samples...\n")

# Create PNG file for high-quality output
png("EPIL_EPIP_individual_samples_regression.png", width = 16, height = 12, units = "in", res = 300)

# Set up 3x4 panel layout with margins
par(mfrow = c(3, 4), 
    mar = c(4, 4, 3, 1),  # bottom, left, top, right margins
    oma = c(2, 2, 4, 2))  # outer margins

# Generate plots for all metals
for(i in 1:length(metal_configs)) {
  metal_symbol <- names(metal_configs)[i]
  config <- metal_configs[[i]]
  
  cat("Processing", metal_symbol, "...\n")
  
  # Check if tissue column exists
  if(!config$tissue %in% names(metal_data)) {
    cat("Warning: Column", config$tissue, "not found for", metal_symbol, "\n")
    create_regression_plot(NULL, metal_symbol, config$name)
    next
  }
  
  # Prepare data
  plot_data <- prepare_epil_epip_data(metal_data, config$tissue)
  
  # Create plot
  create_regression_plot(plot_data, metal_symbol, config$name)
}

# Add overall title
mtext("EPIL vs EPIP Metal Concentrations - Individual Sample Pairs", 
      outer = TRUE, cex = 1.5, font = 2, line = 1)

# Add footnote
mtext("Blue points = normal data, Gray points = outliers (excluded from regression). Red line = robust regression with 95% CI.", 
      outer = TRUE, cex = 0.9, line = 0, side = 1)

dev.off()

cat("Plot saved as 'EPIL_EPIP_individual_samples_regression.png'\n")

# Also create a summary table of regression statistics
cat("\n=== INDIVIDUAL SAMPLE REGRESSION STATISTICS ===\n")

regression_summary <- data.frame(
  Metal = character(),
  N_pairs = numeric(),
  N_clean = numeric(),
  N_outliers = numeric(),
  R = numeric(),
  R_squared = numeric(),
  P_value = numeric(),
  Slope = numeric(),
  Intercept = numeric(),
  stringsAsFactors = FALSE
)

for(i in 1:length(metal_configs)) {
  metal_symbol <- names(metal_configs)[i]
  config <- metal_configs[[i]]
  
  if(!config$tissue %in% names(metal_data)) next
  
  plot_data <- prepare_epil_epip_data(metal_data, config$tissue)
  
  if(is.null(plot_data) || nrow(plot_data) < 10) {
    regression_summary <- rbind(regression_summary, data.frame(
      Metal = metal_symbol,
      N_pairs = ifelse(is.null(plot_data), 0, nrow(plot_data)),
      N_clean = NA,
      N_outliers = NA,
      R = NA,
      R_squared = NA,
      P_value = NA,
      Slope = NA,
      Intercept = NA
    ))
    next
  }
  
  # Recalculate outlier detection for summary
  lm_initial <- lm(conc_EPIP ~ conc_EPIL, data = plot_data)
  cooks_dist <- cooks.distance(lm_initial)
  std_residuals <- rstandard(lm_initial)
  
  cooks_threshold <- 4 / nrow(plot_data)
  outliers_detected <- (cooks_dist > cooks_threshold) | (abs(std_residuals) > 3)
  
  # Calculate final statistics using cleaned data
  if(sum(!outliers_detected) >= 5) {
    clean_data <- plot_data[!outliers_detected, ]
    cor_test <- cor.test(clean_data$conc_EPIL, clean_data$conc_EPIP, method = "pearson")
    lm_fit <- lm(conc_EPIP ~ conc_EPIL, data = clean_data)
  } else {
    # Use full dataset if too many outliers
    cor_test <- cor.test(plot_data$conc_EPIL, plot_data$conc_EPIP, method = "pearson")
    lm_fit <- lm(conc_EPIP ~ conc_EPIL, data = plot_data)
    outliers_detected <- rep(FALSE, nrow(plot_data))
  }
  
  regression_summary <- rbind(regression_summary, data.frame(
    Metal = metal_symbol,
    N_pairs = nrow(plot_data),
    R = round(cor_test$estimate, 3),
    R_squared = round(summary(lm_fit)$r.squared, 3),
    P_value = round(cor_test$p.value, 4),
    Slope = round(coef(lm_fit)[2], 4),
    Intercept = round(coef(lm_fit)[1], 4)
  ))
}

print(regression_summary, row.names = FALSE)

# Save summary table
write.csv(regression_summary, "EPIL_EPIP_individual_samples_summary.csv", row.names = FALSE)

cat("\nRegression summary saved as 'EPIL_EPIP_individual_samples_summary.csv'\n")

# Compare sample sizes and show data structure
cat("\n=== DATA STRUCTURE SUMMARY ===\n")
for(i in 1:min(3, length(metal_configs))) {  # Show first 3 metals as examples
  metal_symbol <- names(metal_configs)[i]
  config <- metal_configs[[i]]
  
  if(!config$tissue %in% names(metal_data)) next
  
  plot_data <- prepare_epil_epip_data(metal_data, config$tissue)
  if(is.null(plot_data)) next
  
  n_sites <- length(unique(paste(plot_data$SITE, plot_data$SAMPLING_DATE)))
  n_pairs <- nrow(plot_data)
  
  cat(sprintf("%s: %d site-date combinations → %d individual sample pairs\n", 
              metal_symbol, n_sites, n_pairs))
}

cat("\nNote: Each site-date combination creates multiple pairs if there are multiple samples per compartment.\n")

# Linear Regression Plots: FILA vs EPIP Metal Concentrations
# 3x4 panel layout for all 12 elements using ALL individual samples (no averaging)

# Load required libraries
library(dplyr)
library(tidyr)

# Read data
metal_data <- read.csv("2_incremental/7_AI_ANALYSIS/metal_analysis_subset.csv")

# Define metal configurations
metal_configs <- list(
  As = list(tissue = "As", name = "Arsenic"),
  Cd = list(tissue = "Cd", name = "Cadmium"),
  Co = list(tissue = "Co", name = "Cobalt"),
  Cr = list(tissue = "Cr", name = "Chromium"),
  Cu = list(tissue = "Cu", name = "Copper"),
  Fe = list(tissue = "Fer.1", name = "Iron"),
  Mn = list(tissue = "Mn", name = "Manganese"),
  Mo = list(tissue = "Mo", name = "Molybdenum"),
  Ni = list(tissue = "Ni", name = "Nickel"),
  Pb = list(tissue = "Pb", name = "Lead"),
  Se = list(tissue = "Se", name = "Selenium"),
  Zn = list(tissue = "Zn", name = "Zinc")
)

# Function to prepare data for FILA vs EPIP comparison using ALL individual samples
prepare_FILA_epip_data <- function(data, tissue_col) {
  
  # Clean and filter data - keep all individual samples
  clean_data <- data %>%
    filter(!is.na(.data[[tissue_col]]) & 
             SAMPLE_DESCRIPTOR %in% c("FILA", "EPIP") &
             !is.na(SITE) & !is.na(SAMPLING_DATE) &
             is.finite(.data[[tissue_col]]))
  
  if(nrow(clean_data) == 0) return(NULL)
  
  # Create unique sample ID for tracking
  clean_data <- clean_data %>%
    mutate(sample_id = paste(SITE, SAMPLING_DATE, SAMPLE_DESCRIPTOR, row_number(), sep = "_")) %>%
    select(SITE, SAMPLING_DATE, SAMPLE_DESCRIPTOR, sample_id, metal_conc = all_of(tissue_col))
  
  # Separate FILA and EPIP samples
  FILA_data <- clean_data %>%
    filter(SAMPLE_DESCRIPTOR == "FILA") %>%
    select(SITE, SAMPLING_DATE, sample_id, FILA_conc = metal_conc)
  
  epip_data <- clean_data %>%
    filter(SAMPLE_DESCRIPTOR == "EPIP") %>%
    select(SITE, SAMPLING_DATE, sample_id, epip_conc = metal_conc)
  
  # For each site-date combination, create all possible FILA-EPIP pairs
  site_dates <- unique(clean_data[, c("SITE", "SAMPLING_DATE")])
  
  all_pairs <- data.frame()
  
  for(i in 1:nrow(site_dates)) {
    site <- site_dates$SITE[i]
    date <- site_dates$SAMPLING_DATE[i]
    
    FILA_subset <- FILA_data %>% filter(SITE == site & SAMPLING_DATE == date)
    epip_subset <- epip_data %>% filter(SITE == site & SAMPLING_DATE == date)
    
    if(nrow(FILA_subset) > 0 & nrow(epip_subset) > 0) {
      # Create all combinations of FILA and EPIP samples for this site-date
      for(j in 1:nrow(FILA_subset)) {
        for(k in 1:nrow(epip_subset)) {
          all_pairs <- rbind(all_pairs, data.frame(
            SITE = site,
            SAMPLING_DATE = date,
            FILA_sample_id = FILA_subset$sample_id[j],
            epip_sample_id = epip_subset$sample_id[k],
            conc_FILA = FILA_subset$FILA_conc[j],
            conc_EPIP = epip_subset$epip_conc[k]
          ))
        }
      }
    }
  }
  
  return(all_pairs)
}

# Function to create individual regression plot with automatic outlier detection
create_regression_plot <- function(data, metal_symbol, metal_name) {
  
  if(is.null(data) || nrow(data) < 10) {  # Need more points for outlier detection
    # Create empty plot for insufficient data
    plot(1, 1, type = "n", xlab = "", ylab = "", main = paste(metal_symbol, "-", metal_name),
         axes = FALSE)
    text(1, 1, "Insufficient\nData", cex = 1.2)
    return()
  }
  
  # Calculate initial regression to identify outliers
  lm_fit_initial <- lm(conc_EPIP ~ conc_FILA, data = data)
  
  # Calculate Cook's distance to identify influential outliers
  cooks_dist <- cooks.distance(lm_fit_initial)
  
  # Define outliers as points with Cook's distance > 4/n (common threshold)
  n_points <- nrow(data)
  cooks_threshold <- 4 / n_points
  outliers <- cooks_dist > cooks_threshold
  
  # Also identify extreme outliers using standardized residuals
  std_residuals <- rstandard(lm_fit_initial)
  extreme_outliers <- abs(std_residuals) > 3  # Points beyond 3 standard deviations
  
  # Combine both criteria
  all_outliers <- outliers | extreme_outliers
  
  # Calculate statistics for both full and cleaned datasets
  # Full dataset
  cor_test_full <- cor.test(data$conc_FILA, data$conc_EPIP, method = "pearson")
  r_full <- round(cor_test_full$estimate, 3)
  r2_full <- round(summary(lm_fit_initial)$r.squared, 3)
  p_full <- cor_test_full$p.value
  
  # Dataset without outliers (if we have enough points left)
  if(sum(!all_outliers) >= 5) {
    data_clean <- data[!all_outliers, ]
    cor_test_clean <- cor.test(data_clean$conc_FILA, data_clean$conc_EPIP, method = "pearson")
    lm_fit_clean <- lm(conc_EPIP ~ conc_FILA, data = data_clean)
    
    r_clean <- round(cor_test_clean$estimate, 3)
    r2_clean <- round(summary(lm_fit_clean)$r.squared, 3)
    p_clean <- cor_test_clean$p.value
    n_clean <- nrow(data_clean)
    n_outliers <- sum(all_outliers)
  } else {
    # Too many outliers removed, use full dataset
    data_clean <- data
    lm_fit_clean <- lm_fit_initial
    cor_test_clean <- cor_test_full
    
    r_clean <- r_full
    r2_clean <- r2_full
    p_clean <- p_full
    n_clean <- n_points
    n_outliers <- 0
    all_outliers <- rep(FALSE, n_points)
  }
  
  # Format p-values
  p_formatted_full <- ifelse(p_full < 0.001, "< 0.001", paste("=", round(p_full, 3)))
  p_formatted_clean <- ifelse(p_clean < 0.001, "< 0.001", paste("=", round(p_clean, 3)))
  
  # Create the plot
  # First plot normal points in blue
  plot(data$conc_FILA[!all_outliers], data$conc_EPIP[!all_outliers],
       xlab = "FILA Concentration (mg/g)",
       ylab = "EPIP Concentration (mg/g)",
       main = paste(metal_symbol, "-", metal_name),
       pch = 16, col = "steelblue", cex = 0.8,
       cex.main = 1.1, cex.lab = 1.0, cex.axis = 0.9,
       xlim = range(data$conc_FILA), ylim = range(data$conc_EPIP))
  
  # Add outliers in gray
  if(any(all_outliers)) {
    points(data$conc_FILA[all_outliers], data$conc_EPIP[all_outliers],
           pch = 16, col = "gray60", cex = 0.8)
  }
  
  # Add regression line using CLEANED data (without outliers)
  abline(lm_fit_clean, col = "red", lwd = 2)
  
  # Add confidence interval for cleaned regression
  x_seq <- seq(min(data$conc_FILA), max(data$conc_FILA), length.out = 100)
  pred_interval <- predict(lm_fit_clean, newdata = data.frame(conc_FILA = x_seq), 
                           interval = "confidence", level = 0.95)
  
  # Add confidence bands
  lines(x_seq, pred_interval[,"lwr"], col = "red", lty = 2, lwd = 1)
  lines(x_seq, pred_interval[,"upr"], col = "red", lty = 2, lwd = 1)
  
  # Create statistics text
  if(n_outliers > 0) {
    stats_text <- paste0(
      "All data (n=", n_points, "):\n",
      "r = ", r_full, ", R² = ", r2_full, "\n",
      "p ", p_formatted_full, "\n\n",
      "Clean data (n=", n_clean, "):\n",
      "r = ", r_clean, ", R² = ", r2_clean, "\n",
      "p ", p_formatted_clean, "\n",
      "Outliers removed: ", n_outliers
    )
  } else {
    stats_text <- paste0(
      "r = ", r_clean, "\n",
      "R² = ", r2_clean, "\n",
      "p ", p_formatted_clean, "\n",
      "n = ", n_clean, "\n",
      "No outliers detected"
    )
  }
  
  # Position text in upper left corner
  x_range <- range(data$conc_FILA)
  y_range <- range(data$conc_EPIP)
  
  text(x_range[1] + 0.05 * diff(x_range), 
       y_range[2] - 0.05 * diff(y_range),
       stats_text, 
       adj = c(0, 1), 
       cex = 0.7, 
       col = "black",
       font = 1)
  
  # Add legend if outliers were detected
  if(any(all_outliers)) {
    legend("bottomright", 
           legend = c("Normal data", "Outliers", "Regression (clean)"),
           col = c("steelblue", "gray60", "red"),
           pch = c(16, 16, NA),
           lty = c(NA, NA, 1),
           lwd = c(NA, NA, 2),
           cex = 0.7,
           bg = "white")
  }
}

# Set up the plotting layout
cat("Generating FILA vs EPIP regression plots using ALL individual samples...\n")

# Create PNG file for high-quality output
png("FILA_EPIP_individual_samples_regression.png", width = 16, height = 12, units = "in", res = 300)

# Set up 3x4 panel layout with margins
par(mfrow = c(3, 4), 
    mar = c(4, 4, 3, 1),  # bottom, left, top, right margins
    oma = c(2, 2, 4, 2))  # outer margins

# Generate plots for all metals
for(i in 1:length(metal_configs)) {
  metal_symbol <- names(metal_configs)[i]
  config <- metal_configs[[i]]
  
  cat("Processing", metal_symbol, "...\n")
  
  # Check if tissue column exists
  if(!config$tissue %in% names(metal_data)) {
    cat("Warning: Column", config$tissue, "not found for", metal_symbol, "\n")
    create_regression_plot(NULL, metal_symbol, config$name)
    next
  }
  
  # Prepare data
  plot_data <- prepare_FILA_epip_data(metal_data, config$tissue)
  
  # Create plot
  create_regression_plot(plot_data, metal_symbol, config$name)
}

# Add overall title
mtext("FILA vs EPIP Metal Concentrations - Individual Sample Pairs", 
      outer = TRUE, cex = 1.5, font = 2, line = 1)

# Add footnote
mtext("Blue points = normal data, Gray points = outliers (excluded from regression). Red line = robust regression with 95% CI.", 
      outer = TRUE, cex = 0.9, line = 0, side = 1)

dev.off()

cat("Plot saved as 'FILA_EPIP_individual_samples_regression.png'\n")

# Also create a summary table of regression statistics
cat("\n=== INDIVIDUAL SAMPLE REGRESSION STATISTICS ===\n")

regression_summary <- data.frame(
  Metal = character(),
  N_pairs = numeric(),
  N_clean = numeric(),
  N_outliers = numeric(),
  R = numeric(),
  R_squared = numeric(),
  P_value = numeric(),
  Slope = numeric(),
  Intercept = numeric(),
  stringsAsFactors = FALSE
)

for(i in 1:length(metal_configs)) {
  metal_symbol <- names(metal_configs)[i]
  config <- metal_configs[[i]]
  
  if(!config$tissue %in% names(metal_data)) next
  
  plot_data <- prepare_FILA_epip_data(metal_data, config$tissue)
  
  if(is.null(plot_data) || nrow(plot_data) < 10) {
    regression_summary <- rbind(regression_summary, data.frame(
      Metal = metal_symbol,
      N_pairs = ifelse(is.null(plot_data), 0, nrow(plot_data)),
      N_clean = NA,
      N_outliers = NA,
      R = NA,
      R_squared = NA,
      P_value = NA,
      Slope = NA,
      Intercept = NA
    ))
    next
  }
  
  # Recalculate outlier detection for summary
  lm_initial <- lm(conc_EPIP ~ conc_FILA, data = plot_data)
  cooks_dist <- cooks.distance(lm_initial)
  std_residuals <- rstandard(lm_initial)
  
  cooks_threshold <- 4 / nrow(plot_data)
  outliers_detected <- (cooks_dist > cooks_threshold) | (abs(std_residuals) > 3)
  
  # Calculate final statistics using cleaned data
  if(sum(!outliers_detected) >= 5) {
    clean_data <- plot_data[!outliers_detected, ]
    cor_test <- cor.test(clean_data$conc_FILA, clean_data$conc_EPIP, method = "pearson")
    lm_fit <- lm(conc_EPIP ~ conc_FILA, data = clean_data)
  } else {
    # Use full dataset if too many outliers
    cor_test <- cor.test(plot_data$conc_FILA, plot_data$conc_EPIP, method = "pearson")
    lm_fit <- lm(conc_EPIP ~ conc_FILA, data = plot_data)
    outliers_detected <- rep(FALSE, nrow(plot_data))
  }
  
  regression_summary <- rbind(regression_summary, data.frame(
    Metal = metal_symbol,
    N_pairs = nrow(plot_data),
    R = round(cor_test$estimate, 3),
    R_squared = round(summary(lm_fit)$r.squared, 3),
    P_value = round(cor_test$p.value, 4),
    Slope = round(coef(lm_fit)[2], 4),
    Intercept = round(coef(lm_fit)[1], 4)
  ))
}

print(regression_summary, row.names = FALSE)

# Save summary table
write.csv(regression_summary, "FILA_EPIP_individual_samples_summary.csv", row.names = FALSE)

cat("\nRegression summary saved as 'FILA_EPIP_individual_samples_summary.csv'\n")

# Compare sample sizes and show data structure
cat("\n=== DATA STRUCTURE SUMMARY ===\n")
for(i in 1:min(3, length(metal_configs))) {  # Show first 3 metals as examples
  metal_symbol <- names(metal_configs)[i]
  config <- metal_configs[[i]]
  
  if(!config$tissue %in% names(metal_data)) next
  
  plot_data <- prepare_FILA_epip_data(metal_data, config$tissue)
  if(is.null(plot_data)) next
  
  n_sites <- length(unique(paste(plot_data$SITE, plot_data$SAMPLING_DATE)))
  n_pairs <- nrow(plot_data)
  
  cat(sprintf("%s: %d site-date combinations → %d individual sample pairs\n", 
              metal_symbol, n_sites, n_pairs))
}

cat("\nNote: Each site-date combination creates multiple pairs if there are multiple samples per compartment.\n")

# Linear Regression Plots: EPIL vs FILA Metal Concentrations
# 3x4 panel layout for all 12 elements using ALL individual samples (no averaging)

# Load required libraries
library(dplyr)
library(tidyr)

# Read data
metal_data <- read.csv("2_incremental/7_AI_ANALYSIS/metal_analysis_subset.csv")

# Define metal configurations
metal_configs <- list(
  As = list(tissue = "As", name = "Arsenic"),
  Cd = list(tissue = "Cd", name = "Cadmium"),
  Co = list(tissue = "Co", name = "Cobalt"),
  Cr = list(tissue = "Cr", name = "Chromium"),
  Cu = list(tissue = "Cu", name = "Copper"),
  Fe = list(tissue = "Fer.1", name = "Iron"),
  Mn = list(tissue = "Mn", name = "Manganese"),
  Mo = list(tissue = "Mo", name = "Molybdenum"),
  Ni = list(tissue = "Ni", name = "Nickel"),
  Pb = list(tissue = "Pb", name = "Lead"),
  Se = list(tissue = "Se", name = "Selenium"),
  Zn = list(tissue = "Zn", name = "Zinc")
)

# Function to prepare data for EPIL vs FILA comparison using ALL individual samples
prepare_epil_FILA_data <- function(data, tissue_col) {
  
  # Clean and filter data - keep all individual samples
  clean_data <- data %>%
    filter(!is.na(.data[[tissue_col]]) & 
             SAMPLE_DESCRIPTOR %in% c("EPIL", "FILA") &
             !is.na(SITE) & !is.na(SAMPLING_DATE) &
             is.finite(.data[[tissue_col]]))
  
  if(nrow(clean_data) == 0) return(NULL)
  
  # Create unique sample ID for tracking
  clean_data <- clean_data %>%
    mutate(sample_id = paste(SITE, SAMPLING_DATE, SAMPLE_DESCRIPTOR, row_number(), sep = "_")) %>%
    select(SITE, SAMPLING_DATE, SAMPLE_DESCRIPTOR, sample_id, metal_conc = all_of(tissue_col))
  
  # Separate EPIL and FILA samples
  epil_data <- clean_data %>%
    filter(SAMPLE_DESCRIPTOR == "EPIL") %>%
    select(SITE, SAMPLING_DATE, sample_id, epil_conc = metal_conc)
  
  FILA_data <- clean_data %>%
    filter(SAMPLE_DESCRIPTOR == "FILA") %>%
    select(SITE, SAMPLING_DATE, sample_id, FILA_conc = metal_conc)
  
  # For each site-date combination, create all possible EPIL-FILA pairs
  site_dates <- unique(clean_data[, c("SITE", "SAMPLING_DATE")])
  
  all_pairs <- data.frame()
  
  for(i in 1:nrow(site_dates)) {
    site <- site_dates$SITE[i]
    date <- site_dates$SAMPLING_DATE[i]
    
    epil_subset <- epil_data %>% filter(SITE == site & SAMPLING_DATE == date)
    FILA_subset <- FILA_data %>% filter(SITE == site & SAMPLING_DATE == date)
    
    if(nrow(epil_subset) > 0 & nrow(FILA_subset) > 0) {
      # Create all combinations of EPIL and FILA samples for this site-date
      for(j in 1:nrow(epil_subset)) {
        for(k in 1:nrow(FILA_subset)) {
          all_pairs <- rbind(all_pairs, data.frame(
            SITE = site,
            SAMPLING_DATE = date,
            epil_sample_id = epil_subset$sample_id[j],
            FILA_sample_id = FILA_subset$sample_id[k],
            conc_EPIL = epil_subset$epil_conc[j],
            conc_FILA = FILA_subset$FILA_conc[k]
          ))
        }
      }
    }
  }
  
  return(all_pairs)
}

# Function to create individual regression plot with automatic outlier detection
create_regression_plot <- function(data, metal_symbol, metal_name) {
  
  if(is.null(data) || nrow(data) < 10) {  # Need more points for outlier detection
    # Create empty plot for insufficient data
    plot(1, 1, type = "n", xlab = "", ylab = "", main = paste(metal_symbol, "-", metal_name),
         axes = FALSE)
    text(1, 1, "Insufficient\nData", cex = 1.2)
    return()
  }
  
  # Calculate initial regression to identify outliers
  lm_fit_initial <- lm(conc_FILA ~ conc_EPIL, data = data)
  
  # Calculate Cook's distance to identify influential outliers
  cooks_dist <- cooks.distance(lm_fit_initial)
  
  # Define outliers as points with Cook's distance > 4/n (common threshold)
  n_points <- nrow(data)
  cooks_threshold <- 4 / n_points
  outliers <- cooks_dist > cooks_threshold
  
  # Also identify extreme outliers using standardized residuals
  std_residuals <- rstandard(lm_fit_initial)
  extreme_outliers <- abs(std_residuals) > 3  # Points beyond 3 standard deviations
  
  # Combine both criteria
  all_outliers <- outliers | extreme_outliers
  
  # Calculate statistics for both full and cleaned datasets
  # Full dataset
  cor_test_full <- cor.test(data$conc_EPIL, data$conc_FILA, method = "pearson")
  r_full <- round(cor_test_full$estimate, 3)
  r2_full <- round(summary(lm_fit_initial)$r.squared, 3)
  p_full <- cor_test_full$p.value
  
  # Dataset without outliers (if we have enough points left)
  if(sum(!all_outliers) >= 5) {
    data_clean <- data[!all_outliers, ]
    cor_test_clean <- cor.test(data_clean$conc_EPIL, data_clean$conc_FILA, method = "pearson")
    lm_fit_clean <- lm(conc_FILA ~ conc_EPIL, data = data_clean)
    
    r_clean <- round(cor_test_clean$estimate, 3)
    r2_clean <- round(summary(lm_fit_clean)$r.squared, 3)
    p_clean <- cor_test_clean$p.value
    n_clean <- nrow(data_clean)
    n_outliers <- sum(all_outliers)
  } else {
    # Too many outliers removed, use full dataset
    data_clean <- data
    lm_fit_clean <- lm_fit_initial
    cor_test_clean <- cor_test_full
    
    r_clean <- r_full
    r2_clean <- r2_full
    p_clean <- p_full
    n_clean <- n_points
    n_outliers <- 0
    all_outliers <- rep(FALSE, n_points)
  }
  
  # Format p-values
  p_formatted_full <- ifelse(p_full < 0.001, "< 0.001", paste("=", round(p_full, 3)))
  p_formatted_clean <- ifelse(p_clean < 0.001, "< 0.001", paste("=", round(p_clean, 3)))
  
  # Create the plot
  # First plot normal points in blue
  plot(data$conc_EPIL[!all_outliers], data$conc_FILA[!all_outliers],
       xlab = "EPIL Concentration (mg/g)",
       ylab = "FILA Concentration (mg/g)",
       main = paste(metal_symbol, "-", metal_name),
       pch = 16, col = "steelblue", cex = 0.8,
       cex.main = 1.1, cex.lab = 1.0, cex.axis = 0.9,
       xlim = range(data$conc_EPIL), ylim = range(data$conc_FILA))
  
  # Add outliers in gray
  if(any(all_outliers)) {
    points(data$conc_EPIL[all_outliers], data$conc_FILA[all_outliers],
           pch = 16, col = "gray60", cex = 0.8)
  }
  
  # Add regression line using CLEANED data (without outliers)
  abline(lm_fit_clean, col = "red", lwd = 2)
  
  # Add confidence interval for cleaned regression
  x_seq <- seq(min(data$conc_EPIL), max(data$conc_EPIL), length.out = 100)
  pred_interval <- predict(lm_fit_clean, newdata = data.frame(conc_EPIL = x_seq), 
                           interval = "confidence", level = 0.95)
  
  # Add confidence bands
  lines(x_seq, pred_interval[,"lwr"], col = "red", lty = 2, lwd = 1)
  lines(x_seq, pred_interval[,"upr"], col = "red", lty = 2, lwd = 1)
  
  # Create statistics text
  if(n_outliers > 0) {
    stats_text <- paste0(
      "All data (n=", n_points, "):\n",
      "r = ", r_full, ", R² = ", r2_full, "\n",
      "p ", p_formatted_full, "\n\n",
      "Clean data (n=", n_clean, "):\n",
      "r = ", r_clean, ", R² = ", r2_clean, "\n",
      "p ", p_formatted_clean, "\n",
      "Outliers removed: ", n_outliers
    )
  } else {
    stats_text <- paste0(
      "r = ", r_clean, "\n",
      "R² = ", r2_clean, "\n",
      "p ", p_formatted_clean, "\n",
      "n = ", n_clean, "\n",
      "No outliers detected"
    )
  }
  
  # Position text in upper left corner
  x_range <- range(data$conc_EPIL)
  y_range <- range(data$conc_FILA)
  
  text(x_range[1] + 0.05 * diff(x_range), 
       y_range[2] - 0.05 * diff(y_range),
       stats_text, 
       adj = c(0, 1), 
       cex = 0.7, 
       col = "black",
       font = 1)
  
  # Add legend if outliers were detected
  if(any(all_outliers)) {
    legend("bottomright", 
           legend = c("Normal data", "Outliers", "Regression (clean)"),
           col = c("steelblue", "gray60", "red"),
           pch = c(16, 16, NA),
           lty = c(NA, NA, 1),
           lwd = c(NA, NA, 2),
           cex = 0.7,
           bg = "white")
  }
}

# Set up the plotting layout
cat("Generating EPIL vs FILA regression plots using ALL individual samples...\n")

# Create PNG file for high-quality output
png("EPIL_FILA_individual_samples_regression.png", width = 16, height = 12, units = "in", res = 300)

# Set up 3x4 panel layout with margins
par(mfrow = c(3, 4), 
    mar = c(4, 4, 3, 1),  # bottom, left, top, right margins
    oma = c(2, 2, 4, 2))  # outer margins

# Generate plots for all metals
for(i in 1:length(metal_configs)) {
  metal_symbol <- names(metal_configs)[i]
  config <- metal_configs[[i]]
  
  cat("Processing", metal_symbol, "...\n")
  
  # Check if tissue column exists
  if(!config$tissue %in% names(metal_data)) {
    cat("Warning: Column", config$tissue, "not found for", metal_symbol, "\n")
    create_regression_plot(NULL, metal_symbol, config$name)
    next
  }
  
  # Prepare data
  plot_data <- prepare_epil_FILA_data(metal_data, config$tissue)
  
  # Create plot
  create_regression_plot(plot_data, metal_symbol, config$name)
}

# Add overall title
mtext("EPIL vs FILA Metal Concentrations - Individual Sample Pairs", 
      outer = TRUE, cex = 1.5, font = 2, line = 1)

# Add footnote
mtext("Blue points = normal data, Gray points = outliers (excluded from regression). Red line = robust regression with 95% CI.", 
      outer = TRUE, cex = 0.9, line = 0, side = 1)

dev.off()

cat("Plot saved as 'EPIL_FILA_individual_samples_regression.png'\n")

# Also create a summary table of regression statistics
cat("\n=== INDIVIDUAL SAMPLE REGRESSION STATISTICS ===\n")

regression_summary <- data.frame(
  Metal = character(),
  N_pairs = numeric(),
  N_clean = numeric(),
  N_outliers = numeric(),
  R = numeric(),
  R_squared = numeric(),
  P_value = numeric(),
  Slope = numeric(),
  Intercept = numeric(),
  stringsAsFactors = FALSE
)

for(i in 1:length(metal_configs)) {
  metal_symbol <- names(metal_configs)[i]
  config <- metal_configs[[i]]
  
  if(!config$tissue %in% names(metal_data)) next
  
  plot_data <- prepare_epil_FILA_data(metal_data, config$tissue)
  
  if(is.null(plot_data) || nrow(plot_data) < 10) {
    regression_summary <- rbind(regression_summary, data.frame(
      Metal = metal_symbol,
      N_pairs = ifelse(is.null(plot_data), 0, nrow(plot_data)),
      N_clean = NA,
      N_outliers = NA,
      R = NA,
      R_squared = NA,
      P_value = NA,
      Slope = NA,
      Intercept = NA
    ))
    next
  }
  
  # Recalculate outlier detection for summary
  lm_initial <- lm(conc_FILA ~ conc_EPIL, data = plot_data)
  cooks_dist <- cooks.distance(lm_initial)
  std_residuals <- rstandard(lm_initial)
  
  cooks_threshold <- 4 / nrow(plot_data)
  outliers_detected <- (cooks_dist > cooks_threshold) | (abs(std_residuals) > 3)
  
  # Calculate final statistics using cleaned data
  if(sum(!outliers_detected) >= 5) {
    clean_data <- plot_data[!outliers_detected, ]
    cor_test <- cor.test(clean_data$conc_EPIL, clean_data$conc_FILA, method = "pearson")
    lm_fit <- lm(conc_FILA ~ conc_EPIL, data = clean_data)
  } else {
    # Use full dataset if too many outliers
    cor_test <- cor.test(plot_data$conc_EPIL, plot_data$conc_FILA, method = "pearson")
    lm_fit <- lm(conc_FILA ~ conc_EPIL, data = plot_data)
    outliers_detected <- rep(FALSE, nrow(plot_data))
  }
  
  regression_summary <- rbind(regression_summary, data.frame(
    Metal = metal_symbol,
    N_pairs = nrow(plot_data),
    R = round(cor_test$estimate, 3),
    R_squared = round(summary(lm_fit)$r.squared, 3),
    P_value = round(cor_test$p.value, 4),
    Slope = round(coef(lm_fit)[2], 4),
    Intercept = round(coef(lm_fit)[1], 4)
  ))
}

print(regression_summary, row.names = FALSE)

# Save summary table
write.csv(regression_summary, "EPIL_FILA_individual_samples_summary.csv", row.names = FALSE)

cat("\nRegression summary saved as 'EPIL_FILA_individual_samples_summary.csv'\n")

# Compare sample sizes and show data structure
cat("\n=== DATA STRUCTURE SUMMARY ===\n")
for(i in 1:min(3, length(metal_configs))) {  # Show first 3 metals as examples
  metal_symbol <- names(metal_configs)[i]
  config <- metal_configs[[i]]
  
  if(!config$tissue %in% names(metal_data)) next
  
  plot_data <- prepare_epil_FILA_data(metal_data, config$tissue)
  if(is.null(plot_data)) next
  
  n_sites <- length(unique(paste(plot_data$SITE, plot_data$SAMPLING_DATE)))
  n_pairs <- nrow(plot_data)
  
  cat(sprintf("%s: %d site-date combinations → %d individual sample pairs\n", 
              metal_symbol, n_sites, n_pairs))
}

cat("\nNote: Each site-date combination creates multiple pairs if there are multiple samples per compartment.\n")