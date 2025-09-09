# Comprehensive Compartment Regression Analysis
# EPIL vs EPIP, FILA vs EPIP, and EPIL vs FILA with outlier detection

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

# Define compartment pairs for analysis
compartment_pairs <- list(
  list(x_comp = "EPIL", y_comp = "EPIP", name = "EPIL_vs_EPIP"),
  list(x_comp = "FILA", y_comp = "EPIP", name = "FILA_vs_EPIP"),
  list(x_comp = "EPIL", y_comp = "FILA", name = "EPIL_vs_FILA")
)

# Function to prepare data for any compartment pair comparison
prepare_compartment_data <- function(data, tissue_col, x_comp, y_comp) {
  
  # Clean and filter data
  clean_data <- data %>%
    filter(!is.na(.data[[tissue_col]]) & 
             SAMPLE_DESCRIPTOR %in% c(x_comp, y_comp) &
             !is.na(SITE) & !is.na(SAMPLING_DATE) &
             is.finite(.data[[tissue_col]]))
  
  if(nrow(clean_data) == 0) return(NULL)
  
  # Create unique sample ID for tracking
  clean_data <- clean_data %>%
    mutate(sample_id = paste(SITE, SAMPLING_DATE, SAMPLE_DESCRIPTOR, row_number(), sep = "_")) %>%
    select(SITE, SAMPLING_DATE, SAMPLE_DESCRIPTOR, sample_id, metal_conc = all_of(tissue_col))
  
  # Separate X and Y compartment samples
  x_data <- clean_data %>%
    filter(SAMPLE_DESCRIPTOR == x_comp) %>%
    select(SITE, SAMPLING_DATE, sample_id, x_conc = metal_conc)
  
  y_data <- clean_data %>%
    filter(SAMPLE_DESCRIPTOR == y_comp) %>%
    select(SITE, SAMPLING_DATE, sample_id, y_conc = metal_conc)
  
  # Create all possible pairs for each site-date combination
  site_dates <- unique(clean_data[, c("SITE", "SAMPLING_DATE")])
  all_pairs <- data.frame()
  
  for(i in 1:nrow(site_dates)) {
    site <- site_dates$SITE[i]
    date <- site_dates$SAMPLING_DATE[i]
    
    x_subset <- x_data %>% filter(SITE == site & SAMPLING_DATE == date)
    y_subset <- y_data %>% filter(SITE == site & SAMPLING_DATE == date)
    
    if(nrow(x_subset) > 0 & nrow(y_subset) > 0) {
      for(j in 1:nrow(x_subset)) {
        for(k in 1:nrow(y_subset)) {
          all_pairs <- rbind(all_pairs, data.frame(
            SITE = site,
            SAMPLING_DATE = date,
            x_sample_id = x_subset$sample_id[j],
            y_sample_id = y_subset$sample_id[k],
            x_conc = x_subset$x_conc[j],
            y_conc = y_subset$y_conc[k]
          ))
        }
      }
    }
  }
  
  return(all_pairs)
}

# Function to analyze outliers and regression
analyze_regression_with_outliers <- function(data, metal_symbol, x_comp, y_comp) {
  
  if(is.null(data) || nrow(data) < 10) {
    return(list(
      metal = metal_symbol,
      comparison = paste(x_comp, "vs", y_comp),
      n_total = ifelse(is.null(data), 0, nrow(data)),
      n_clean = NA,
      n_outliers = NA,
      r_full = NA,
      r2_full = NA,
      p_full = NA,
      r_clean = NA,
      r2_clean = NA,
      p_clean = NA,
      slope_clean = NA,
      intercept_clean = NA,
      outlier_samples = NA
    ))
  }
  
  # Initial regression
  lm_fit_initial <- lm(y_conc ~ x_conc, data = data)
  
  # Outlier detection
  cooks_dist <- cooks.distance(lm_fit_initial)
  std_residuals <- rstandard(lm_fit_initial)
  
  cooks_threshold <- 4 / nrow(data)
  all_outliers <- (cooks_dist > cooks_threshold) | (abs(std_residuals) > 3)
  
  # Full dataset statistics
  cor_test_full <- cor.test(data$x_conc, data$y_conc, method = "pearson")
  r_full <- cor_test_full$estimate
  r2_full <- summary(lm_fit_initial)$r.squared
  p_full <- cor_test_full$p.value
  
  # Clean dataset statistics
  if(sum(!all_outliers) >= 5) {
    data_clean <- data[!all_outliers, ]
    cor_test_clean <- cor.test(data_clean$x_conc, data_clean$y_conc, method = "pearson")
    lm_fit_clean <- lm(y_conc ~ x_conc, data = data_clean)
    
    r_clean <- cor_test_clean$estimate
    r2_clean <- summary(lm_fit_clean)$r.squared
    p_clean <- cor_test_clean$p.value
    slope_clean <- coef(lm_fit_clean)[2]
    intercept_clean <- coef(lm_fit_clean)[1]
    n_clean <- nrow(data_clean)
  } else {
    # Use full dataset if too many outliers
    r_clean <- r_full
    r2_clean <- r2_full
    p_clean <- p_full
    slope_clean <- coef(lm_fit_initial)[2]
    intercept_clean <- coef(lm_fit_initial)[1]
    n_clean <- nrow(data)
    all_outliers <- rep(FALSE, nrow(data))
  }
  
  # Identify outlier samples
  outlier_info <- if(any(all_outliers)) {
    outlier_data <- data[all_outliers, ]
    paste0(outlier_data$SITE, "_", outlier_data$SAMPLING_DATE, " (", 
           round(outlier_data$x_conc, 4), ",", round(outlier_data$y_conc, 4), ")",
           collapse = "; ")
  } else {
    "None"
  }
  
  return(list(
    metal = metal_symbol,
    comparison = paste(x_comp, "vs", y_comp),
    n_total = nrow(data),
    n_clean = n_clean,
    n_outliers = sum(all_outliers),
    r_full = r_full,
    r2_full = r2_full,
    p_full = p_full,
    r_clean = r_clean,
    r2_clean = r2_clean,
    p_clean = p_clean,
    slope_clean = slope_clean,
    intercept_clean = intercept_clean,
    outlier_samples = outlier_info
  ))
}

# Function to create regression plot
create_regression_plot <- function(data, metal_symbol, metal_name, x_comp, y_comp) {
  
  if(is.null(data) || nrow(data) < 10) {
    plot(1, 1, type = "n", xlab = "", ylab = "", 
         main = paste(metal_symbol, "-", metal_name),
         axes = FALSE)
    text(1, 1, "Insufficient\nData", cex = 1.2)
    return()
  }
  
  # Outlier detection
  lm_fit_initial <- lm(y_conc ~ x_conc, data = data)
  cooks_dist <- cooks.distance(lm_fit_initial)
  std_residuals <- rstandard(lm_fit_initial)
  
  cooks_threshold <- 4 / nrow(data)
  all_outliers <- (cooks_dist > cooks_threshold) | (abs(std_residuals) > 3)
  
  # Statistics
  cor_test_full <- cor.test(data$x_conc, data$y_conc, method = "pearson")
  
  if(sum(!all_outliers) >= 5) {
    data_clean <- data[!all_outliers, ]
    cor_test_clean <- cor.test(data_clean$x_conc, data_clean$y_conc, method = "pearson")
    lm_fit_clean <- lm(y_conc ~ x_conc, data = data_clean)
  } else {
    data_clean <- data
    lm_fit_clean <- lm_fit_initial
    cor_test_clean <- cor_test_full
    all_outliers <- rep(FALSE, nrow(data))
  }
  
  # Create plot
  plot(data$x_conc[!all_outliers], data$y_conc[!all_outliers],
       xlab = paste(x_comp, "Concentration (mg/g)"),
       ylab = paste(y_comp, "Concentration (mg/g)"),
       main = paste(metal_symbol, "-", metal_name),
       pch = 16, col = "steelblue", cex = 0.8,
       cex.main = 1.1, cex.lab = 1.0, cex.axis = 0.9,
       xlim = range(data$x_conc), ylim = range(data$y_conc))
  
  # Add outliers in gray
  if(any(all_outliers)) {
    points(data$x_conc[all_outliers], data$y_conc[all_outliers],
           pch = 16, col = "gray60", cex = 0.8)
  }
  
  # Add regression line and confidence intervals
  abline(lm_fit_clean, col = "red", lwd = 2)
  
  x_seq <- seq(min(data$x_conc), max(data$x_conc), length.out = 100)
  pred_interval <- predict(lm_fit_clean, newdata = data.frame(x_conc = x_seq), 
                           interval = "confidence", level = 0.95)
  lines(x_seq, pred_interval[,"lwr"], col = "red", lty = 2, lwd = 1)
  lines(x_seq, pred_interval[,"upr"], col = "red", lty = 2, lwd = 1)
  
  # Add statistics text
  r_clean <- round(cor_test_clean$estimate, 3)
  r2_clean <- round(summary(lm_fit_clean)$r.squared, 3)
  p_clean <- ifelse(cor_test_clean$p.value < 0.001, "< 0.001", 
                    paste("=", round(cor_test_clean$p.value, 3)))
  
  stats_text <- paste0("r = ", r_clean, "\nR² = ", r2_clean, "\np ", p_clean, 
                       "\nn = ", nrow(data_clean))
  if(any(all_outliers)) {
    stats_text <- paste0(stats_text, "\nOutliers: ", sum(all_outliers))
  }
  
  x_range <- range(data$x_conc)
  y_range <- range(data$y_conc)
  text(x_range[1] + 0.05 * diff(x_range), 
       y_range[2] - 0.05 * diff(y_range),
       stats_text, adj = c(0, 1), cex = 0.7, col = "black", font = 1)
}

# Main analysis loop
cat("Running comprehensive compartment regression analysis...\n")

all_results <- list()
result_counter <- 1

for(pair_info in compartment_pairs) {
  x_comp <- pair_info$x_comp
  y_comp <- pair_info$y_comp
  pair_name <- pair_info$name
  
  cat("\n=== ANALYZING", pair_name, "===\n")
  
  # Create plots
  png(paste0(pair_name, "_regression_plots.png"), width = 16, height = 12, units = "in", res = 300)
  par(mfrow = c(3, 4), mar = c(4, 4, 3, 1), oma = c(2, 2, 4, 2))
  
  for(i in 1:length(metal_configs)) {
    metal_symbol <- names(metal_configs)[i]
    config <- metal_configs[[i]]
    
    if(!config$tissue %in% names(metal_data)) {
      create_regression_plot(NULL, metal_symbol, config$name, x_comp, y_comp)
      next
    }
    
    # Prepare data
    plot_data <- prepare_compartment_data(metal_data, config$tissue, x_comp, y_comp)
    
    # Analyze regression
    result <- analyze_regression_with_outliers(plot_data, metal_symbol, x_comp, y_comp)
    all_results[[result_counter]] <- result
    result_counter <- result_counter + 1
    
    # Create plot
    create_regression_plot(plot_data, metal_symbol, config$name, x_comp, y_comp)
  }
  
  # Add title and save
  mtext(paste(x_comp, "vs", y_comp, "Metal Concentrations - Individual Sample Pairs"), 
        outer = TRUE, cex = 1.5, font = 2, line = 1)
  mtext("Blue = normal data, Gray = outliers. Red line = robust regression with 95% CI.", 
        outer = TRUE, cex = 0.9, line = 0, side = 1)
  
  dev.off()
  cat("Plot saved as", paste0(pair_name, "_regression_plots.png"), "\n")
}

# Create comprehensive scorecard
cat("\n=== CREATING COMPREHENSIVE SCORECARD ===\n")

scorecard <- data.frame(
  Metal = character(),
  Comparison = character(),
  N_Total = numeric(),
  N_Clean = numeric(),
  N_Outliers = numeric(),
  R_Full = numeric(),
  R2_Full = numeric(),
  P_Full = numeric(),
  R_Clean = numeric(),
  R2_Clean = numeric(),
  P_Clean = numeric(),
  Slope_Clean = numeric(),
  Intercept_Clean = numeric(),
  Outlier_Samples = character(),
  stringsAsFactors = FALSE
)

for(result in all_results) {
  scorecard <- rbind(scorecard, data.frame(
    Metal = result$metal,
    Comparison = result$comparison,
    N_Total = result$n_total,
    N_Clean = result$n_clean,
    N_Outliers = result$n_outliers,
    R_Full = round(result$r_full, 3),
    R2_Full = round(result$r2_full, 3),
    P_Full = round(result$p_full, 4),
    R_Clean = round(result$r_clean, 3),
    R2_Clean = round(result$r2_clean, 3),
    P_Clean = round(result$p_clean, 4),
    Slope_Clean = round(result$slope_clean, 4),
    Intercept_Clean = round(result$intercept_clean, 4),
    Outlier_Samples = result$outlier_samples,
    stringsAsFactors = FALSE
  ))
}

# Display and save scorecard
print(scorecard, row.names = FALSE)
write.csv(scorecard, "comprehensive_compartment_regression_scorecard.csv", row.names = FALSE)

# Create summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")

summary_stats <- scorecard %>%
  group_by(Comparison) %>%
  summarise(
    N_Metals = n(),
    Avg_N_Total = round(mean(N_Total, na.rm = TRUE), 0),
    Avg_R_Clean = round(mean(abs(R_Clean), na.rm = TRUE), 3),
    Avg_R2_Clean = round(mean(R2_Clean, na.rm = TRUE), 3),
    Strong_Correlations = sum(abs(R_Clean) > 0.5 & P_Clean < 0.05, na.rm = TRUE),
    Total_Outliers = sum(N_Outliers, na.rm = TRUE),
    Avg_Outliers_Per_Metal = round(mean(N_Outliers, na.rm = TRUE), 1),
    .groups = 'drop'
  )

print(summary_stats, row.names = FALSE)

# Identify most problematic samples (frequent outliers)
cat("\n=== FREQUENT OUTLIER SAMPLES ===\n")

outlier_samples <- scorecard %>%
  filter(Outlier_Samples != "None" & !is.na(Outlier_Samples)) %>%
  select(Metal, Comparison, Outlier_Samples)

if(nrow(outlier_samples) > 0) {
  # Parse outlier strings to identify frequent problematic samples
  all_outlier_sites <- c()
  for(i in 1:nrow(outlier_samples)) {
    outliers <- strsplit(outlier_samples$Outlier_Samples[i], "; ")[[1]]
    sites <- sapply(strsplit(outliers, " \\("), function(x) x[1])
    all_outlier_sites <- c(all_outlier_sites, sites)
  }
  
  outlier_frequency <- table(all_outlier_sites)
  frequent_outliers <- sort(outlier_frequency, decreasing = TRUE)
  
  cat("Sample sites flagged as outliers (frequency):\n")
  print(frequent_outliers)
} else {
  cat("No outliers detected across all regressions.\n")
}

# Save all results
write.csv(summary_stats, "compartment_regression_summary.csv", row.names = FALSE)
write.csv(outlier_samples, "outlier_samples_detailed.csv", row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Files saved:\n")
cat("- comprehensive_compartment_regression_scorecard.csv\n")
cat("- compartment_regression_summary.csv\n") 
cat("- outlier_samples_detailed.csv\n")
cat("- EPIL_vs_EPIP_regression_plots.png\n")
cat("- FILA_vs_EPIP_regression_plots.png\n")
cat("- EPIL_vs_FILA_regression_plots.png\n")