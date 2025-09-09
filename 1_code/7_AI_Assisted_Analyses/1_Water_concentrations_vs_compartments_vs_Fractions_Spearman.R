# Fixed Multi-Metal Analysis with Nested Water Fractions and Fraction Percentages
# Load required libraries
library(dplyr)

# Read subset data
metal_data <- read.csv("2_incremental/7_AI_ANALYSIS/metal_analysis_subset_clean.csv")

# Define metal configurations
metal_configs <- list(
  As = list(tissue = "As", water = c("UF_As_Oshift", "F_As_Oshift", "W_As_Oshift")),
  Cd = list(tissue = "Cd", water = c("UF_Cd114", "F_Cd114", "W_Cd114")),
  Co = list(tissue = "Co", water = c("UF_Co", "F_Co", "W_Co")),
  Cr = list(tissue = "Cr", water = c("UF_Cr_NH3", "F_Cr_NH3", "W_Cr_NH3")),
  Cu = list(tissue = "Cu", water = c("UF_Cu_NH3", "F_Cu_NH3", "W_Cu_NH3")),
  Fe = list(tissue = "Fer.1", water = c("UF_Fe_NH3", "F_Fe_NH3", "W_Fe_NH3")),
  Mn = list(tissue = "Mn", water = c("UF_Mn", "F_Mn", "W_Mn")),
  Mo = list(tissue = "Mo", water = c("UF_Mo_Oshift", "F_Mo_Oshift", "W_Mo_Oshift")),
  Ni = list(tissue = "Ni", water = c("UF_Ni", "F_Ni", "W_Ni")),
  Pb = list(tissue = "Pb", water = c("UF_Pb_NH3", "F_Pb_NH3", "W_Pb_NH3")),
  Se = list(tissue = "Se", water = c("UF_Se82_Oshift", "F_Se82_Oshift", "W_Se82_Oshift")),
  Zn = list(tissue = "Zn", water = c("UF_Zn", "F_Zn", "W_Zn"))
)

main_compartments <- c("EPIL", "EPIP", "FILA")

# Function to analyze single metal
analyze_metal_with_percentages <- function(data, metal_name, config) {
  
  tissue_col <- config$tissue
  water_cols <- config$water[config$water %in% names(data)]
  
  # Clean data
  clean_data <- data %>%
    filter(!is.na(.data[[tissue_col]]) & 
             SAMPLE_DESCRIPTOR %in% main_compartments &
             !is.na(SITE) & !is.na(SAMPLING_DATE))
  
  if(nrow(clean_data) == 0) return(NULL)
  
  # Average by site/date/compartment
  avg_cols <- c("SITE", "SAMPLING_DATE", "SAMPLE_DESCRIPTOR", tissue_col, water_cols)
  existing_cols <- avg_cols[avg_cols %in% names(clean_data)]
  
  averaged_data <- clean_data %>%
    select(all_of(existing_cols)) %>%
    group_by(SITE, SAMPLING_DATE, SAMPLE_DESCRIPTOR) %>%
    summarise(across(all_of(c(tissue_col, water_cols)), 
                     ~mean(.x, na.rm = TRUE)), 
              n_samples = n(),
              .groups = 'drop') %>%
    filter(is.finite(.data[[tissue_col]]))
  
  # Calculate derived fractions and percentages
  uf_col <- water_cols[grepl("^UF_", water_cols)][1]
  f_col <- water_cols[grepl("^F_", water_cols)][1] 
  w_col <- water_cols[grepl("^W_", water_cols)][1]
  
  if(!is.na(uf_col) && !is.na(f_col) && !is.na(w_col)) {
    averaged_data <- averaged_data %>%
      mutate(
        colloidal_fraction = .data[[f_col]] - .data[[uf_col]],
        particulate_fraction = .data[[w_col]] - .data[[f_col]],
        pct_dissolved = (.data[[uf_col]] / .data[[w_col]]) * 100,
        pct_colloidal = ((.data[[f_col]] - .data[[uf_col]]) / .data[[w_col]]) * 100,
        pct_particulate = ((.data[[w_col]] - .data[[f_col]]) / .data[[w_col]]) * 100
      )
  }
  
  # Calculate tissue statistics
  tissue_stats <- averaged_data %>%
    group_by(SAMPLE_DESCRIPTOR) %>%
    summarise(n_site_date = n(),
              mean_metal = mean(.data[[tissue_col]], na.rm = TRUE),
              .groups = 'drop')
  
  # Calculate correlations
  correlation_results <- data.frame()
  
  for(comp in main_compartments) {
    comp_data <- averaged_data %>% filter(SAMPLE_DESCRIPTOR == comp)
    if(nrow(comp_data) < 5) next
    
    for(water_col in water_cols) {
      pairs <- comp_data %>%
        filter(!is.na(.data[[tissue_col]]) & 
                 !is.na(.data[[water_col]]) & 
                 is.finite(.data[[water_col]])) %>%
        select(all_of(c(tissue_col, water_col)))
      
      if(nrow(pairs) >= 5) {
        cor_result <- cor.test(pairs[[tissue_col]], pairs[[water_col]], method = "spearman")
        
        fraction_type <- case_when(
          grepl("^UF_", water_col) ~ "UF",
          grepl("^F_", water_col) ~ "F", 
          grepl("^W_", water_col) ~ "W",
          TRUE ~ "Unknown"
        )
        
        correlation_results <- rbind(correlation_results, 
                                     data.frame(
                                       Compartment = comp,
                                       Water_Fraction = water_col,
                                       Fraction_Type = fraction_type,
                                       Correlation = cor_result$estimate,
                                       N = nrow(pairs)
                                     ))
      }
    }
  }
  
  # Fraction hierarchy analysis
  fraction_hierarchy <- if(nrow(correlation_results) > 0) {
    correlation_results %>%
      filter(Fraction_Type %in% c("UF", "F", "W")) %>%
      group_by(Fraction_Type) %>%
      summarise(mean_abs_r = mean(abs(Correlation), na.rm = TRUE), .groups = 'drop') %>%
      arrange(desc(mean_abs_r))
  } else {
    data.frame(Fraction_Type = c("UF", "F", "W"), mean_abs_r = c(NA, NA, NA))
  }
  
  # Calculate fraction percentages
  fraction_percentages <- if("pct_dissolved" %in% names(averaged_data)) {
    averaged_data %>%
      filter(is.finite(pct_dissolved) & is.finite(pct_colloidal) & is.finite(pct_particulate) &
               pct_dissolved >= 0 & pct_colloidal >= 0 & pct_particulate >= 0) %>%
      summarise(
        mean_pct_dissolved = mean(pct_dissolved, na.rm = TRUE),
        mean_pct_colloidal = mean(pct_colloidal, na.rm = TRUE), 
        mean_pct_particulate = mean(pct_particulate, na.rm = TRUE),
        .groups = 'drop'
      )
  } else {
    data.frame(mean_pct_dissolved = NA, mean_pct_colloidal = NA, mean_pct_particulate = NA)
  }
  
  # Determine bioavailability pattern
  if(nrow(fraction_hierarchy) >= 3 && any(!is.na(fraction_hierarchy$mean_abs_r))) {
    hierarchy_string <- paste(fraction_hierarchy$Fraction_Type[!is.na(fraction_hierarchy$mean_abs_r)], collapse = " > ")
    
    # Check for truly dissolved pattern
    dissolved_pct <- if(nrow(fraction_percentages) > 0 && !is.na(fraction_percentages$mean_pct_dissolved)) {
      fraction_percentages$mean_pct_dissolved
    } else { NA }
    
    # Pattern classification
    if(!is.na(dissolved_pct) && dissolved_pct > 90) {
      cor_range <- max(fraction_hierarchy$mean_abs_r, na.rm = TRUE) - min(fraction_hierarchy$mean_abs_r, na.rm = TRUE)
      if(cor_range < 0.15) {
        bioavail_pattern <- "Truly Dissolved (>90%)"
      } else {
        bioavail_pattern <- "Dissolved + Colloidal Uptake"
      }
    } else if(hierarchy_string == "F > UF > W" || hierarchy_string == "F > W > UF") {
      bioavail_pattern <- "Dissolved + Colloidal Uptake"
    } else if(hierarchy_string == "UF > F > W" || hierarchy_string == "UF > W > F") {
      bioavail_pattern <- "Dissolved-Only Uptake"
    } else if(hierarchy_string == "W > F > UF" || hierarchy_string == "W > UF > F") {
      bioavail_pattern <- "Particulate-Enhanced Uptake"
    } else {
      bioavail_pattern <- "Mixed/Complex Uptake"
    }
  } else {
    hierarchy_string <- "No data"
    bioavail_pattern <- "No Data"
  }
  
  # Calculate strongest correlations by compartment
  strongest_cors <- if(nrow(correlation_results) > 0) {
    correlation_results %>%
      filter(Fraction_Type %in% c("UF", "F", "W")) %>%
      group_by(Compartment) %>%
      summarise(strongest_r = max(abs(Correlation), na.rm = TRUE), .groups = 'drop')
  } else {
    data.frame(Compartment = main_compartments, strongest_r = NA)
  }
  
  # Calculate BCF
  uf_col <- water_cols[grepl("^UF_", water_cols)][1]
  bcf_data <- if(!is.na(uf_col)) {
    averaged_data %>%
      filter(!is.na(.data[[tissue_col]]) & 
               !is.na(.data[[uf_col]]) & 
               .data[[uf_col]] > 0) %>%
      mutate(bcf = .data[[tissue_col]] / .data[[uf_col]]) %>%
      group_by(SAMPLE_DESCRIPTOR) %>%
      summarise(mean_bcf = mean(bcf, na.rm = TRUE), .groups = 'drop')
  } else {
    data.frame(SAMPLE_DESCRIPTOR = main_compartments, mean_bcf = NA)
  }
  
  # Return results
  list(
    total_combinations = nrow(averaged_data),
    tissue_stats = tissue_stats,
    strongest_cors = strongest_cors,
    bcf_data = bcf_data,
    fraction_hierarchy = fraction_hierarchy,
    hierarchy_string = hierarchy_string,
    bioavail_pattern = bioavail_pattern,
    n_strong_correlations = sum(abs(correlation_results$Correlation[correlation_results$Fraction_Type %in% c("UF", "F", "W")]) > 0.3, na.rm = TRUE),
    tissue_hierarchy = paste(tissue_stats$SAMPLE_DESCRIPTOR[order(-tissue_stats$mean_metal)], collapse = " > "),
    fraction_percentages = fraction_percentages
  )
}

# Run analysis
results <- list()

cat("Running enhanced multi-metal analysis...\n")

for(metal_name in names(metal_configs)) {
  config <- metal_configs[[metal_name]]
  
  if(!config$tissue %in% names(metal_data)) {
    cat("Skipping", metal_name, "- tissue column not found\n")
    next
  }
  
  cat("Analyzing", metal_name, "...\n")
  result <- analyze_metal_with_percentages(metal_data, metal_name, config)
  
  if(!is.null(result)) {
    results[[metal_name]] <- result
  }
}

# Create enhanced summary
create_enhanced_summary <- function(results) {
  
  summary_df <- data.frame(
    Metric = c(
      "Site/date combinations",
      "EPIL combinations", 
      "EPIP combinations",
      "FILA combinations",
      "EPIL mean concentration",
      "EPIP mean concentration", 
      "FILA mean concentration",
      "EPIL strongest correlation",
      "EPIP strongest correlation",
      "FILA strongest correlation",
      "EPIL mean BCF",
      "EPIP mean BCF",
      "FILA mean BCF",
      "Correlations > |0.3|",
      "Tissue hierarchy",
      "Fraction hierarchy",
      "Bioavailability pattern",
      "UF mean correlation",
      "F mean correlation", 
      "W mean correlation",
      "% Truly dissolved (UF)",
      "% Colloidal (F-UF)", 
      "% Particulate (W-F)"
    ),
    stringsAsFactors = FALSE
  )
  
  for(metal_name in names(results)) {
    result <- results[[metal_name]]
    
    # Extract data
    epil_stats <- result$tissue_stats[result$tissue_stats$SAMPLE_DESCRIPTOR == "EPIL", ]
    epip_stats <- result$tissue_stats[result$tissue_stats$SAMPLE_DESCRIPTOR == "EPIP", ]
    fila_stats <- result$tissue_stats[result$tissue_stats$SAMPLE_DESCRIPTOR == "FILA", ]
    
    epil_cor <- result$strongest_cors$strongest_r[result$strongest_cors$Compartment == "EPIL"]
    epip_cor <- result$strongest_cors$strongest_r[result$strongest_cors$Compartment == "EPIP"]
    fila_cor <- result$strongest_cors$strongest_r[result$strongest_cors$Compartment == "FILA"]
    
    epil_bcf <- result$bcf_data$mean_bcf[result$bcf_data$SAMPLE_DESCRIPTOR == "EPIL"]
    epip_bcf <- result$bcf_data$mean_bcf[result$bcf_data$SAMPLE_DESCRIPTOR == "EPIP"]
    fila_bcf <- result$bcf_data$mean_bcf[result$bcf_data$SAMPLE_DESCRIPTOR == "FILA"]
    
    # Extract fraction data
    frac_hier <- result$fraction_hierarchy
    uf_cor <- frac_hier$mean_abs_r[frac_hier$Fraction_Type == "UF"]
    f_cor <- frac_hier$mean_abs_r[frac_hier$Fraction_Type == "F"]
    w_cor <- frac_hier$mean_abs_r[frac_hier$Fraction_Type == "W"]
    
    if(length(uf_cor) == 0) uf_cor <- NA
    if(length(f_cor) == 0) f_cor <- NA
    if(length(w_cor) == 0) w_cor <- NA
    
    # Extract percentages
    frac_pcts <- result$fraction_percentages
    pct_dissolved <- if(nrow(frac_pcts) > 0 && !is.na(frac_pcts$mean_pct_dissolved)) {
      round(frac_pcts$mean_pct_dissolved, 1)
    } else { NA }
    pct_colloidal <- if(nrow(frac_pcts) > 0 && !is.na(frac_pcts$mean_pct_colloidal)) {
      round(frac_pcts$mean_pct_colloidal, 1)
    } else { NA }
    pct_particulate <- if(nrow(frac_pcts) > 0 && !is.na(frac_pcts$mean_pct_particulate)) {
      round(frac_pcts$mean_pct_particulate, 1)
    } else { NA }
    
    summary_df[[metal_name]] <- c(
      result$total_combinations,
      ifelse(nrow(epil_stats) > 0, epil_stats$n_site_date, 0),
      ifelse(nrow(epip_stats) > 0, epip_stats$n_site_date, 0),
      ifelse(nrow(fila_stats) > 0, fila_stats$n_site_date, 0),
      ifelse(nrow(epil_stats) > 0, round(epil_stats$mean_metal, 5), NA),
      ifelse(nrow(epip_stats) > 0, round(epip_stats$mean_metal, 5), NA),
      ifelse(nrow(fila_stats) > 0, round(fila_stats$mean_metal, 5), NA),
      ifelse(length(epil_cor) > 0 && !is.na(epil_cor), round(epil_cor, 3), NA),
      ifelse(length(epip_cor) > 0 && !is.na(epip_cor), round(epip_cor, 3), NA),
      ifelse(length(fila_cor) > 0 && !is.na(fila_cor), round(fila_cor, 3), NA),
      ifelse(length(epil_bcf) > 0 && !is.na(epil_bcf), round(epil_bcf, 3), NA),
      ifelse(length(epip_bcf) > 0 && !is.na(epip_bcf), round(epip_bcf, 3), NA),
      ifelse(length(fila_bcf) > 0 && !is.na(fila_bcf), round(fila_bcf, 3), NA),
      result$n_strong_correlations,
      result$tissue_hierarchy,
      result$hierarchy_string,
      result$bioavail_pattern,
      ifelse(!is.na(uf_cor), round(uf_cor, 3), NA),
      ifelse(!is.na(f_cor), round(f_cor, 3), NA),
      ifelse(!is.na(w_cor), round(w_cor, 3), NA),
      pct_dissolved,
      pct_colloidal,
      pct_particulate
    )
  }
  
  summary_df
}

# Generate results
enhanced_summary <- create_enhanced_summary(results)

# Display results
cat("\n=== ENHANCED ANALYSIS WITH FRACTION PERCENTAGES ===\n")
print(enhanced_summary, row.names = FALSE)

# Create bioavailability patterns summary
bioavail_patterns <- data.frame(
  Metal = names(results),
  Pattern = sapply(results, function(x) x$bioavail_pattern),
  Hierarchy = sapply(results, function(x) x$hierarchy_string),
  Pct_Dissolved = sapply(results, function(x) {
    if(nrow(x$fraction_percentages) > 0 && !is.na(x$fraction_percentages$mean_pct_dissolved)) {
      round(x$fraction_percentages$mean_pct_dissolved, 1)
    } else { NA }
  }),
  Pct_Colloidal = sapply(results, function(x) {
    if(nrow(x$fraction_percentages) > 0 && !is.na(x$fraction_percentages$mean_pct_colloidal)) {
      round(x$fraction_percentages$mean_pct_colloidal, 1)
    } else { NA }
  }),
  Pct_Particulate = sapply(results, function(x) {
    if(nrow(x$fraction_percentages) > 0 && !is.na(x$fraction_percentages$mean_pct_particulate)) {
      round(x$fraction_percentages$mean_pct_particulate, 1)
    } else { NA }
  }),
  stringsAsFactors = FALSE
)

cat("\n=== BIOAVAILABILITY PATTERNS WITH FRACTION PERCENTAGES ===\n")
print(bioavail_patterns, row.names = FALSE)

# Save results
write.csv(enhanced_summary, "2_incremental/7_AI_ANALYSIS/enhanced_summary_with_percentages.csv", row.names = FALSE)
write.csv(bioavail_patterns, "2_incremental/7_AI_ANALYSIS/bioavailability_patterns_with_percentages.csv", row.names = FALSE)

cat("\nAnalysis complete! Files saved:\n")
cat("- enhanced_summary_with_percentages.csv\n")
cat("- bioavailability_patterns_with_percentages.csv\n")