# Hypothesis 2: Compartment Interactions Analysis
# Testing whether compartment abundances affect metal concentrations in other compartments
# Based on corrected nested fraction understanding from Hypothesis 1

# Load required libraries
library(dplyr)
library(tidyr)

# Read subset data
metal_data <- read.csv("2_incremental/7_AI_ANALYSIS/metal_analysis_subset.csv")

# Define metal configurations with detailed speciation understanding
metal_configs <- list(
  As = list(
    tissue = "As", 
    water = c("UF_As_Oshift", "F_As_Oshift", "W_As_Oshift"),
    dissolved_forms = "AsO4³⁻, AsO3³⁻ (truly dissolved)",
    colloidal_forms = "As-Fe/Al oxide nanoparticles, As-organic colloids", 
    particulate_forms = "As-Fe/Al oxides >0.47μm, As-sulfide particles",
    bioavail_interpretation = "F>UF>W suggests dissolved + colloidal uptake"
  ),
  Cd = list(
    tissue = "Cd", 
    water = c("UF_Cd114", "F_Cd114", "W_Cd114"),
    dissolved_forms = "Cd²⁺, CdCl⁺, CdOH⁺ (truly dissolved)",
    colloidal_forms = "Cd-humic/fulvic acid complexes, Cd-Fe/Mn oxide colloids",
    particulate_forms = "CdS particles, Cd-Fe/Mn oxide particles, Cd-clay particles",
    bioavail_interpretation = "F>UF>W suggests dissolved + colloidal uptake"
  ),
  Co = list(
    tissue = "Co", 
    water = c("UF_Co", "F_Co", "W_Co"),
    dissolved_forms = "Co²⁺, CoOH⁺ (truly dissolved)",
    colloidal_forms = "Co-Mn oxide colloids, Co-organic colloids",
    particulate_forms = "Co-Mn/Fe oxide particles, Co-sulfide particles",
    bioavail_interpretation = "Mixed pattern suggests multiple uptake pathways"
  ),
  Cr = list(
    tissue = "Cr", 
    water = c("UF_Cr_NH3", "F_Cr_NH3", "W_Cr_NH3"),
    dissolved_forms = "CrO4²⁻ (Cr VI), Cr³⁺, Cr(OH)2⁺ (truly dissolved)",
    colloidal_forms = "Cr(OH)3 colloids, Cr-organic colloids",
    particulate_forms = "Cr2O3 particles, Cr-Fe oxide particles",
    bioavail_interpretation = "Mixed pattern, Cr VI vs Cr III speciation critical"
  ),
  Cu = list(
    tissue = "Cu", 
    water = c("UF_Cu_NH3", "F_Cu_NH3", "W_Cu_NH3"),
    dissolved_forms = "Cu²⁺, CuOH⁺, small Cu-organic complexes (truly dissolved)",
    colloidal_forms = "Cu-humic/fulvic acid colloids, Cu-Fe oxide colloids",
    particulate_forms = "CuS particles, Cu-Fe/Mn oxide particles",
    bioavail_interpretation = "F>UF>W suggests dissolved + colloidal uptake"
  ),
  Fe = list(
    tissue = "Fer.1", 
    water = c("UF_Fe_NH3", "F_Fe_NH3", "W_Fe_NH3"),
    dissolved_forms = "Fe³⁺, Fe(OH)2⁺ (truly dissolved, limited)",
    colloidal_forms = "Fe(OH)3 colloids (dominant form)",
    particulate_forms = "Fe(OH)3 particles, Fe2O3 particles, Fe-clay particles",
    bioavail_interpretation = "Weak correlations suggest biological regulation"
  ),
  Mn = list(
    tissue = "Mn", 
    water = c("UF_Mn", "F_Mn", "W_Mn"),
    dissolved_forms = "Mn²⁺, MnOH⁺ (truly dissolved)",
    colloidal_forms = "MnO2 colloids, Mn-organic colloids",
    particulate_forms = "MnO2 particles, Mn-Fe oxide particles",
    bioavail_interpretation = "UF>F>W suggests dissolved-only uptake preference"
  ),
  Mo = list(
    tissue = "Mo", 
    water = c("UF_Mo_Oshift", "F_Mo_Oshift", "W_Mo_Oshift"),
    dissolved_forms = "MoO4²⁻, HMoO4⁻ (highly soluble)",
    colloidal_forms = "Mo-organic colloids (minimal)",
    particulate_forms = "Mo-Fe/Mn oxide particles (very limited)",
    bioavail_interpretation = "All fractions similar - molybdate highly bioavailable"
  ),
  Ni = list(
    tissue = "Ni", 
    water = c("UF_Ni", "F_Ni", "W_Ni"),
    dissolved_forms = "Ni²⁺, NiOH⁺ (truly dissolved)",
    colloidal_forms = "Ni-organic colloids, Ni-Fe oxide colloids",
    particulate_forms = "NiS particles, Ni-Fe/Mn oxide particles",
    bioavail_interpretation = "Weak correlations suggest biological regulation"
  ),
  Pb = list(
    tissue = "Pb", 
    water = c("UF_Pb_NH3", "F_Pb_NH3", "W_Pb_NH3"),
    dissolved_forms = "Pb²⁺, PbOH⁺, PbCO3 (truly dissolved)",
    colloidal_forms = "Pb-humic acid colloids, Pb-Fe/Mn oxide colloids",
    particulate_forms = "PbS particles, PbCO3 particles, Pb-Fe oxide particles",
    bioavail_interpretation = "UF>F>W suggests dissolved-only uptake preference"
  ),
  Se = list(
    tissue = "Se", 
    water = c("UF_Se82_Oshift", "F_Se82_Oshift", "W_Se82_Oshift"),
    dissolved_forms = "SeO4²⁻, SeO3²⁻, HSe⁻ (truly dissolved)",
    colloidal_forms = "Se-organic colloids, Se-Fe oxide colloids",
    particulate_forms = "Se⁰ particles, Se-Fe/S compound particles",
    bioavail_interpretation = "UF>F>W suggests dissolved-only uptake preference"
  ),
  Zn = list(
    tissue = "Zn", 
    water = c("UF_Zn", "F_Zn", "W_Zn"),
    dissolved_forms = "Zn²⁺, ZnOH⁺, ZnCO3 (truly dissolved)",
    colloidal_forms = "Zn-organic colloids, Zn-carbonate colloids",
    particulate_forms = "ZnS particles, ZnCO3 particles, Zn-Fe oxide particles",
    bioavail_interpretation = "W>F>UF suggests particulate uptake pathways"
  )
)

main_compartments <- c("EPIL", "EPIP", "FILA")

# Function to analyze compartment interactions for a single metal
analyze_compartment_interactions <- function(data, metal_name, config) {
  
  tissue_col <- config$tissue
  
  # Clean and prepare data
  clean_data <- data %>%
    filter(!is.na(.data[[tissue_col]]) & 
             SAMPLE_DESCRIPTOR %in% main_compartments &
             !is.na(SITE) & !is.na(SAMPLING_DATE) &
             !is.na(TOTAL.MG.M2.RITCHIE)) %>%  # Biomass proxy for abundance
    select(SITE, SAMPLING_DATE, SAMPLE_DESCRIPTOR, 
           all_of(tissue_col), TOTAL.MG.M2.RITCHIE)
  
  if(nrow(clean_data) == 0) return(NULL)
  
  # Average by site/date/compartment to eliminate pseudoreplication
  averaged_data <- clean_data %>%
    group_by(SITE, SAMPLING_DATE, SAMPLE_DESCRIPTOR) %>%
    summarise(
      metal_conc = mean(.data[[tissue_col]], na.rm = TRUE),
      biomass = mean(TOTAL.MG.M2.RITCHIE, na.rm = TRUE),
      n_samples = n(),
      .groups = 'drop'
    ) %>%
    filter(is.finite(metal_conc) & is.finite(biomass))
  
  # Reshape to wide format for interaction analysis
  # Each row = site/date, columns = compartment concentrations and biomasses
  wide_data <- averaged_data %>%
    pivot_wider(
      id_cols = c(SITE, SAMPLING_DATE),
      names_from = SAMPLE_DESCRIPTOR,
      values_from = c(metal_conc, biomass),
      names_sep = "_"
    ) %>%
    # Only keep complete cases (all 3 compartments present)
    filter(!is.na(metal_conc_EPIL) & !is.na(metal_conc_EPIP) & !is.na(metal_conc_FILA) &
             !is.na(biomass_EPIL) & !is.na(biomass_EPIP) & !is.na(biomass_FILA))
  
  if(nrow(wide_data) < 5) return(NULL)  # Need minimum observations
  
  # Calculate interaction effects
  interaction_results <- data.frame()
  
  # Test each compartment as target vs others as predictors
  compartment_pairs <- list(
    list(target = "EPIL", predictors = c("EPIP", "FILA")),
    list(target = "EPIP", predictors = c("EPIL", "FILA")),
    list(target = "FILA", predictors = c("EPIL", "EPIP"))
  )
  
  for(pair in compartment_pairs) {
    target_comp <- pair$target
    target_conc <- paste0("metal_conc_", target_comp)
    target_biomass <- paste0("biomass_", target_comp)
    
    for(pred_comp in pair$predictors) {
      pred_conc <- paste0("metal_conc_", pred_comp)
      pred_biomass <- paste0("biomass_", pred_comp)
      
      # Test biomass effect (abundance interaction)
      if(all(c(target_conc, pred_biomass) %in% names(wide_data))) {
        
        # Pearson correlation between target metal and predictor biomass
        pairs_data <- wide_data %>%
          filter(!is.na(.data[[target_conc]]) & !is.na(.data[[pred_biomass]]) &
                   is.finite(.data[[target_conc]]) & is.finite(.data[[pred_biomass]]))
        
        if(nrow(pairs_data) >= 5) {
          cor_result <- cor.test(pairs_data[[target_conc]], 
                                 pairs_data[[pred_biomass]], 
                                 method = "pearson")
          
          # Calculate percent effect using biomass quartiles
          q1_biomass <- quantile(pairs_data[[pred_biomass]], 0.25, na.rm = TRUE)
          q3_biomass <- quantile(pairs_data[[pred_biomass]], 0.75, na.rm = TRUE)
          
          low_biomass_metal <- mean(pairs_data[[target_conc]][pairs_data[[pred_biomass]] <= q1_biomass], na.rm = TRUE)
          high_biomass_metal <- mean(pairs_data[[target_conc]][pairs_data[[pred_biomass]] >= q3_biomass], na.rm = TRUE)
          
          percent_effect <- if(!is.na(low_biomass_metal) && low_biomass_metal > 0) {
            ((high_biomass_metal - low_biomass_metal) / low_biomass_metal) * 100
          } else { NA }
          
          interaction_type <- case_when(
            is.na(percent_effect) ~ "Indeterminate",
            abs(percent_effect) < 5 ~ "Neutral",
            percent_effect > 5 ~ "Facilitation", 
            percent_effect < -5 ~ "Competition",
            TRUE ~ "Weak"
          )
          
          interaction_results <- rbind(interaction_results, data.frame(
            Target_Compartment = target_comp,
            Predictor_Compartment = pred_comp,
            Interaction_Type = "Biomass_Effect",
            Correlation = cor_result$estimate,
            P_Value = cor_result$p.value,
            N = nrow(pairs_data),
            Percent_Effect = percent_effect,
            Interaction_Pattern = interaction_type,
            Effect_Direction = ifelse(is.na(cor_result$estimate), "None", 
                                      ifelse(cor_result$estimate > 0, "Positive", "Negative")),
            Chemical_Context = config$bioavail_interpretation,
            stringsAsFactors = FALSE
          ))
        }
      }
      
      # Test metal concentration effect (chemical interaction)
      if(all(c(target_conc, pred_conc) %in% names(wide_data))) {
        
        pairs_data <- wide_data %>%
          filter(!is.na(.data[[target_conc]]) & !is.na(.data[[pred_conc]]) &
                   is.finite(.data[[target_conc]]) & is.finite(.data[[pred_conc]]))
        
        if(nrow(pairs_data) >= 5) {
          cor_result <- cor.test(pairs_data[[target_conc]], 
                                 pairs_data[[pred_conc]], 
                                 method = "pearson")
          
          # Calculate percent effect using concentration quartiles
          q1_conc <- quantile(pairs_data[[pred_conc]], 0.25, na.rm = TRUE)
          q3_conc <- quantile(pairs_data[[pred_conc]], 0.75, na.rm = TRUE)
          
          low_conc_metal <- mean(pairs_data[[target_conc]][pairs_data[[pred_conc]] <= q1_conc], na.rm = TRUE)
          high_conc_metal <- mean(pairs_data[[target_conc]][pairs_data[[pred_conc]] >= q3_conc], na.rm = TRUE)
          
          percent_effect <- if(!is.na(low_conc_metal) && low_conc_metal > 0) {
            ((high_conc_metal - low_conc_metal) / low_conc_metal) * 100
          } else { NA }
          
          interaction_type <- case_when(
            is.na(percent_effect) ~ "Indeterminate",
            abs(percent_effect) < 5 ~ "Neutral",
            percent_effect > 5 ~ "Facilitation",
            percent_effect < -5 ~ "Competition", 
            TRUE ~ "Weak"
          )
          
          interaction_results <- rbind(interaction_results, data.frame(
            Target_Compartment = target_comp,
            Predictor_Compartment = pred_comp,
            Interaction_Type = "Metal_Conc_Effect",
            Correlation = cor_result$estimate,
            P_Value = cor_result$p.value,
            N = nrow(pairs_data),
            Percent_Effect = percent_effect,
            Interaction_Pattern = interaction_type,
            Effect_Direction = ifelse(is.na(cor_result$estimate), "None",
                                      ifelse(cor_result$estimate > 0, "Positive", "Negative")),
            Chemical_Context = config$bioavail_interpretation,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Calculate summary statistics
  n_complete_triads <- nrow(wide_data)
  
  significant_interactions <- interaction_results %>%
    filter(P_Value < 0.05 & abs(Correlation) > 0.3)
  
  competition_effects <- interaction_results %>%
    filter(Interaction_Pattern == "Competition")
  
  facilitation_effects <- interaction_results %>%
    filter(Interaction_Pattern == "Facilitation")
  
  strongest_competition <- if(nrow(competition_effects) > 0) {
    competition_effects[which.min(competition_effects$Percent_Effect), ]
  } else { NULL }
  
  strongest_facilitation <- if(nrow(facilitation_effects) > 0) {
    facilitation_effects[which.max(facilitation_effects$Percent_Effect), ]
  } else { NULL }
  
  # Return results with chemical context
  list(
    n_complete_triads = n_complete_triads,
    interaction_results = interaction_results,
    n_significant = nrow(significant_interactions),
    n_competition = nrow(competition_effects),
    n_facilitation = nrow(facilitation_effects),
    strongest_competition = strongest_competition,
    strongest_facilitation = strongest_facilitation,
    mean_biomass_by_comp = wide_data %>%
      summarise(
        mean_EPIL_biomass = mean(biomass_EPIL, na.rm = TRUE),
        mean_EPIP_biomass = mean(biomass_EPIP, na.rm = TRUE),
        mean_FILA_biomass = mean(biomass_FILA, na.rm = TRUE)
      ),
    chemical_info = config[c("dissolved_forms", "colloidal_forms", "particulate_forms", "bioavail_interpretation")]
  )
}

# Run analysis for all metals
cat("Running Hypothesis 2: Compartment Interaction Analysis...\n")

interaction_results <- list()

for(metal_name in names(metal_configs)) {
  config <- metal_configs[[metal_name]]
  
  if(!config$tissue %in% names(metal_data)) {
    cat("Skipping", metal_name, "- tissue column not found\n")
    next
  }
  
  cat("Analyzing compartment interactions for", metal_name, "...\n")
  result <- analyze_compartment_interactions(metal_data, metal_name, config)
  
  if(!is.null(result)) {
    interaction_results[[metal_name]] <- result
  }
}

# Create comprehensive summary with chemical context
create_interaction_summary <- function(results) {
  
  summary_data <- data.frame()
  
  for(metal_name in names(results)) {
    result <- results[[metal_name]]
    
    # Extract key metrics
    strongest_comp <- if(!is.null(result$strongest_competition)) {
      paste0(result$strongest_competition$Predictor_Compartment, "→", 
             result$strongest_competition$Target_Compartment, " (",
             round(result$strongest_competition$Percent_Effect, 1), "%)")
    } else { "None" }
    
    strongest_facil <- if(!is.null(result$strongest_facilitation)) {
      paste0(result$strongest_facilitation$Predictor_Compartment, "→",
             result$strongest_facilitation$Target_Compartment, " (+",
             round(result$strongest_facilitation$Percent_Effect, 1), "%)")
    } else { "None" }
    
    summary_data <- rbind(summary_data, data.frame(
      Metal = metal_name,
      N_Complete_Triads = result$n_complete_triads,
      N_Significant_Interactions = result$n_significant,
      N_Competition_Effects = result$n_competition,
      N_Facilitation_Effects = result$n_facilitation,
      Strongest_Competition = strongest_comp,
      Strongest_Facilitation = strongest_facil,
      Mean_EPIL_Biomass = round(result$mean_biomass_by_comp$mean_EPIL_biomass, 2),
      Mean_EPIP_Biomass = round(result$mean_biomass_by_comp$mean_EPIP_biomass, 2),
      Mean_FILA_Biomass = round(result$mean_biomass_by_comp$mean_FILA_biomass, 2),
      Dissolved_Forms = result$chemical_info$dissolved_forms,
      Bioavail_Interpretation = result$chemical_info$bioavail_interpretation,
      stringsAsFactors = FALSE
    ))
  }
  
  summary_data
}

# Create detailed interaction matrix with chemical context
create_interaction_matrix <- function(results) {
  
  all_interactions <- data.frame()
  
  for(metal_name in names(results)) {
    result <- results[[metal_name]]
    
    if(nrow(result$interaction_results) > 0) {
      metal_interactions <- result$interaction_results %>%
        mutate(Metal = metal_name) %>%
        select(Metal, Target_Compartment, Predictor_Compartment, 
               Interaction_Type, Correlation, P_Value, N,
               Percent_Effect, Interaction_Pattern, Effect_Direction,
               Chemical_Context)
      
      all_interactions <- rbind(all_interactions, metal_interactions)
    }
  }
  
  all_interactions
}

# Generate results
interaction_summary <- create_interaction_summary(interaction_results)
interaction_matrix <- create_interaction_matrix(interaction_results)

# Display results
cat("\n=== HYPOTHESIS 2: COMPARTMENT INTERACTION SUMMARY ===\n")
print(interaction_summary, row.names = FALSE)

cat("\n=== DETAILED INTERACTION MATRIX ===\n")
# Show only significant interactions (p < 0.05, |r| > 0.3)
significant_interactions <- interaction_matrix %>%
  filter(P_Value < 0.05 & abs(Correlation) > 0.3) %>%
  arrange(Metal, Target_Compartment, desc(abs(Percent_Effect)))

if(nrow(significant_interactions) > 0) {
  print(significant_interactions, row.names = FALSE)
} else {
  cat("No significant interactions found (p < 0.05, |r| > 0.3)\n")
}

# Chemical speciation analysis of interactions
cat("\n=== INTERACTIONS BY CHEMICAL SPECIATION ===\n")
speciation_interaction_analysis <- significant_interactions %>%
  mutate(
    Primary_Chemical_Form = case_when(
      grepl("MoO4", Chemical_Context) ~ "Molybdate (highly soluble)",
      grepl("AsO", Chemical_Context) ~ "Arsenate/Arsenite oxyanions", 
      grepl("dissolved \\+ colloidal", Chemical_Context) ~ "Mixed dissolved + colloidal",
      grepl("dissolved-only", Chemical_Context) ~ "Dissolved-only preference",
      grepl("particulate", Chemical_Context) ~ "Particulate pathways",
      grepl("biological regulation", Chemical_Context) ~ "Biologically regulated",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(Primary_Chemical_Form) %>%
  summarise(
    N_Interactions = n(),
    Metals = paste(unique(Metal), collapse = ", "),
    Avg_Effect_Size = round(mean(abs(Percent_Effect), na.rm = TRUE), 1),
    Max_Effect = round(max(abs(Percent_Effect), na.rm = TRUE), 1),
    Interaction_Types = paste(unique(Interaction_Pattern), collapse = ", "),
    .groups = 'drop'
  ) %>%
  arrange(desc(Avg_Effect_Size))

print(speciation_interaction_analysis, row.names = FALSE)

# Competition/Facilitation patterns by metal speciation
cat("\n=== COMPETITION vs FACILITATION BY SPECIATION ===\n")
pattern_by_speciation <- interaction_matrix %>%
  filter(Interaction_Pattern %in% c("Competition", "Facilitation")) %>%
  mutate(
    Speciation_Category = case_when(
      grepl("molybdate", Chemical_Context) ~ "Highly soluble oxyanion",
      grepl("dissolved \\+ colloidal", Chemical_Context) ~ "Mixed uptake strategy",
      grepl("dissolved-only", Chemical_Context) ~ "Dissolved-only strategy", 
      grepl("particulate", Chemical_Context) ~ "Particulate strategy",
      grepl("biological regulation", Chemical_Context) ~ "Biologically regulated",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(Speciation_Category, Interaction_Pattern) %>%
  summarise(
    Count = n(),
    Metals = paste(unique(Metal), collapse = ", "),
    Mean_Effect = round(mean(abs(Percent_Effect), na.rm = TRUE), 1),
    .groups = 'drop'
  ) %>%
  pivot_wider(names_from = Interaction_Pattern, 
              values_from = c(Count, Mean_Effect),
              values_fill = 0)

print(pattern_by_speciation, row.names = FALSE)

# Save results
write.csv(interaction_summary, "2_incremental/7_AI_ANALYSIS/hypothesis2_interaction_summary.csv", row.names = FALSE)
write.csv(interaction_matrix, "2_incremental/7_AI_ANALYSIS/hypothesis2_interaction_matrix.csv", row.names = FALSE)
write.csv(significant_interactions, "2_incremental/7_AI_ANALYSIS/hypothesis2_significant_interactions.csv", row.names = FALSE)
write.csv(speciation_interaction_analysis, "2_incremental/7_AI_ANALYSIS/hypothesis2_speciation_analysis.csv", row.names = FALSE)

cat("\nHypothesis 2 Analysis Complete! Files saved:\n")
cat("- hypothesis2_interaction_summary.csv\n")
cat("- hypothesis2_interaction_matrix.csv\n") 
cat("- hypothesis2_significant_interactions.csv\n")
cat("- hypothesis2_speciation_analysis.csv\n")

# Print key findings summary
cat("\n=== KEY FINDINGS SUMMARY ===\n")
total_triads <- sum(interaction_summary$N_Complete_Triads)
total_significant <- sum(interaction_summary$N_Significant_Interactions)
total_competition <- sum(interaction_summary$N_Competition_Effects)
total_facilitation <- sum(interaction_summary$N_Facilitation_Effects)

cat("Total complete site/date triads (all 3 compartments present):", total_triads, "\n")
cat("Total significant interactions (p<0.05, |r|>0.3):", total_significant, "\n")
cat("Competition effects:", total_competition, "\n")
cat("Facilitation effects:", total_facilitation, "\n")

if(nrow(significant_interactions) > 0) {
  cat("\nMost notable interactions:\n")
  notable <- significant_interactions %>%
    filter(abs(Percent_Effect) > 10) %>%
    arrange(desc(abs(Percent_Effect))) %>%
    head(5)
  
  for(i in 1:nrow(notable)) {
    cat(sprintf("- %s (%s): %s %s affects %s (%+.1f%%, r=%.3f, p=%.3f)\n",
                notable$Metal[i],
                notable$Chemical_Context[i],
                notable$Predictor_Compartment[i],
                notable$Interaction_Type[i],
                notable$Target_Compartment[i],
                notable$Percent_Effect[i],
                notable$Correlation[i],
                notable$P_Value[i]))
  }
}