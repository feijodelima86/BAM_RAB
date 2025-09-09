# Load required libraries
library(reshape2)

# Assume your data is already loaded as:
algae_data <- read.csv("2_incremental/AI_Subsets/Cd_redux.csv")
water_data <- read.csv("2_incremental/AI_Subsets/integrated_parsed.csv")

# =============================================================================
# DIAGNOSE AND HANDLE DATE FORMATS
# =============================================================================

# Check the structure and sample values of date columns
cat("=== DIAGNOSING DATE FORMATS ===\n")
cat("Algae data SAMPLING_DATE:\n")
cat("Class:", class(algae_data$SAMPLING_DATE), "\n")
cat("Sample values:", head(algae_data$SAMPLING_DATE), "\n\n")

cat("Water data Date:\n") 
cat("Class:", class(water_data$Date), "\n")
cat("Sample values:", head(water_data$Date), "\n\n")

# Function to convert dates only if needed
convert_date_if_needed <- function(date_col) {
  # If already Date, return as-is
  if(inherits(date_col, "Date")) {
    cat("Column already in Date format - no conversion needed\n")
    return(date_col)
  }
  
  # Convert factor to character if needed
  if(is.factor(date_col)) {
    date_col <- as.character(date_col)
  }
  
  # Try multiple common formats
  formats_to_try <- c("%m/%d/%Y", "%m/%d/%y", "%Y-%m-%d", "%d/%m/%Y")
  
  for(fmt in formats_to_try) {
    result <- as.Date(date_col, format = fmt)
    if(sum(!is.na(result)) > 0) {
      cat("Successfully converted using format:", fmt, "\n")
      return(result)
    }
  }
  
  # If none work, try automatic parsing
  cat("Trying automatic date parsing...\n")
  return(as.Date(date_col))
}

# Convert dates only if needed
cat("Converting algae SAMPLING_DATE:\n")
algae_data$SAMPLING_DATE <- convert_date_if_needed(algae_data$SAMPLING_DATE)

cat("Converting water Date:\n")
water_data$Date <- convert_date_if_needed(water_data$Date)

# Verify conversion worked
cat("\n=== AFTER CONVERSION ===\n")
cat("Algae dates sample:", head(algae_data$SAMPLING_DATE), "\n")
cat("Water dates sample:", head(water_data$Date), "\n")

# =============================================================================
# DATA FRAME 1: Algae samples with repeated water quality data
# =============================================================================

# Simple merge on site and date
df1 <- merge(algae_data, water_data, 
             by.x = c("SITE", "SAMPLING_DATE"), 
             by.y = c("Site", "Date"), 
             all.x = TRUE)

print("Data Frame 1 created - algae samples with water quality data repeated")
print(paste("Dimensions:", nrow(df1), "rows x", ncol(df1), "columns"))

# =============================================================================
# DATA FRAME 2: Water chemistry with algal attributes by compartment
# =============================================================================

# Step 1: Identify numeric algal columns to aggregate
algal_numeric_cols <- c("OM.AREA.g.m2", "ASH.AREA.g.m2", "PHICOCYANIN.MG.M2", 
                        "TOTAL.MG.M2.RITCHIE", "Cd", "Cd.1", "UF_Cd114", 
                        "F_Cd114", "W_Cd114", "SEP", "AUTO")

# Step 2: Aggregate algae data by Site, Date, and Compartment (taking means)
algae_agg <- aggregate(algae_data[algal_numeric_cols], 
                       by = list(Site = algae_data$SITE, 
                                 Date = algae_data$SAMPLING_DATE, 
                                 Compartment = algae_data$SAMPLE_DESCRIPTOR), 
                       FUN = mean, na.rm = TRUE)

# Step 3: Reshape to wide format - each algal attribute gets columns for each compartment
algae_wide <- reshape(algae_agg, 
                      timevar = "Compartment", 
                      idvar = c("Site", "Date"), 
                      direction = "wide")

# Step 4: Clean up column names (remove "." and replace with "_")
names(algae_wide) <- gsub("\\.", "_", names(algae_wide))

# Step 5: Merge with water chemistry data
df2 <- merge(water_data, algae_wide, 
             by = c("Site", "Date"), 
             all.x = TRUE)

print("Data Frame 2 created - water chemistry with algal attributes by compartment")
print(paste("Dimensions:", nrow(df2), "rows x", ncol(df2), "columns"))

# =============================================================================
# Optional: View structure of both data frames
# =============================================================================

cat("\n=== DATA FRAME 1 STRUCTURE ===\n")
str(df1)

cat("\n=== DATA FRAME 2 STRUCTURE ===\n")
str(df2)

cat("\n=== DATA FRAME 2 COLUMN NAMES (showing algal attributes by compartment) ===\n")
algal_cols_df2 <- names(df2)[grep("_EPIP|_EPIL|_FILA", names(df2))]
print(algal_cols_df2)


