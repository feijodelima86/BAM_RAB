# Load required libraries
library(corrplot)
library(dplyr)

# Read the data (assuming it's saved as 'data.csv')
data <- read.csv("2_incremental/AI_Subsets/integrated_parsed.csv")

# Select only numeric variables (excluding dates and text)
numeric_vars <- data %>%
  select(where(is.numeric)) %>%
  # Remove any variables that might be IDs or row numbers if needed
  select(-matches("^X$|^X\\d+$"))  # This removes columns like X, X1, X2 that might be row indices

# Calculate correlation matrix
cor_matrix <- cor(numeric_vars, use = "complete.obs")

# Create correlogram
corrplot(cor_matrix, 
         method = "color",           # Use colors to show correlations
         type = "upper",            # Show upper triangle only
         order = "hclust",          # Order variables by hierarchical clustering
         tl.cex = 0.7,             # Text label size
         tl.col = "black",         # Text label color
         tl.srt = 45,              # Text label rotation
         diag = FALSE)             # Don't show diagonal

# Alternative: Show correlation coefficients as numbers
corrplot(cor_matrix, 
         method = "number",
         type = "upper",
         order = "hclust",
         tl.cex = 0.7,
         tl.col = "black",
         tl.srt = 45,
         number.cex = 0.6,
         diag = FALSE)

# Optional: If you want a more detailed correlogram with both colors and numbers
corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         order = "hclust",
         tl.cex = 0.7,
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "white",    # Add correlation coefficients in white
         number.cex = 0.5,
         diag = FALSE)