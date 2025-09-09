# Separate code chunk: Temporal trends of the three Cd fractions (dissolved, colloidal, suspended) over time
# Assumes all objects from the conservative transport model are available

# Determine which date format to use and create proper date objects for plotting
if(exists("date_char", where=reach_data)) {
  reach_data$plot_date <- as.Date(reach_data$date_char)
  date_label <- "Date"
} else {
  # Convert Excel serial dates if needed
  reach_data$plot_date <- as.Date(reach_data$date, origin="1900-01-01") - 2
  date_label <- "Date"
}

# Create summary data by date for each fraction
dates <- sort(unique(reach_data$plot_date))

# Calculate mean deviations by date for each fraction
dissolved_by_date <- aggregate(dissolved_rel_dev ~ plot_date, data=reach_data, FUN=function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x[!is.na(x)])))
colloidal_by_date <- aggregate(colloidal_rel_dev ~ plot_date, data=reach_data, FUN=function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x[!is.na(x)])))
suspended_by_date <- aggregate(suspended_rel_dev ~ plot_date, data=reach_data, FUN=function(x) c(mean=mean(x, na.rm=T), sd=sd(x, na.rm=T), n=length(x[!is.na(x)])))

# Extract mean values for plotting
dissolved_means <- dissolved_by_date$dissolved_rel_dev[,1]
colloidal_means <- colloidal_by_date$colloidal_rel_dev[,1]
suspended_means <- suspended_by_date$suspended_rel_dev[,1]

# Extract standard deviations for error bars
dissolved_sds <- dissolved_by_date$dissolved_rel_dev[,2]
colloidal_sds <- colloidal_by_date$colloidal_rel_dev[,2]
suspended_sds <- suspended_by_date$suspended_rel_dev[,2]

# Set up plotting layout (2 rows, 2 columns)
par(mfrow=c(2,2))

# Plot 1: All three fractions on one plot - Mean trends over time
y_range <- range(c(dissolved_means, colloidal_means, suspended_means), na.rm=T)
y_range <- y_range + c(-0.1, 0.1) * diff(y_range)  # Add some padding

plot(dissolved_by_date$plot_date, dissolved_means, 
     type="o", pch=16, col="blue", lwd=2,
     xlab=date_label, ylab="Mean Relative Deviation (%)",
     main="Temporal Trends: All Cd Fractions", ylim=y_range)

lines(colloidal_by_date$plot_date, colloidal_means, type="o", pch=17, col="green", lwd=2)
lines(suspended_by_date$plot_date, suspended_means, type="o", pch=18, col="red", lwd=2)

abline(h=0, lty=2, col="gray")
legend("topright", c("Dissolved", "Colloidal", "Suspended"), 
       col=c("blue", "green", "red"), pch=c(16, 17, 18), lwd=2)

# Plot 2: Dissolved fraction with error bars
plot(dissolved_by_date$plot_date, dissolved_means, 
     type="o", pch=16, col="blue", lwd=2,
     xlab=date_label, ylab="Relative Deviation (%)",
     main="Dissolved Cd: Mean ± SD",
     ylim=range(c(dissolved_means - dissolved_sds, dissolved_means + dissolved_sds), na.rm=T))

# Add error bars
for(i in 1:length(dissolved_means)) {
  if(!is.na(dissolved_means[i]) & !is.na(dissolved_sds[i])) {
    lines(c(dissolved_by_date$plot_date[i], dissolved_by_date$plot_date[i]), 
          c(dissolved_means[i] - dissolved_sds[i], dissolved_means[i] + dissolved_sds[i]), 
          col="blue", lwd=1)
  }
}
abline(h=0, lty=2, col="gray")

# Plot 3: Colloidal fraction with error bars
plot(colloidal_by_date$plot_date, colloidal_means, 
     type="o", pch=17, col="green", lwd=2,
     xlab=date_label, ylab="Relative Deviation (%)",
     main="Colloidal Cd: Mean ± SD",
     ylim=range(c(colloidal_means - colloidal_sds, colloidal_means + colloidal_sds), na.rm=T))

# Add error bars
for(i in 1:length(colloidal_means)) {
  if(!is.na(colloidal_means[i]) & !is.na(colloidal_sds[i])) {
    lines(c(colloidal_by_date$plot_date[i], colloidal_by_date$plot_date[i]), 
          c(colloidal_means[i] - colloidal_sds[i], colloidal_means[i] + colloidal_sds[i]), 
          col="green", lwd=1)
  }
}
abline(h=0, lty=2, col="gray")

# Plot 4: Suspended fraction with error bars
plot(suspended_by_date$plot_date, suspended_means, 
     type="o", pch=18, col="red", lwd=2,
     xlab=date_label, ylab="Relative Deviation (%)",
     main="Suspended Cd: Mean ± SD",
     ylim=range(c(suspended_means - suspended_sds, suspended_means + suspended_sds), na.rm=T))

# Add error bars
for(i in 1:length(suspended_means)) {
  if(!is.na(suspended_means[i]) & !is.na(suspended_sds[i])) {
    lines(c(suspended_by_date$plot_date[i], suspended_by_date$plot_date[i]), 
          c(suspended_means[i] - suspended_sds[i], suspended_means[i] + suspended_sds[i]), 
          col="red", lwd=1)
  }
}
abline(h=0, lty=2, col="gray")

# Print temporal analysis summary
cat("\n=== TEMPORAL TRENDS BY CD FRACTION ===\n")

# Calculate overall temporal trends (slope analysis)
if(length(dates) > 2) {
  
  # Linear trends for each fraction
  time_numeric <- as.numeric(dissolved_by_date$plot_date)
  
  if(length(dissolved_means[!is.na(dissolved_means)]) > 2) {
    dissolved_trend <- lm(dissolved_means ~ time_numeric)
    dissolved_slope <- coef(dissolved_trend)[2] * 30  # Convert to % per month
    dissolved_p <- summary(dissolved_trend)$coefficients[2,4]
    
    cat("\nDISSOLVED FRACTION:\n")
    cat("Temporal trend:", round(dissolved_slope, 2), "% per month\n")
    cat("Trend p-value:", round(dissolved_p, 3), ifelse(dissolved_p < 0.05, " *", ""), "\n")
    cat("Mean deviation:", round(mean(dissolved_means, na.rm=T), 1), "±", round(sd(dissolved_means, na.rm=T), 1), "%\n")
  }
  
  if(length(colloidal_means[!is.na(colloidal_means)]) > 2) {
    colloidal_trend <- lm(colloidal_means ~ time_numeric)
    colloidal_slope <- coef(colloidal_trend)[2] * 30
    colloidal_p <- summary(colloidal_trend)$coefficients[2,4]
    
    cat("\nCOLLOIDAL FRACTION:\n")
    cat("Temporal trend:", round(colloidal_slope, 2), "% per month\n")
    cat("Trend p-value:", round(colloidal_p, 3), ifelse(colloidal_p < 0.05, " *", ""), "\n")
    cat("Mean deviation:", round(mean(colloidal_means, na.rm=T), 1), "±", round(sd(colloidal_means, na.rm=T), 1), "%\n")
  }
  
  if(length(suspended_means[!is.na(suspended_means)]) > 2) {
    suspended_trend <- lm(suspended_means ~ time_numeric)
    suspended_slope <- coef(suspended_trend)[2] * 30
    suspended_p <- summary(suspended_trend)$coefficients[2,4]
    
    cat("\nSUSPENDED FRACTION:\n")
    cat("Temporal trend:", round(suspended_slope, 2), "% per month\n")
    cat("Trend p-value:", round(suspended_p, 3), ifelse(suspended_p < 0.05, " *", ""), "\n")
    cat("Mean deviation:", round(mean(suspended_means, na.rm=T), 1), "±", round(sd(suspended_means, na.rm=T), 1), "%\n")
  }
}

# Seasonal analysis (early vs late season)
early_season <- reach_data$plot_date <= as.Date("2021-07-15")
late_season <- reach_data$plot_date >= as.Date("2021-08-15")

cat("\n=== SEASONAL COMPARISON ===\n")

if(sum(early_season, na.rm=T) > 0 & sum(late_season, na.rm=T) > 0) {
  
  early_dissolved <- mean(reach_data$dissolved_rel_dev[early_season], na.rm=T)
  late_dissolved <- mean(reach_data$dissolved_rel_dev[late_season], na.rm=T)
  
  early_colloidal <- mean(reach_data$colloidal_rel_dev[early_season], na.rm=T)
  late_colloidal <- mean(reach_data$colloidal_rel_dev[late_season], na.rm=T)
  
  early_suspended <- mean(reach_data$suspended_rel_dev[early_season], na.rm=T)
  late_suspended <- mean(reach_data$suspended_rel_dev[late_season], na.rm=T)
  
  cat("Early season (≤July 15) vs Late season (≥Aug 15):\n")
  cat("Dissolved:", round(early_dissolved, 1), "% vs", round(late_dissolved, 1), "%\n")
  cat("Colloidal:", round(early_colloidal, 1), "% vs", round(late_colloidal, 1), "%\n")
  cat("Suspended:", round(early_suspended, 1), "% vs", round(late_suspended, 1), "%\n")
  
  # Identify largest seasonal changes
  seasonal_changes <- abs(c(late_dissolved - early_dissolved, 
                            late_colloidal - early_colloidal, 
                            late_suspended - early_suspended))
  fraction_names <- c("Dissolved", "Colloidal", "Suspended")
  max_change_idx <- which.max(seasonal_changes)
  
  cat("\nLargest seasonal change:", fraction_names[max_change_idx], 
      "(" , round(seasonal_changes[max_change_idx], 1), "% difference)\n")
}

# Reset plotting parameters
par(mfrow=c(1,1))
