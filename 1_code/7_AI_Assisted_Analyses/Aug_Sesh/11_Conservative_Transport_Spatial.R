# Separate code chunk: Temporal Deviations for All Three Cd Fractions
# Assumes all objects from the conservative transport model are available

# Set up plotting layout (1 row, 3 columns)
par(mfrow=c(1,3))

# Determine which date format to use
if(exists("date_char", where=reach_data)) {
  date_column <- reach_data$date_char
  date_label <- "Date"
} else {
  date_column <- reach_data$date  
  date_label <- "Date"
}

# Plot 1: Dissolved Cd Deviations by Date
valid_dissolved <- reach_data[!is.na(reach_data$dissolved_rel_dev) & is.finite(reach_data$dissolved_rel_dev), ]
if(nrow(valid_dissolved) > 0) {
  if(exists("date_char", where=reach_data)) {
    boxplot(valid_dissolved$dissolved_rel_dev ~ valid_dissolved$date_char, 
            main="Dissolved Cd Deviations", 
            xlab=date_label, 
            ylab="Relative Deviation (%)",
            las=2, col="lightblue")
  } else {
    boxplot(valid_dissolved$dissolved_rel_dev ~ valid_dissolved$date, 
            main="Dissolved Cd Deviations", 
            xlab=date_label, 
            ylab="Relative Deviation (%)",
            las=2, col="lightblue")
  }
  abline(h=0, lty=2, col="gray")
  
  # Add sample size info
  n_dissolved <- length(valid_dissolved$dissolved_rel_dev)
  mtext(paste("n =", n_dissolved), side=3, line=0, adj=1, cex=0.8)
} else {
  plot(1, type="n", xlab=date_label, ylab="Relative Deviation (%)", main="Dissolved Cd Deviations")
  text(1, 1, "No valid data", cex=1.2)
}

# Plot 2: Colloidal Cd Deviations by Date
valid_colloidal <- reach_data[!is.na(reach_data$colloidal_rel_dev) & is.finite(reach_data$colloidal_rel_dev), ]
if(nrow(valid_colloidal) > 0) {
  if(exists("date_char", where=reach_data)) {
    boxplot(valid_colloidal$colloidal_rel_dev ~ valid_colloidal$date_char, 
            main="Colloidal Cd Deviations", 
            xlab=date_label, 
            ylab="Relative Deviation (%)",
            las=2, col="lightgreen")
  } else {
    boxplot(valid_colloidal$colloidal_rel_dev ~ valid_colloidal$date, 
            main="Colloidal Cd Deviations", 
            xlab=date_label, 
            ylab="Relative Deviation (%)",
            las=2, col="lightgreen")
  }
  abline(h=0, lty=2, col="gray")
  
  # Add sample size info
  n_colloidal <- length(valid_colloidal$colloidal_rel_dev)
  mtext(paste("n =", n_colloidal), side=3, line=0, adj=1, cex=0.8)
} else {
  plot(1, type="n", xlab=date_label, ylab="Relative Deviation (%)", main="Colloidal Cd Deviations")
  text(1, 1, "No valid data", cex=1.2)
}

# Plot 3: Suspended Cd Deviations by Date
valid_suspended <- reach_data[!is.na(reach_data$suspended_rel_dev) & is.finite(reach_data$suspended_rel_dev), ]
if(nrow(valid_suspended) > 0) {
  if(exists("date_char", where=reach_data)) {
    boxplot(valid_suspended$suspended_rel_dev ~ valid_suspended$date_char, 
            main="Suspended Cd Deviations", 
            xlab=date_label, 
            ylab="Relative Deviation (%)",
            las=2, col="lightcoral")
  } else {
    boxplot(valid_suspended$suspended_rel_dev ~ valid_suspended$date, 
            main="Suspended Cd Deviations", 
            xlab=date_label, 
            ylab="Relative Deviation (%)",
            las=2, col="lightcoral")
  }
  abline(h=0, lty=2, col="gray")
  
  # Add sample size info
  n_suspended <- length(valid_suspended$suspended_rel_dev)
  mtext(paste("n =", n_suspended), side=3, line=0, adj=1, cex=0.8)
} else {
  plot(1, type="n", xlab=date_label, ylab="Relative Deviation (%)", main="Suspended Cd Deviations")
  text(1, 1, "No valid data", cex=1.2)
}

# Add overall title
mtext("Temporal Patterns in Conservative Transport Deviations", side=3, line=-2, outer=TRUE, cex=1.2, font=2)

# Print summary statistics for temporal analysis
cat("\n=== TEMPORAL DEVIATION SUMMARY ===\n")
cat("Sample sizes by fraction:\n")
cat("Dissolved:", nrow(valid_dissolved), "data points\n")
cat("Colloidal:", nrow(valid_colloidal), "data points\n") 
cat("Suspended:", nrow(valid_suspended), "data points\n")

if(nrow(valid_dissolved) > 0) {
  cat("\nDissolved fraction temporal stats:\n")
  cat("Mean absolute deviation:", round(mean(abs(valid_dissolved$dissolved_rel_dev)), 1), "%\n")
  cat("Range:", round(range(valid_dissolved$dissolved_rel_dev), 1), "%\n")
}

if(nrow(valid_colloidal) > 0) {
  cat("\nColloidal fraction temporal stats:\n")
  cat("Mean absolute deviation:", round(mean(abs(valid_colloidal$colloidal_rel_dev)), 1), "%\n")
  cat("Range:", round(range(valid_colloidal$colloidal_rel_dev), 1), "%\n")
}

if(nrow(valid_suspended) > 0) {
  cat("\nSuspended fraction temporal stats:\n")
  cat("Mean absolute deviation:", round(mean(abs(valid_suspended$suspended_rel_dev)), 1), "%\n")
  cat("Range:", round(range(valid_suspended$suspended_rel_dev), 1), "%\n")
}

# Reset plotting parameters
par(mfrow=c(1,1))