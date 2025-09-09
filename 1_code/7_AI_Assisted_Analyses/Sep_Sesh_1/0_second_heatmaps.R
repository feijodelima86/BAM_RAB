# Simple 3-panel heatmap: Flux, Stock, Percent Storage
# Sites ordered by distance (y-axis), Time (x-axis)

# Set up 3-panel layout
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

# Extract and sort dates from matrix rownames
dates_num <- as.numeric(rownames(reach_uf_flux))
date_order <- order(dates_num)
dates_sorted <- dates_num[date_order]

# Reorder matrices by sorted dates
flux_sorted <- reach_uf_flux[date_order, ]
stock_sorted <- reach_stock[date_order, ]
percent_sorted <- percent_storage[date_order, ]

# Apply log+1 transformation to all values
flux_log <- log(flux_sorted + 1)
stock_log <- log(stock_sorted + 1)
percent_log <- log(percent_sorted + 1)

blue_yellow <- colorRampPalette(c("blue", "yellow"))(20)

reach_positions <- 1:ncol(flux_sorted)

# Panel 1: Flux
image(dates_sorted, reach_positions, flux_log, 
      col = blue_yellow, 
      main = "log(UF Cd Flux + 1)",
      xlab = "Date", ylab = "River Reach",
      yaxt = "n")
axis(2, at = reach_positions, labels = colnames(flux_sorted), las = 2)

# Panel 2: Stock  
image(dates_sorted, reach_positions, stock_log,
      col = blue_yellow,
      main = "log(Cd Stock + 1)", 
      xlab = "Date", ylab = "",
      yaxt = "n")
axis(2, at = reach_positions, labels = colnames(stock_sorted), las = 2)

# Panel 3: Percent Storage
image(dates_sorted, reach_positions, percent_log,
      col = blue_yellow,
      main = "log(% Storage + 1)",
      xlab = "Date", ylab = "",
      yaxt = "n") 
axis(2, at = reach_positions, labels = colnames(percent_sorted), las = 2)

# Reset layout
par(mfrow = c(1, 1))