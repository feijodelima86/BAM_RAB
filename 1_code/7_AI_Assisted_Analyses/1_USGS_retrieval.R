# Install and load required packages
if (!require(dataRetrieval)) {
  install.packages("dataRetrieval")
  library(dataRetrieval)
}
if (!require(dplyr)) {
  install.packages("dplyr") 
  library(dplyr)
}

# ===============================================================
# DEFINE USGS STATIONS AND SAMPLING SITES
# ===============================================================

# Main stem Clark Fork stations
main_stem_stations <- c(
  "12323750",  # Silver Bow Creek at Warm Springs MT (WS)
  "12323800",  # Clark Fork near Galen MT
  "12324200",  # Clark Fork at Deer Lodge MT (DL)
  "12324400",  # Clark Fork ab Little Blackfoot R nr Garrison MT (GR)
  "12324680",  # Clark Fork at Goldcreek MT (GC)
  "12331800",  # Clark Fork near Drummond MT
  "12340500"   # Clark Fork above Missoula MT
)

# Tributary stations
tributary_stations <- c(
  "12323850",  # Lost Creek near Galen, MT
  "12324590",  # Little Blackfoot River near Garrison MT
  "12331500",  # Flint Creek near Drummond MT
  "12334510"   # Rock Creek near Clinton MT
)

# All stations combined
all_stations <- c(main_stem_stations, tributary_stations)

# Station metadata
station_info <- data.frame(
  station_id = all_stations,
  station_type = c(rep("main_stem", length(main_stem_stations)), 
                   rep("tributary", length(tributary_stations))),
  sampling_site = c("WS", "Between_WS_DL", "DL", "GR", "GC", "Between_GC_BG", "Between_BN_Missoula",
                    "Lost_Creek", "Little_Blackfoot", "Flint_Creek", "Rock_Creek"),
  confluence_location = c("WS", "Before_DL", "DL", "GR", "GC", "After_GC", "Before_Missoula",
                          "At_WS", "After_GR", "After_GC", "Between_BG_BN")
)

# ===============================================================
# DOWNLOAD USGS DATA
# ===============================================================

# Date range for your study
start_date <- "2021-06-01"
end_date <- "2021-10-31"

# Parameters to download
# 00060 = Discharge (cfs)
# 00010 = Temperature (°C)
# 00095 = Specific conductance (μS/cm)
# 00400 = pH (if available)
parameters <- c("00060", "00010", "00095", "00400")

# Function to download data with error handling
download_usgs_data <- function(station_ids, start_date, end_date, parameters) {
  all_data <- list()
  
  for(station in station_ids) {
    cat("Downloading data for station:", station, "\n")
    
    tryCatch({
      # Get daily values
      daily_data <- readNWISdv(siteNumbers = station,
                               parameterCd = parameters,
                               startDate = start_date,
                               endDate = end_date)
      
      # Get site information
      site_info <- readNWISsite(siteNumbers = station)
      
      # Add site metadata to data
      if(nrow(daily_data) > 0) {
        daily_data$station_name <- site_info$station_nm[1]
        daily_data$latitude <- site_info$dec_lat_va[1]
        daily_data$longitude <- site_info$dec_long_va[1]
        daily_data$drainage_area_sqmi <- site_info$drain_area_va[1]
        
        all_data[[station]] <- daily_data
        cat("  Success! Downloaded", nrow(daily_data), "days of data\n")
      } else {
        cat("  No data available for this station\n")
      }
      
    }, error = function(e) {
      cat("  Error downloading station", station, ":", e$message, "\n")
    })
  }
  
  # Combine all station data
  if(length(all_data) > 0) {
    combined_data <- bind_rows(all_data, .id = "station_id")
    return(combined_data)
  } else {
    return(NULL)
  }
}

# Download all data
cat("Downloading USGS data for Clark Fork system...\n")
usgs_data <- download_usgs_data(all_stations, start_date, end_date, parameters)

# Add station metadata
usgs_data <- usgs_data %>%
  left_join(station_info, by = "station_id")

# ===============================================================
# DATA SUMMARY AND QUALITY CHECK (FIXED)
# ===============================================================

# Check what columns actually exist in the data
cat("\n=== AVAILABLE COLUMNS IN DATASET ===\n")
print(names(usgs_data))

# Function to safely check for parameter availability
safe_param_check <- function(data, param_column) {
  if(param_column %in% names(data)) {
    return(sum(!is.na(data[[param_column]])))
  } else {
    return(0)
  }
}

# Check what data we got (safely)
cat("\n=== DATA DOWNLOAD SUMMARY ===\n")
data_summary <- usgs_data %>%
  group_by(station_id, station_name, station_type) %>%
  summarise(
    n_days = n(),
    discharge_available = safe_param_check(pick(everything()), "X_00060_00003"),
    temp_available = safe_param_check(pick(everything()), "X_00010_00001"),
    conductance_available = safe_param_check(pick(everything()), "X_00095_00003"),
    ph_available = safe_param_check(pick(everything()), "X_00400_00003"),
    date_start = min(Date),
    date_end = max(Date),
    .groups = 'drop'
  )

print(data_summary)

# Alternative approach - check each parameter manually
cat("\n=== PARAMETER AVAILABILITY CHECK ===\n")
available_params <- names(usgs_data)[grepl("^X_", names(usgs_data))]
cat("Available parameter columns:\n")
print(available_params)

# Manual summary for available parameters only
for(param in available_params) {
  param_summary <- usgs_data %>%
    group_by(station_id, station_type) %>%
    summarise(
      available_days = sum(!is.na(.data[[param]])),
      .groups = 'drop'
    )
  
  cat("\n", param, "availability:\n")
  print(param_summary)
}

# ===============================================================
# CALCULATE FLOWS FOR YOUR SAMPLING SITES (UPDATED)
# ===============================================================

# Function to safely extract available parameters
calculate_site_flows_safe <- function(usgs_data) {
  
  # Check what discharge column exists
  discharge_col <- names(usgs_data)[grepl("X_00060", names(usgs_data))][1]
  temp_col <- names(usgs_data)[grepl("X_00010", names(usgs_data))][1]
  
  if(is.na(discharge_col)) {
    cat("No discharge data found!\n")
    return(NULL)
  }
  
  cat("Using discharge column:", discharge_col, "\n")
  if(!is.na(temp_col)) {
    cat("Using temperature column:", temp_col, "\n")
  } else {
    cat("No temperature data available\n")
  }
  
  # Extract main discharge data
  discharge_data <- usgs_data %>%
    filter(!is.na(.data[[discharge_col]])) %>%
    select(station_id, Date, station_type, sampling_site, 
           discharge_cfs = all_of(discharge_col))
  
  # Add temperature if available
  if(!is.na(temp_col)) {
    discharge_data <- discharge_data %>%
      left_join(
        usgs_data %>% select(station_id, Date, temp_c = all_of(temp_col)),
        by = c("station_id", "Date")
      )
  }
  
  # Get flows for your direct sampling sites
  site_flows <- discharge_data %>%
    filter(sampling_site %in% c("WS", "DL", "GR", "GC")) %>%
    arrange(Date, station_id)
  
  # Calculate tributary contributions
  tributary_flows <- discharge_data %>%
    filter(station_type == "tributary") %>%
    select(Date, station_id, sampling_site, discharge_cfs) %>%
    rename(tributary_flow = discharge_cfs, tributary = sampling_site)
  
  return(list(site_flows = site_flows, tributary_flows = tributary_flows))
}

# Calculate flows safely
flow_data <- calculate_site_flows_safe(usgs_data)

if(!is.null(flow_data)) {
  cat("\n=== FLOW DATA FOR YOUR SAMPLING SITES ===\n")
  print(head(flow_data$site_flows))
  
  cat("\n=== TRIBUTARY FLOW DATA ===\n") 
  print(head(flow_data$tributary_flows))
} else {
  cat("Could not extract flow data - check station availability\n")
}

# ===============================================================
# VISUALIZE DISCHARGE PATTERNS (UPDATED)
# ===============================================================

# Get the actual discharge column name
discharge_col <- names(usgs_data)[grepl("X_00060", names(usgs_data))][1]

if(!is.na(discharge_col)) {
  # Plot discharge over time for main stem sites
  main_stem_data <- usgs_data %>%
    filter(station_type == "main_stem" & !is.na(.data[[discharge_col]]))
  
  if(nrow(main_stem_data) > 0) {
    # Create discharge plot
    par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))
    
    # Plot main stem discharge
    plot(as.Date(main_stem_data$Date), main_stem_data[[discharge_col]],
         type = "n", xlab = "Date", ylab = "Discharge (cfs)",
         main = "Clark Fork Main Stem Discharge")
    
    # Add lines for each station
    colors <- rainbow(length(unique(main_stem_data$station_id)))
    stations <- unique(main_stem_data$station_id)
    
    for(i in 1:length(stations)) {
      station_data <- main_stem_data[main_stem_data$station_id == stations[i], ]
      if(nrow(station_data) > 0) {
        lines(as.Date(station_data$Date), station_data[[discharge_col]], 
              col = colors[i], lwd = 2)
      }
    }
    
    # Create legend with station info
    legend_text <- paste(stations, "-", 
                         main_stem_data$sampling_site[match(stations, main_stem_data$station_id)])
    legend("topright", legend = legend_text, col = colors, lwd = 2, cex = 0.7)
    
    cat("\n=== PLOTS CREATED ===\n")
    cat("Discharge plot shows data for", length(stations), "main stem stations\n")
  } else {
    cat("No main stem discharge data available for plotting\n")
  }
} else {
  cat("No discharge data columns found in dataset\n")
}

cat("\n=== DIAGNOSIS COMPLETE ===\n")
cat("Check the parameter availability above to see what data is actually available\n")

write.csv(usgs_data, "2_incremental/AI_Subsets/usgs_data.csv")
