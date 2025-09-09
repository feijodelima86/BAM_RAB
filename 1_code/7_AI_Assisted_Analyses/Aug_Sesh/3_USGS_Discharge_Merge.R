# ===============================================================
# MERGE USGS DISCHARGE DATA WITH SPECIATION_SPATIAL
# INCLUDING TRIBUTARY FLOW ADDITIONS
# ===============================================================

library(dplyr)

usgs_data <-read.csv("2_incremental/AI_Subsets/usgs_data.csv")

# Parse dates in speciation_spatial first
speciation_spatial$Date_parsed <- as.Date(speciation_spatial$Date, format = "%m/%d/%Y")

# Check what dates we have in speciation_spatial
print("Sampling dates in speciation_spatial:")
print(unique(speciation_spatial$Date_parsed))

# Extract main stem discharge data
main_stem_discharge <- usgs_data %>%
  filter(station_id %in% c("12323750", "12324200", "12324400", "12324680")) %>%
  select(Date, station_id, discharge_cfs = X_00060_00003) %>%
  filter(!is.na(discharge_cfs)) %>%
  mutate(
    Date = as.Date(Date),
    discharge_m3s = discharge_cfs * 0.0283168,
    site_name = case_when(
      station_id == "12323750" ~ "WS",  
      station_id == "12324200" ~ "DL",    
      station_id == "12324400" ~ "GR",      
      station_id == "12324680" ~ "GC",    
      TRUE ~ NA_character_
    )
  )

# Extract tributary discharge data
tributary_discharge <- usgs_data %>%
  filter(station_id %in% c("12324590", "12331500", "12334510")) %>%
  select(Date, station_id, discharge_cfs = X_00060_00003) %>%
  filter(!is.na(discharge_cfs)) %>%
  mutate(
    Date = as.Date(Date),
    discharge_m3s = discharge_cfs * 0.0283168,
    tributary_name = case_when(
      station_id == "12324590" ~ "Little_Blackfoot",  
      station_id == "12331500" ~ "Flint_Creek",       
      station_id == "12334510" ~ "Rock_Creek",        
      TRUE ~ NA_character_
    )
  )

# Check date ranges
print("Date range in main stem discharge data:")
print(range(main_stem_discharge$Date))
print("Date range in tributary discharge data:")
print(range(tributary_discharge$Date))

# Function to get discharge with flexible date matching
get_discharge_value <- function(discharge_df, target_date, station_name, col_name = "discharge_m3s") {
  # Try exact match first
  exact_match <- discharge_df %>%
    filter(Date == target_date & site_name == station_name)
  
  if(nrow(exact_match) > 0) {
    return(exact_match[[col_name]][1])
  }
  
  # If no exact match, find closest date within 2 days
  nearby <- discharge_df %>%
    filter(site_name == station_name) %>%
    mutate(date_diff = abs(as.numeric(Date - target_date))) %>%
    filter(date_diff <= 2) %>%
    arrange(date_diff)
  
  if(nrow(nearby) > 0) {
    return(nearby[[col_name]][1])
  }
  
  return(NA_real_)
}

# Function to get tributary discharge
get_trib_discharge <- function(trib_df, target_date, trib_name) {
  # Try exact match first
  exact_match <- trib_df %>%
    filter(Date == target_date & tributary_name == trib_name)
  
  if(nrow(exact_match) > 0) {
    return(exact_match$discharge_m3s[1])
  }
  
  # If no exact match, find closest date within 2 days
  nearby <- trib_df %>%
    filter(tributary_name == trib_name) %>%
    mutate(date_diff = abs(as.numeric(Date - target_date))) %>%
    filter(date_diff <= 2) %>%
    arrange(date_diff)
  
  if(nrow(nearby) > 0) {
    return(nearby$discharge_m3s[1])
  }
  
  return(0)  # Return 0 if no tributary data (assume no flow rather than NA)
}

# Calculate discharge for each site
speciation_spatial$discharge_m3s <- NA_real_

for(i in 1:nrow(speciation_spatial)) {
  site <- speciation_spatial$Site[i]
  date <- speciation_spatial$Date_parsed[i]
  
  if(site == "WS") {
    # WS (Site 1): Use Warm Springs flow directly
    speciation_spatial$discharge_m3s[i] <- get_discharge_value(main_stem_discharge, date, "WS")
    
  } else if(site == "DL") {
    # DL (Site 6): Use Deer Lodge flow
    speciation_spatial$discharge_m3s[i] <- get_discharge_value(main_stem_discharge, date, "DL")
    
  } else if(site == "GR") {
    # GR (Site 8): Use Garrison flow (before Little Blackfoot)
    speciation_spatial$discharge_m3s[i] <- get_discharge_value(main_stem_discharge, date, "GR")
    
  } else if(site == "GC") {
    # GC (Site 9): Use Gold Creek flow (before Flint Creek)
    speciation_spatial$discharge_m3s[i] <- get_discharge_value(main_stem_discharge, date, "GC")
    
  } else if(site == "BG") {
    # BG (Site 10 - Bear Gulch): Gold Creek + Flint Creek
    gc_flow <- get_discharge_value(main_stem_discharge, date, "GC")
    flint_flow <- get_trib_discharge(tributary_discharge, date, "Flint_Creek")
    
    if(!is.na(gc_flow)) {
      speciation_spatial$discharge_m3s[i] <- gc_flow + flint_flow
    }
    
  } else if(site == "BN") {
    # BN (Site 11 - Bonita): Gold Creek + Flint Creek + Rock Creek
    gc_flow <- get_discharge_value(main_stem_discharge, date, "GC")
    flint_flow <- get_trib_discharge(tributary_discharge, date, "Flint_Creek")
    rock_flow <- get_trib_discharge(tributary_discharge, date, "Rock_Creek")
    
    if(!is.na(gc_flow)) {
      speciation_spatial$discharge_m3s[i] <- gc_flow + flint_flow + rock_flow
    }
    
  } else {
    # For any other site names (sites 2-5, 7, 12-13), assign based on location
    # You may need to add these if they exist in your data
    print(paste("Unhandled site:", site))
  }
}

# Report results
na_count <- sum(is.na(speciation_spatial$discharge_m3s))
successful <- sum(!is.na(speciation_spatial$discharge_m3s))
print(paste("Successfully matched discharge for", successful, "rows"))
print(paste("Failed to match discharge for", na_count, "rows"))

if(na_count > 0) {
  print("Sites with missing discharge:")
  missing_sites <- speciation_spatial %>%
    filter(is.na(discharge_m3s)) %>%
    select(Site, Date) %>%
    unique()
  print(missing_sites)
}

# Clean up - remove temporary column
speciation_spatial <- speciation_spatial %>%
  select(-Date_parsed)