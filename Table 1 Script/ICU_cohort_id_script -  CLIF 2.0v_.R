
"
Script Name: ICU_cohort_id_script.R
Purpose: Identify the cohort for CLIF concept paper POCs and generate table 1
Authors: Vaishvik Chaudhari, Kaveri Chhikara, Rachel Baccile
Inputs: 1. CLIF-2.0 tables - patient, hospitalization, ADT, vitals, labs, 
                          respiratory support, medication admin continuous, 
                          patient assessments

           required_meds        = norepinephrine, epinephrine, phenylephrine, vasopressin, 
                                   dopamine, angiotensin, dobutamine, milrinone
           required_vitals      = weight_kg, sbp, dbp, map, spo2
           required_labs        = creatinine, bilirubin_total, po2_arterial, platelet_count
           required_assessments = gcs_total
        2. config.json elements
Outputs: 1. icu_data.csv (intermediate dataset for further analysis); 
         2. table1_meds_<site>.csv; 
         3. table1_peep_fio2_<site>.csv; 
         4. table1_mode_category_<site>.csv; 
         5. table1_sofa_<site>.csv; 
         6. table1_<site>.csv
         7. histograms of vent settings and SOFA components
"

########################## SETUP ###############################################
packages <- c("jsonlite", "duckdb", "lubridate", "data.table",
              "tidyverse", "dplyr","table1",'rvest', "readr", "ggplot2",
              "arrow", "fst", "lightgbm", "caret", "Metrics", "patchwork",
              "ROCR", "pROC", "collapse")

install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

sapply(packages, install_if_missing)

con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = ":memory:")

load_config <- function() {
  json_path <- file.path( "config.json")
  
  if (file.exists(json_path)) {
    config <- fromJSON(json_path)
    print("Loaded configuration from config.json")
  } else {
    stop("Configuration file not found. Please create config.json based on the config_template.")
  }
  
  return(config)
}

# Load the configuration
config <- load_config()

tables_location <- config$clif2_path
site <-config$site
file_type <- paste0(".", config$filetype)

# Check if the output directory exists; if not, create it
if (!dir.exists("output")) {
  dir.create("output")
}

# Check if the output directory exists; if not, create it
if (!dir.exists("output/graphs")) {
  dir.create("output/graphs")
}

############################ FUNCTIONS  ########################################
read_data <- function(file_path) {
  if (grepl("\\.csv$", file_path)) {
    return(read.csv(file_path))
  } else if (grepl("\\.parquet$", file_path)) {
    return(arrow::read_parquet(file_path))
  } else if (grepl("\\.fst$", file_path)) {
    return(fst::read.fst(file_path))
  } else {
    stop("Unsupported file format")
  }
}

generate_histograms <- function(df, numeric_vars, bins = 30, 
                                color = "darkblue", 
                                output_path = "output/histograms.png",
                                max_cols = 3) {
  # Create individual histograms for each numeric variable
  plots <- lapply(numeric_vars, function(var) {
    ggplot(df, aes_string(x = var)) +
      geom_histogram(bins = bins, fill = color, color = "black", alpha = 0.7) +
      labs(
        title = paste("Histogram of", var),
        x = var,
        y = "Frequency"
      ) +
      theme_minimal()
  })
  # Arrange the plots in a grid
  combined_plot <- wrap_plots(plots, ncol = max_cols)
  # Save the combined plot
  ggsave(
    filename = output_path,
    plot = combined_plot,
    width = min(50, max_cols * 5),   # Adjust width based on columns
    height = min(50, ceiling(length(numeric_vars) / max_cols) * 5),  # Adjust height based on rows
    limitsize = FALSE               # Override the size limit
  )
  message("Histograms saved to: ", output_path)
}


get_conversion_factor <- function(med_category, med_dose_unit, weight_kg) {
  # 1) Retrieve medication info (if missing, return NA)
  med_info <- med_unit_info[[med_category]]
  if (is.null(med_info)) {
    return(NA_real_)
  }
  # 2) Convert the incoming unit to lowercase
  med_dose_unit <- tolower(med_dose_unit)
  # 3) Check if it's an acceptable unit; if not, return NA
  if (!(med_dose_unit %in% med_info$acceptable_units)) {
    return(NA_real_)
  }
  # 4) Determine conversion factor for each med+unit
  factor <- NA_real_
  # Group 1: norepinephrine, epinephrine, phenylephrine, dopamine, metaraminol, dobutamine
  #   required_unit: "mcg/kg/min"
  if (med_category %in% c("norepinephrine", "epinephrine", "phenylephrine",
                          "dopamine", "milrinone", "dobutamine")) {
    if (med_dose_unit == "mcg/kg/min") {
      factor <- 1
    } else if (med_dose_unit == "mcg/kg/hr") {
      factor <- 1 / 60
    } else if (med_dose_unit == "mg/kg/hr") {
      factor <- 1000 / 60
    } else if (med_dose_unit == "mcg/min") {
      factor <- 1 / weight_kg
    } else if (med_dose_unit == "mg/hr") {
      factor <- (1000 / 60) / weight_kg
    } else {
      return(NA_real_)
    }
    # Group 2: angiotensin
    #   required_unit: "mcg/kg/min"
  } else if (med_category == "angiotensin") {
    
    if (med_dose_unit == "ng/kg/min") {
      factor <- 1 / 1000
    } else if (med_dose_unit == "ng/kg/hr") {
      factor <- 1 / 1000 / 60
    } else {
      return(NA_real_)
    }
    # Group 3: vasopressin
    #   required_unit: "units/min"
  } else if (med_category == "vasopressin") {
    
    if (med_dose_unit == "units/min") {
      factor <- 1
    } else if (med_dose_unit == "units/hr") {
      factor <- 1 / 60
    } else if (med_dose_unit == "milliunits/min") {
      factor <- 1 / 1000
    } else if (med_dose_unit == "milliunits/hr") {
      factor <- 1 / 1000 / 60
    } else {
      return(NA_real_)
    }
    # If none of the above, return NA
  } else {
    return(NA_real_)
  }
  return(factor)
}


calc_pao2 <- function(s) {
  s <- s / 100
  a <- (11700) / ((1 / s) - 1)
  b <- sqrt((50^3) + (a^2))
  pao2 <- ((b + a)^(1/3)) - ((b - a)^(1/3))
  return(pao2)
}
############################     LOAD DATA         ############################
# Read data using the function and assign to variables
location <- read_data(paste0(tables_location, "/clif_adt", 
                             file_type))
encounter <- read_data(paste0(tables_location, "/clif_hospitalization", 
                              file_type))
demog <- read_data(paste0(tables_location, "/clif_patient", 
                          file_type))
ventilator <- read_data(paste0(tables_location, "/clif_respiratory_support", 
                               file_type))

############################ COHORT IDENTIFICATION ############################
# Rename columns in the location data frame
location <- location %>%
  rename(encounter_id = hospitalization_id)

# Rename columns in the encounter data frame
encounter <- encounter %>%
  rename(encounter_id = hospitalization_id)

# Rename columns in the demog data frame
demog <- demog %>%
  rename(
    race = race_category,
    ethnicity = ethnicity_category,
    sex = sex_category
  )

ventilator <- ventilator %>%
  rename(
    encounter_id = hospitalization_id
  )

# First join operation
join <- location %>%
  select(encounter_id, location_category, in_dttm, out_dttm)

# Second join operation to get 'icu_data'
icu_data <- join %>%
  left_join(encounter %>% select(patient_id, encounter_id, age_at_admission, 
                                 discharge_category,admission_dttm), 
            by = "encounter_id") %>%
  mutate(
    admission_dttm = ymd_hms(admission_dttm), # Convert to POSIXct, adjust the function as per your date format
    in_dttm = ymd_hms(in_dttm), # Convert to POSIXct, adjust the function as per your date format
    out_dttm = ymd_hms(out_dttm)
  )

# Filter rows where location is ICU and in_dttm is within 48 hours of admission_dttm
icu_data <- icu_data %>%
  mutate(location_category = ifelse(location_category == "procedural", "OR", location_category)) %>%
  mutate(location_category = toupper(location_category))

icu_48hr_check <- icu_data %>%
  filter(location_category == "ICU",
         in_dttm >= admission_dttm,
         in_dttm <= admission_dttm + lubridate::hours(48),
         lubridate::year(admission_dttm) >= 2020,
         lubridate::year(admission_dttm) <= 2021,
         age_at_admission >= 18,
         !is.na(age_at_admission)) %>%
  distinct(encounter_id) %>%
  pull(encounter_id)
  
# Filter icu_data to only include rows with encounter_ids in icu_48hr_check and within 72 hours of admission
icu_data <- icu_data %>%
  filter(encounter_id %in% icu_48hr_check,
         in_dttm <= admission_dttm + hours(72)) %>%
  arrange(in_dttm) %>%
  mutate(RANK = rank(in_dttm, ties.method = "first")) %>%
  arrange(encounter_id, in_dttm) %>%
  group_by(encounter_id) %>%
  mutate(RANK = rank(in_dttm, ties.method = "first"))

# Compute minimum rank for ICU locations
min_icu <- icu_data %>%
  filter(location_category == "ICU") %>%
  group_by(encounter_id) %>%
  summarize(min_icu = min(RANK))

# Merge the minimum ICU rank back into the original dataset
icu_data <- icu_data %>%
  left_join(min_icu, by = "encounter_id")

# Filter based on rank being at least the minimum ICU rank
icu_data <- icu_data %>%
  filter(RANK >= min_icu) %>%
  arrange(in_dttm)

# Change 'OR' to 'ICU' in location_category
icu_data <- icu_data %>%
  mutate(location_category = ifelse(location_category == "OR", "ICU", location_category))

# Create a new group_id based on changes in location_category
icu_data <- icu_data %>%
  group_by(encounter_id) %>%
  mutate(group_id = cumsum(location_category != lag(location_category, 
                                                    default = first(location_category)))) %>%
  ungroup()

icu_data <- icu_data %>%
  group_by(patient_id , encounter_id, location_category, group_id) %>%
  summarize(
    min_in_dttm = min(in_dttm),
    max_out_dttm = max(out_dttm),
    admission_dttm = first(admission_dttm),
    age = first(age_at_admission),
    dispo = first(discharge_category),
    .groups = 'drop'
  )

# Compute minimum group_id for each encounter_id where location_category is 'ICU'
min_icu <- icu_data %>%
  filter(location_category == "ICU") %>%
  group_by(encounter_id) %>%
  summarize(min_icu = min(group_id), .groups = 'drop')

# Merge the minimum ICU group_id back into the original dataset
icu_data <- left_join(icu_data, min_icu, by = "encounter_id")

# Filter based on group_id matching min_icu and duration condition
icu_data <- icu_data %>%
  filter(min_icu == group_id,
         interval(min_in_dttm, max_out_dttm) >= dhours(24)) %>%
  arrange(min_in_dttm)

  # Add 24 hours to the 'min_in_dttm' column
icu_data <- icu_data %>%
  mutate(after_24hr = min_in_dttm + hours(24))

# Select specific columns
icu_data <- icu_data %>%
  select(patient_id, encounter_id, min_in_dttm, after_24hr,max_out_dttm, 
         age, dispo)

# Merge with demographic data and select specific columns
icu_data <- icu_data %>%
  left_join(demog, by = "patient_id") %>%
  select(encounter_id, min_in_dttm, after_24hr,max_out_dttm, age, dispo, 
         sex, ethnicity, race)

ventilator_imv <- ventilator %>%
  filter(device_category =="IMV")%>%
  select(encounter_id) %>% distinct() %>% deframe()


# Remove rows with missing 'sex' and create new variables
icu_data <- icu_data %>%
  filter(!is.na(sex)) %>%
  mutate(
    isfemale = as.integer(tolower(sex) == "female"),
    Mortality  = as.integer(grepl("dead|expired|death|died", dispo, ignore.case = TRUE)),
    site = site,
    Ventilator = ifelse(encounter_id %in% ventilator_imv, 1, 0)
  )

# Define race and ethnicity mappings using case_when
icu_data <- icu_data %>%
  mutate(
    race = case_when(
      race == "White" ~ "White",
      race == "Black or African American" ~ "Black",
      race == "Black or African-American" ~ "Black",
      race == "Asian" ~ "Asian",
      race %in% c("Other", "Unknown", "Did Not Encounter", "Refusal", 
                  "American Indian or Alaska Native", 
                  "Native Hawaiian or Other Pacific Islander") ~ "Others",
      TRUE ~ "Others"  # Default case for NA and any other unexpected values
    ),
    ethnicity = case_when(
      ethnicity == "Hispanic" ~ "Hispanic or Latino",
      ethnicity == "Hispanic or Latino" ~ "Hispanic or Latino",
      TRUE ~ "Not Hispanic"  # Default case for NA and any other unexpected values
    )
  ) %>%
  mutate(sex = as.factor(sex),
         race = as.factor(race),
         ethnicity = as.factor(ethnicity),
         Mortality = as.factor(Mortality),
         Ventilator = as.factor(Ventilator))

# Calculate the difference in hours
icu_data$ICU_stay_hrs <- as.numeric(difftime(icu_data$max_out_dttm, 
                                             icu_data$min_in_dttm, 
                                             units = "secs")) / 3600

write.csv(icu_data, paste0( "output/ICU_cohort", '.csv'), row.names = FALSE)

rm(location, encounter)
gc()
############################VASOPRESSORS and VITALS############################
meds <- read_data(paste0(tables_location, "/clif_medication_admin_continuous", 
                         file_type)) 
vitals <- read_data(paste0(tables_location, "/clif_vitals", 
                           file_type))

required_meds <- c("norepinephrine", "epinephrine", "phenylephrine",
                   "vasopressin", "dopamine", "angiotensin", 
                   "dobutamine", "milrinone")
required_vitals <- c("weight_kg", "sbp", "dbp", "map")

cohort_ids <- icu_data |> 
  select(encounter_id) |> 
  distinct()

total_encounters <- nrow(cohort_ids)

vitals_weight_dt <- vitals |>
  rename(encounter_id = hospitalization_id) |>
  filter(encounter_id %in% cohort_ids$encounter_id) |>
  mutate(vital_value = as.numeric(vital_value)) |> 
  filter(vital_category == "weight_kg") |>
  select(encounter_id, recorded_dttm, weight_kg = vital_value) |>
  left_join(
    icu_data %>%
      select(encounter_id, min_in_dttm, after_24hr),
    by = "encounter_id") |>
  # filter to the first 24 hours
  mutate(recorded_dttm = ymd_hms(recorded_dttm)) |>  # Convert to POSIXct
  filter(
    recorded_dttm >= min_in_dttm,
    recorded_dttm <= after_24hr) |>
  # optionally keep only columns you need
  select(encounter_id, recorded_dttm, weight_kg)

meds_dt <- meds |>
  rename(encounter_id = hospitalization_id) |>
  filter(encounter_id %in% cohort_ids$encounter_id) |>
  filter(med_category %in% required_meds) |>
  select(encounter_id, admin_dttm, med_category, med_dose, med_dose_unit) |>
  # join to icu_data to get min_in_dttm and after_24hr
  left_join(
    icu_data %>%
      select(encounter_id, min_in_dttm, after_24hr),
    by = "encounter_id") |>
  # keep only rows in the first 24 hours of ICU
  mutate(admin_dttm = ymd_hms(admin_dttm)) |>
  filter(
    admin_dttm >= min_in_dttm,
    admin_dttm <= after_24hr) |>
  select(encounter_id, admin_dttm, med_category, med_dose, med_dose_unit)

rm(meds)
gc()

# Convert both to data.table
setDT(vitals_weight_dt)
setDT(meds_dt)
setkey(vitals_weight_dt, encounter_id, recorded_dttm)
setkey(meds_dt, encounter_id, admin_dttm)

# Perform rolling join
# forward picks the most recent weight at or before admin_time so 
# we don't jump forward in time to a future weight.
meds_with_weights_dt <- meds_dt[
  vitals_weight_dt,
  on   = .(encounter_id, admin_dttm = recorded_dttm),
  roll = Inf 
]

# remove rows with no med dose
meds_with_weights_dt <- meds_with_weights_dt[
  !is.na(med_dose)
]
# Convert med_dose_unit to lowercase
meds_with_weights_dt[, med_dose_unit := tolower(med_dose_unit)]

# Define medications and their unit conversion information
med_unit_info <- list(
  norepinephrine = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mg/kg/hr", "mcg/min", "mg/hr")
  ),
  epinephrine = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mg/kg/hr", "mcg/min", "mg/hr")
  ),
  phenylephrine = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mg/kg/hr", "mcg/min", "mg/hr")
  ),
  vasopressin = list(
    required_unit = "units/min",
    acceptable_units = c("units/min", "units/hr", "milliunits/min", "milliunits/hr")
  ),
  dopamine = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mg/kg/hr", "mcg/min", "mg/hr")
  ),
  angiotensin = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("ng/kg/min", "ng/kg/hr")
  ),
  dobutamine = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mg/kg/hr", "mcg/min", "mg/hr")
  ),
  milrinone = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mg/kg/hr", "mcg/min", "mg/hr")
  )
)

# Apply conversion logic
meds_with_weights_dt[, conversion_factor := mapply(
  get_conversion_factor,
  med_category,
  med_dose_unit,
  weight_kg
)]

meds_with_weights_dt[, med_dose_converted := 
                       fifelse(is.na(conversion_factor),
                               NA_real_,
                               round(med_dose * conversion_factor, 3)  # Round to 3 decimal places
                       )
]

## excluding zeroes- zeores correspond to when the med was paused or stopped
## not releavnt for calculating the median of the med dose 
med_summaries <- as_tibble(meds_with_weights_dt) |> 
  group_by(encounter_id, med_category) |> 
  summarize(
    median_dose = median(med_dose_converted, na.rm = TRUE),
    iqr_lower   = quantile(med_dose_converted, 0.25, na.rm = TRUE),
    iqr_upper   = quantile(med_dose_converted, 0.75, na.rm = TRUE),
    .groups     = "drop"
  )

table1_meds <- med_summaries |> 
  group_by(med_category) |> 
  summarize(
    overall_median_dose = median(median_dose[median_dose != 0], na.rm = TRUE),
    overall_iqr_lower   = quantile(median_dose[median_dose != 0], 0.25, na.rm = TRUE),
    overall_iqr_upper   = quantile(median_dose[median_dose != 0], 0.75, na.rm = TRUE),
    n_encounters        = n_distinct(encounter_id),
    pct_encounters      = 100 * n_encounters / total_encounters,
    .groups             = "drop"
  ) |> # the join ensures all required meds are in the final table
  right_join(
    tibble(med_category = required_meds), 
    by = "med_category"
  )

write.csv(table1_meds, paste0( "output/table1_meds_",site, '.csv'), 
          row.names = FALSE)

######################VENTILATOR SUMMARY########################################
# unique encounters on vent
enc_on_vent <- icu_data |> 
  filter(Ventilator == 1) |> 
  select(encounter_id) |> distinct()

vent_encounters <- nrow(enc_on_vent)

ventilator_filtered <- ventilator |> 
  filter(encounter_id %in% enc_on_vent$encounter_id) |> 
  left_join(
    icu_data %>%
      select(encounter_id, min_in_dttm, after_24hr),
    by = "encounter_id"
  ) |> 
  mutate(
    recorded_dttm = ymd_hms(recorded_dttm),
    fio2_set = as.numeric(fio2_set),
    peep_set = as.numeric(peep_set),
    # Convert fio2_set to decimal scale if its mean is greater than 1
    fio2_set = if (mean(fio2_set, na.rm = TRUE) > 1) fio2_set / 100 else fio2_set
  ) |> 
  # Apply thresholds: convert outliers to NA
  mutate(
    peep_set = ifelse(peep_set >= 0 & peep_set <= 30, peep_set, NA_real_),
    fio2_set = ifelse(fio2_set >= 0.21 & fio2_set <= 1, fio2_set, NA_real_)
  ) |> 
  select(
    encounter_id, recorded_dttm, min_in_dttm, after_24hr,
    device_category, mode_category, 
    fio2_set, peep_set
  )


numeric_vars <- c("fio2_set", "peep_set")
generate_histograms(
  df = ventilator_filtered, 
  numeric_vars = numeric_vars, 
  bins = 30, 
  color = "blue", 
  output_path = paste0("output/graphs/histograms_fio2_peep_", site, ".png"),
)

# Summarize peep and fio2 during the first 24 hours of ICU admission
## Encounters not on vent during the first 24 hrs of ICU admission will be excluded
table1_peep_fio2 <- ventilator_filtered |> 
  filter(
    recorded_dttm >= min_in_dttm,
    recorded_dttm <= after_24hr) |>
  summarize(
    median_peep = median(peep_set, na.rm = TRUE),
    iqr_peep_lower = quantile(peep_set, 0.25, na.rm = TRUE),
    iqr_peep_upper = quantile(peep_set, 0.75, na.rm = TRUE),
    
    median_fio2 = median(fio2_set, na.rm = TRUE),
    iqr_fio2_lower = quantile(fio2_set, 0.25, na.rm = TRUE),
    iqr_fio2_upper = quantile(fio2_set, 0.75, na.rm = TRUE)
  )

write.csv(table1_peep_fio2, paste0( "output/table1_peep_fio2_",site, '.csv'), 
          row.names = FALSE)

# Summarise initial mode category for each encounter
# Note: without filling missing values, the number of encounters with info on 
# mode catgeory will be less than the number of encounters identified on IMV.
# Not filtering time for first 24 hrs because we have not applied the waterfall here

required_modes <- c('Assist Control-Volume Control', 
                    'Pressure Control', 
                    'Pressure-Regulated Volume Control', 
                    'SIMV', 'Pressure Support/CPAP', 'Volume Support', 
                    'Other')

mode_cat_filtered <- ventilator_filtered |>
  filter(!is.na(mode_category)& mode_category != "") |> 
  arrange(encounter_id, recorded_dttm) |>  
  group_by(encounter_id) |> 
  slice_head(n=1) |>  # Get the first row for each encounter
  ungroup() 

## MODE CATEGORY WITH VENT ENCOUNTERS
table1_mode_category <- mode_cat_filtered |> 
  group_by(mode_category) |>  
  summarize(
    n_encounters = n(),  # Count unique encounters with initial mode
    .groups = "drop"  
  ) |> 
  mutate(
    pct_encounters = (n_encounters / vent_encounters) * 100  # Calculate percentage
  ) |> 
  right_join(
    tibble(mode_category = required_modes), 
    by = "mode_category"
  )

no_mode_count <- vent_encounters - sum(table1_mode_category$n_encounters, na.rm = TRUE)
no_mode_pct <- (no_mode_count / vent_encounters) * 100

# Append the "No mode documented" row
table1_mode_category <- table1_mode_category |> 
  bind_rows(
    tibble(
      mode_category = "No mode documented",
      n_encounters = no_mode_count,
      pct_encounters = no_mode_pct
    )
  )

table1_mode_category <- table1_mode_category |> 
  bind_rows(
    tibble(
      mode_category = "Total",
      n_encounters = sum(table1_mode_category$n_encounters, na.rm = TRUE),
      pct_encounters = sum(table1_mode_category$pct_encounters, na.rm = TRUE),
    )
  )

## MODE CATEGORY WITH ALL ENCOUNTERS
table1_mode_category_overall <- mode_cat_filtered |> 
  group_by(mode_category) |>  
  summarize(
    n_encounters = n(),  # Count unique encounters with initial mode
    .groups = "drop"  
  ) |> 
  mutate(
    pct_encounters = (n_encounters / total_encounters) * 100  # Calculate percentage
  ) |> 
  right_join(
    tibble(mode_category = required_modes), 
    by = "mode_category"
  )

no_mode_count_overall <- total_encounters - sum(table1_mode_category_overall$n_encounters, na.rm = TRUE)
no_mode_pct_overall <- (no_mode_count_overall / total_encounters) * 100

# Append the "No mode documented" row
table1_mode_category_overall <- table1_mode_category_overall |> 
  bind_rows(
    tibble(
      mode_category = "No mode documented",
      n_encounters = no_mode_count_overall,
      pct_encounters = no_mode_pct_overall
    )
  ) 

table1_mode_category_overall <- table1_mode_category_overall |> 
  bind_rows(
    tibble(
      mode_category = "Total",
      n_encounters = sum(table1_mode_category_overall$n_encounters, na.rm = TRUE),
      pct_encounters = sum(table1_mode_category_overall$pct_encounters, na.rm = TRUE),
    )
  )

write.csv(table1_mode_category, paste0( "output/table1_mode_category_only_vent",site, '.csv'), 
          row.names = FALSE)
write.csv(table1_mode_category_overall, paste0( "output/table1_mode_category_all_encounters",site, '.csv'), 
          row.names = FALSE)

############################SOFA-97#############################################

# Start with Vitals to get MAP and SpO2
cohort_24hr <- icu_data |> 
  select(encounter_id, min_in_dttm, after_24hr) |> 
  distinct()

vitals_sofa_dt <- vitals |>
  rename(encounter_id = hospitalization_id) |>
  mutate(vital_value = as.numeric(vital_value)) |> 
  filter(encounter_id %in% cohort_24hr$encounter_id) |>
  filter(
    (vital_category == "spo2" & vital_value >= 60) |
      (vital_category == "map" & vital_value >= 30)) |>
  select(encounter_id, recorded_dttm, vital_category, vital_value)

# Filter to first 24 hours in ICU
vitals_icu <- cohort_24hr %>%
  left_join(vitals_sofa_dt) %>%
  filter((recorded_dttm >= min_in_dttm) & (recorded_dttm <=after_24hr)) %>%
  select(encounter_id, recorded_dttm, vital_category, vital_value)%>%
  group_by(encounter_id, recorded_dttm, vital_category) %>%
  summarise(worst_value = min(vital_value, na.rm = TRUE)) %>%
  pivot_wider(
    names_from = vital_category,
    values_from = worst_value)
rm(vitals_sofa_dt)

vitals_icu <- vitals_icu %>%
  mutate(pao2_imputed_min = calc_pao2(spo2)) %>%
  # Replace with NA if SpO2 >= 97
  mutate(pao2_imputed_min = ifelse(spo2>=97, NA, pao2_imputed_min)) %>%
  group_by(encounter_id) %>%
  summarise(
    pao2_imputed_min = if (all(is.na(pao2_imputed_min))) NA_real_ else min(pao2_imputed_min, na.rm = TRUE),
    spo2_min = if (all(is.na(spo2))) NA_real_ else min(spo2, na.rm = TRUE),
    map_min = if (all(is.na(map))) NA_real_ else min(map, na.rm = TRUE),
    .groups = "drop")

resp_support_dt <- ventilator |>
  filter(encounter_id %in% cohort_24hr$encounter_id) |>
  select(encounter_id, recorded_dttm, device_category, mode_category, peep_set, tidal_volume_set, 
         lpm_set, fio2_set) |>
  mutate(fio2_set = ifelse((fio2_set<0.21 | fio2_set>1), NA, fio2_set)) 
resp_support_dt <- as_tibble(resp_support_dt)
resp_support_dt$encounter_id <- as.character(resp_support_dt$encounter_id)

resp_support_icu <- cohort_24hr %>%
  left_join(resp_support_dt) %>%
  filter((recorded_dttm >= min_in_dttm) & (recorded_dttm <=after_24hr)) %>%
  # Try to fill in device based on other values
  mutate(device_category = case_when(
    !is.na(device_category) ~ device_category,
    mode_category %in% c("SIMV", "Pressure-regulated Volume Control", "Assist Control-Volume Control") ~ "Vent",
    is.na(device_category) & fio2_set == 0.21 & is.na(lpm_set) & is.na(peep_set) & is.na(tidal_volume_set) ~ "Room Air",
    is.na(device_category) & is.na(fio2_set) & lpm_set == 0 & is.na(peep_set) & is.na(tidal_volume_set) ~ "Room Air",
    is.na(device_category) & is.na(fio2_set) & lpm_set <= 20 & lpm_set > 0 & is.na(peep_set) & is.na(tidal_volume_set) ~ "Nasal Cannula",
    is.na(device_category) & is.na(fio2_set) & lpm_set > 20 & is.na(peep_set) & is.na(tidal_volume_set) ~ "High Flow NC",
    device_category == "Nasal Cannula" & is.na(fio2_set) & lpm_set > 20 ~ "High Flow NC",
    TRUE ~ device_category # Keep original device_category if no condition is met
  )) %>%
  # Try to fill in FiO2 based on other values
  mutate(fio2_combined = case_when(
    !is.na(fio2_set) ~ fio2_set,
    is.na(fio2_set) & device_category == "Room Air" ~ 0.21,
    is.na(fio2_set) & device_category == "Nasal Cannula" ~ (0.24 + (0.04 * lpm_set)),
    TRUE ~ NA_real_
  )) %>%
  select(encounter_id, recorded_dttm, device_category, fio2_combined) %>%
  # Rank devices to get highest level of support
  mutate(device_rank = case_when(
    device_category == 'Vent' ~ 1,
    device_category == 'NIPPV' ~ 2,
    device_category == 'CPAP' ~ 3,
    device_category == 'High Flow NC' ~ 4,
    device_category == 'Face Mask' ~ 5,
    device_category == 'Trach Collar' ~ 6,
    device_category == 'Nasal Cannula' ~ 7,
    device_category == 'Other' ~ 8,
    device_category == 'Room Air' ~ 9,
    TRUE ~ NA_integer_
  )) %>%
  group_by(encounter_id) %>%
  # Get worst FiO2 and device
  summarise(
    device_rank_min = if (all(is.na(device_rank))) NA_real_ else min(device_rank, na.rm = TRUE),
    fio2_max = if (all(is.na(fio2_combined))) NA_real_ else max(fio2_combined, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(resp_support_max = case_when(
    device_rank_min == 1 ~ 'Vent',
    device_rank_min == 2 ~ 'NIPPV',
    device_rank_min == 3 ~ 'CPAP',
    device_rank_min == 4 ~ 'High Flow NC',
    device_rank_min == 5 ~ 'Face Mask',
    device_rank_min == 6 ~ 'Trach Collar',
    device_rank_min == 7 ~ 'Nasal Cannula',
    device_rank_min == 8 ~ 'Other',
    device_rank_min == 9 ~ 'Room Air',
    TRUE ~ NA_character_
  )) %>%
  select(encounter_id, resp_support_max, fio2_max)

## Now Labs
labs <- read_data(paste0(tables_location, 
                         "/clif_labs", 
                         file_type))

labs_dt <- labs |>
  rename(encounter_id = hospitalization_id) |>
  filter(encounter_id %in% cohort_24hr$encounter_id) |>
  filter(
    (lab_category == "creatinine" & lab_value_numeric >= 0 & lab_value_numeric <= 20) |
      (lab_category == "bilirubin_total" & lab_value_numeric >= 0 & lab_value_numeric <= 80)|
      (lab_category == "po2_arterial" & lab_value_numeric >= 30 & lab_value_numeric <= 700)|
      (lab_category == "platelet_count" & lab_value_numeric >= 0 & lab_value_numeric <= 2000)) |>
  filter(!is.na(lab_value_numeric)) |>
  select(encounter_id, lab_order_dttm, lab_category, lab_value_numeric) 

rm(labs)
gc()

# Filter for labs taken within first 24 hours of being in the ICU
labs_icu <- cohort_24hr %>%
  left_join(labs_dt) %>%
  filter((lab_order_dttm >= min_in_dttm) & (lab_order_dttm <=after_24hr)) %>%
  distinct() %>%
  # Get unique values by time
  # Summarize the minimum and maximum values of each lab category
  group_by(encounter_id, min_in_dttm, after_24hr, lab_order_dttm, lab_category) %>% 
  reframe(
    worst_value = case_when(
      lab_category %in% c("creatinine", "bilirubin_total") ~ max(lab_value_numeric, na.rm = TRUE),
      lab_category %in% c("po2_arterial", "platelet_count") ~ min(lab_value_numeric, na.rm = TRUE),
    )) %>% 
  distinct() %>%
  # Pivot wider
  pivot_wider(
    names_from = lab_category,
    values_from = worst_value) %>%
  # Summarize the minimum and maximum values of each lab category for whole 24 hours, return NA if missing
  group_by(encounter_id, min_in_dttm, after_24hr) %>% 
  summarise(
    creatinine_max = if (all(is.na(creatinine))) NA_real_ else max(creatinine, na.rm = TRUE),
    bilirubin_total_max = if (all(is.na(bilirubin_total))) NA_real_ else max(bilirubin_total, na.rm = TRUE),
    po2_arterial_min = if (all(is.na(po2_arterial))) NA_real_ else min(po2_arterial, na.rm = TRUE),
    platelet_count_min = if (all(is.na(platelet_count))) NA_real_ else min(platelet_count, na.rm = TRUE),
    .groups = "drop")

## Now GCS Scores
scores <- read_data(paste0(tables_location, 
                           "/clif_patient_assessments", 
                           file_type))
scores_dt <- scores |>
  rename(encounter_id = hospitalization_id) |>
  filter(encounter_id %in% cohort_24hr$encounter_id) |>
  filter(assessment_category == "gcs_total") |> 
  filter(!is.na(numerical_value)) |>
  select(encounter_id, recorded_dttm, numerical_value)
rm(scores)
gc()

scores_dt$encounter_id <- as.character(scores_dt$encounter_id)
scores_dt$min_gcs_score <- as.numeric(scores_dt$numerical_value)

scores_icu <- cohort_24hr %>%
  left_join(scores_dt) %>%
  filter((recorded_dttm >= min_in_dttm) & (recorded_dttm <=after_24hr)) %>%
  distinct() %>%
  group_by(encounter_id, min_in_dttm, after_24hr) %>%
  summarise(min_gcs_score = min(numerical_value, na.rm = TRUE))

## Now meds
meds_icu <- cohort_24hr %>%
  left_join(meds_with_weights_dt) %>%
  filter((admin_dttm >= min_in_dttm) & (admin_dttm <=after_24hr)) %>%
  select(encounter_id, admin_dttm, med_category, med_dose_converted) %>%
  distinct() %>%
  group_by(encounter_id,  med_category) %>%
  summarise(worst_value = max(med_dose_converted, na.rm = TRUE)) %>%
  pivot_wider(
    names_from = med_category,
    values_from = worst_value)


# Join to Vitals, Respiratory Support, Labs, Meds, and Scores to calculate SOFA
icu_sofa_data <- cohort_24hr %>% 
  left_join(resp_support_icu) %>%
  left_join(vitals_icu) %>%
  left_join(labs_icu) %>%
  left_join(scores_icu) %>%
  left_join(meds_icu) %>%
  mutate(p_f = ifelse(!is.na(fio2_max) & !is.na(po2_arterial_min), po2_arterial_min / fio2_max, NA),
         p_f_imputed = ifelse(!is.na(fio2_max) & !is.na(pao2_imputed_min), pao2_imputed_min / fio2_max, NA),
         s_f = ifelse(!is.na(fio2_max) & !is.na(spo2_min), spo2_min / fio2_max, NA))


# Calculate SOFA
icu_sofa_data <- icu_sofa_data %>%
  mutate(
    sofa_cv_97 = case_when(
      dopamine > 15 | epinephrine > 0.1 | norepinephrine > 0.1 ~ 4,
      dopamine > 5 | (epinephrine <= 0.1 & epinephrine>0) | (norepinephrine <= 0.1 & norepinephrine>0)~ 3,
      (dopamine <= 5 & dopamine>0) | dobutamine>0 ~ 2,
      map_min < 70 ~ 1,
      TRUE ~ 0
    ),
    sofa_coag = case_when(
      platelet_count_min < 20 ~ 4,
      platelet_count_min < 50 ~ 3,
      platelet_count_min < 100 ~ 2,
      platelet_count_min < 150 ~ 1,
      TRUE ~ 0
    ),
    sofa_liver = case_when(
      bilirubin_total_max >= 12.0 ~ 4,
      bilirubin_total_max >= 6.0 ~ 3,
      bilirubin_total_max >= 2.0 ~ 2,
      bilirubin_total_max >= 1.2 ~ 1,
      TRUE ~ 0
    ),
    sofa_renal = case_when(
      creatinine_max >= 5.0 ~ 4,
      creatinine_max >= 3.5 & creatinine_max < 5.0 ~ 3,
      creatinine_max >= 2.0 & creatinine_max < 3.5 ~ 2,
      creatinine_max >= 1.2 & creatinine_max < 2.0 ~ 1,
      TRUE ~ 0
    ),
    sofa_resp_pf = case_when(
      p_f < 100 & (resp_support_max == "NIPPV" | resp_support_max == "CPAP" | resp_support_max == "Vent") ~ 4,
      p_f < 200 & (resp_support_max == "NIPPV"| resp_support_max == "CPAP" | resp_support_max == "Vent") ~ 3,
      p_f < 300 ~ 2,
      p_f < 400 ~ 1,
      TRUE ~ 0
    ),
    sofa_resp_pf_imp = case_when(
      p_f_imputed < 100 & (resp_support_max == "NIPPV" | resp_support_max == "CPAP" | resp_support_max == "Vent") ~ 4,
      p_f_imputed < 200 & (resp_support_max == "NIPPV" | resp_support_max == "CPAP" | resp_support_max == "Vent") ~ 3,
      p_f_imputed < 300 ~ 2,
      p_f_imputed < 400 ~ 1,
      TRUE ~ 0
    ),
    sofa_cns = case_when(
      min_gcs_score < 6 ~ 4,
      min_gcs_score >= 6 & min_gcs_score <= 9 ~ 3,
      min_gcs_score >= 10 & min_gcs_score <= 12 ~ 2,
      min_gcs_score >= 13 & min_gcs_score <= 14 ~ 1,
      TRUE ~ 0
    ))
icu_sofa_data <- icu_sofa_data %>% rowwise() %>%
  mutate(sofa_resp = max(sofa_resp_pf, sofa_resp_pf_imp, na.rm=T))

icu_sofa_data$sofa_97_24hr <- icu_sofa_data$sofa_cv_97 + 
  icu_sofa_data$sofa_coag + 
  icu_sofa_data$sofa_renal + 
  icu_sofa_data$sofa_liver + 
  icu_sofa_data$sofa_resp + 
  icu_sofa_data$sofa_cns

numeric_vars <- names(icu_sofa_data)[sapply(icu_sofa_data, is.numeric)]
generate_histograms(
  df = icu_sofa_data, 
  numeric_vars = numeric_vars, 
  bins = 30, 
  color = "blue", 
  output_path = paste0("output/graphs/histograms_sofa24_hrs_", site, ".png"),
  max_cols = 6
)

table1_sofa <- icu_sofa_data %>%
  select(encounter_id, sofa_97_24hr, sofa_cv_97, sofa_coag, sofa_renal, sofa_liver, sofa_resp, sofa_cns) %>%
  pivot_longer(
    cols = -encounter_id,
    names_to = "sofa_category",
    values_to = "sofa_value") %>%
  group_by(sofa_category) %>%
  summarize(
    overall_median = median(sofa_value, na.rm = TRUE),
    overall_iqr_lower   = quantile(sofa_value, 0.25, na.rm = TRUE),
    overall_iqr_upper   = quantile(sofa_value, 0.75, na.rm = TRUE),
    .groups             = "drop"
  )

write.csv(table1_sofa, paste0( "output/table1_sofa_",site, '.csv'), 
          row.names = FALSE)

############################TABLE ONE############################################
# HTML content (make sure your actual HTML string is correctly input here)
html_content <- table1(~ sex + age + race + ethnicity + Mortality + Ventilator + ICU_stay_hrs, data=icu_data)

# Use rvest to read the HTML table
table <- read_html(html_content) %>%
  html_table(fill = TRUE)

# The first element of the list should be your table
df <- table[[1]]

# Rename 'Overall(N=14598)' to 'fabc(N=14598)' using the site variable
names(df) <- gsub("Overall\\(N=(\\d+)\\)", paste0(site, ' ', "(N=\\1)"), names(df))
write.csv(df, paste0( "output/table1_",site, '.csv'), row.names = FALSE)
