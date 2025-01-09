packages <- c("jsonlite", "duckdb", "lubridate", "data.table",
              "tidyverse", "dplyr","table1",'rvest', "readr", "arrow", "fst", "lightgbm", "caret", "Metrics", "ROCR", "pROC")

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


# Read data using the function and assign to variables
location <- read_data(paste0(tables_location, "/clif_adt", file_type))
encounter <- read_data(paste0(tables_location, "/clif_hospitalization", file_type))
demog <- read_data(paste0(tables_location, "/clif_patient", file_type))
ventilator <- read_data(paste0(tables_location, "/clif_respiratory_support", file_type))

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
  left_join(encounter %>% select(patient_id, encounter_id, age_at_admission, discharge_category,admission_dttm), by = "encounter_id") %>%
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
  mutate(group_id = cumsum(location_category != lag(location_category, default = first(location_category)))) %>%
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
  select(patient_id, encounter_id, min_in_dttm, after_24hr,max_out_dttm, age, dispo)

# Merge with demographic data and select specific columns
icu_data <- icu_data %>%
  left_join(demog, by = "patient_id") %>%
  select(encounter_id, min_in_dttm, after_24hr,max_out_dttm, age, dispo, sex, ethnicity, race)

ventilator <- ventilator %>%
  filter(device_category =="IMV")%>%
  select(encounter_id) %>% distinct() %>% deframe()


# Remove rows with missing 'sex' and create new variables
icu_data <- icu_data %>%
  filter(!is.na(sex)) %>%
  mutate(
    isfemale = as.integer(tolower(sex) == "female"),
    Mortality  = as.integer(grepl("dead|expired|death|died", dispo, ignore.case = TRUE)),
    site = site,
    Ventilator = ifelse(encounter_id %in% ventilator, 1, 0)
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
icu_data$ICU_stay_hrs <- as.numeric(difftime(icu_data$max_out_dttm, icu_data$min_in_dttm, units = "secs")) / 3600

write.csv(icu_data, paste0( "output/ICU_cohort", '.csv'), row.names = FALSE)


############################VASOPRESSORS and VITALS############################

required_meds <- c("norepinephrine", "epinephrine", "phenylephrine",
                   "vasopressin", "dopamine", "angiotensin", "dobutamine")
required_vitals <- c("weight_kg", "sbp", "dbp", "map")

meds <- arrow::open_dataset(paste0(tables_location, 
                                   "/clif_medication_admin_continuous", 
                                   file_type)) 
vitals <- arrow::open_dataset(paste0(tables_location, 
                                     "/clif_vitals", 
                                     file_type))

cohort_ids <- icu_data |> 
  select(encounter_id) |> 
  distinct()

total_encounters <- nrow(cohort_ids)

vitals_weight_dt <- vitals |>
  rename(encounter_id = hospitalization_id) |>
  filter(encounter_id %in% cohort_ids$encounter_id) |>
  filter(vital_category == "weight_kg") |>
  select(encounter_id, recorded_dttm, weight_kg = vital_value) |>
  collect()

meds_dt <- meds |>
  rename(encounter_id = hospitalization_id) |>
  filter(encounter_id %in% cohort_ids$encounter_id) |>
  filter(med_category %in% required_meds) |>
  select(encounter_id, admin_dttm, med_category, med_dose, med_dose_unit) |>
  collect()


# Convert both to data.table
setDT(vitals_weight_dt)
setDT(meds_dt)
setkey(vitals_weight_dt, encounter_id, recorded_dttm)
setkey(meds_dt, encounter_id, admin_dttm)

# Perform rolling join
#forward picks the most recent weight at or before admin_time so we don't jump forward in time to a future weight.
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
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mg/kg/hr", "mcg/min", "mg/hr"),
    conversion_factors = c("mcg/kg/min" = 1, "mcg/kg/hr" = 1/60, "mg/kg/hr" = 1000/60, "mcg/min" = 1, "mg/hr" = 1000/60)
  ),
  epinephrine = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mg/kg/hr", "mcg/min", "mg/hr"),
    conversion_factors = c("mcg/kg/min" = 1, "mcg/kg/hr" = 1/60, "mg/kg/hr" = 1000/60, "mcg/min" = 1, "mg/hr" = 1000/60)
  ),
  phenylephrine = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mg/kg/hr", "mcg/min", "mg/hr"),
    conversion_factors = c("mcg/kg/min" = 1, "mcg/kg/hr" = 1/60, "mg/kg/hr" = 1000/60, "mcg/min" = 1, "mg/hr" = 1000/60)
  ),
  vasopressin = list(
    required_unit = "units/min",
    acceptable_units = c("units/min", "units/hr", "milliunits/min", "milliunits/hr"),
    conversion_factors = c("units/min" = 1, "units/hr" = 1/60, "milliunits/min" = 1/1000, "milliunits/hr" = 1/1000/60)
  ),
  dopamine = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mg/kg/hr", "mcg/min", "mg/hr"),
    conversion_factors = c("mcg/kg/min" = 1, "mcg/kg/hr" = 1/60, "mg/kg/hr" = 1000/60, "mcg/min" = 1, "mg/hr" = 1000/60)
  ),
  angiotensin = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("ng/kg/min", "ng/kg/hr"),
    conversion_factors = c("ng/kg/min" = 1/1000, "ng/kg/hr" = 1/1000/60)
  ),
  dobutamine = list(
    required_unit = "mcg/kg/min",
    acceptable_units = c("mcg/kg/min", "mcg/kg/hr", "mg/kg/hr", "mcg/min", "mg/hr"),
    conversion_factors = c("mcg/kg/min" = 1, "mcg/kg/hr" = 1/60, "mg/kg/hr" = 1000/60, "mcg/min" = 1, "mg/hr" = 1000/60)
  )
)

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
                          "dopamine", "metaraminol", "dobutamine")) {
    
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
med_summaries <- as_tibble(meds_with_weights_dt) %>%
  group_by(encounter_id, med_category) %>%
  summarize(
    median_dose = median(med_dose_converted, na.rm = TRUE),
    iqr_lower   = quantile(med_dose_converted, 0.25, na.rm = TRUE),
    iqr_upper   = quantile(med_dose_converted, 0.75, na.rm = TRUE),
    .groups     = "drop"
  )

cohort_summaries <- med_summaries %>%
  group_by(med_category) %>%
  summarize(
    overall_median_dose = median(median_dose[median_dose != 0], na.rm = TRUE),
    overall_iqr_lower   = quantile(median_dose[median_dose != 0], 0.25, na.rm = TRUE),
    overall_iqr_upper   = quantile(median_dose[median_dose != 0], 0.75, na.rm = TRUE),
    n_encounters        = n_distinct(encounter_id),
    pct_encounters      = 100 * n_encounters / total_encounters,
    .groups             = "drop"
  )

write.csv(cohort_summaries, paste0( "output/table1_meds",site, '.csv'), 
          row.names = FALSE)

############################SOFA-97#############################################





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