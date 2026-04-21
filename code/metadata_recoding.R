# Metadata recoding ----

## Setup ----

setwd(".")
set.seed(0509)

library(tidyverse)

## Co-morbidity group definitions ----

replacements <- list()

### Anemia ----

replacements$anemia <- c(
  "anemia" = "anemia",
  "anemia of 888 disease" = "anemia",
  "anemia of chronic disease" = "anemia",
  "iron deficiency anemia" = "anemia",
  "mild anemia" = "anemia",
  "pernicious anemia" = "anemia",
  "pernivious anemia" = "anemia"
)

### Arthritis ----

replacements$arthritis <- c(
  "arthritis" = "arthritis",
  "arthritis H knee/hang & H hip" = "arthritis",
  "polyarthritis" = "arthritis",
  "Gout" = "arthritis",
  "gout" = "arthritis",
  "gout/ hyperuricemia" = "arthritis",
  "OA" = "arthritis",
  "OA knees" = "arthritis",
  "OA+DDD lumbar spine" = "arthritis",
  "chronic pain OA" = "arthritis",
  "neuropathic pain/ OA" = "arthritis",
  "osteoarthritis" = "arthritis",
  "osteoarthritis (bilat knee replacements)" = "arthritis",
  "osteoarthritis (knees)" = "arthritis"
)

### Atherosclerosis ----

replacements$atherosclerosis <- c(
  "CAD" = "atherosclerosis",
  "CAD?MI" = "atherosclerosis",
  "coronary artery disease" = "atherosclerosis",
  "ischemic HD" = "atherosclerosis",
  "PVD" = "atherosclerosis"
)

### COPD ----

replacements$copd <- c(
  "COPD" = "copd"
)

### Dementia ----

replacements$dementia <- c(
  "Alzheimer’s" = "dementia",
  "agitation/ dementia" = "dementia",
  "dementia" = "dementia",
  "dementia- early onset/ aggressive" = "dementia",
  "mild dementia" = "dementia",
  "vascular dementia" = "dementia",
  "vascular/ Alzheimer's dementia" = "dementia"
)


### Depression ----

replacements$depression <- c(
  "Depression" = "depression",
  "chronic depression/ anxiety" = "depression",
  "depression" = "depression",
  "depression/ anxiety" = "depression"
)

### Diabetes ----

replacements$diabetes <- c(
  "1DDM" = "diabetes",
  "DM" = "diabetes",
  "DMII" = "diabetes",
  "IDDM" = "diabetes",
  "diabetes" = "diabetes",
  "diabetes II (diet only)" = "diabetes",
  "diabetes mellitus" = "diabetes",
  "diabetes millitus" = "diabetes",
  "diabetes/ diabetic retinopathy" = "diabetes",
  "diabetic neuropathy" = "diabetes",
  "diabetic peripheral neuropathy" = "diabetes",
  "type II diabetes" = "diabetes"
)

### Dyslipidemia ----

replacements$dyslipidemia <- c(
  "dislipidemia" = "dyslipidemia",
  "dyslipidemia" = "dyslipidemia",
  "hyperlipidaemia" = "dyslipidemia",
  "hyperlipidemia" = "dyslipidemia"
)

### Hypertension ----

replacements$hypertension <- c(
  "HTN" = "hypertension",
  "HTN cont'd LTIA & HTN" = "hypertension",
  "hypertension" = "hypertension"
)

### Hypothyroidism ----

replacements$hypothyroidism <- c(
  "hypothyroid" = "hypothyroidism",
  "hypothyroidism" = "hypothyroidism"
)

### Osteoporosis ----

replacements$osteoporosis <- c(
  "osteoporosis" = "osteoporosis",
  "osteoporosis #hip #ankle" = "osteoporosis"
)

### Renal disease ----

replacements$renal_disease <- c(
  "CRF" = "renal_disease",
  "Chronic renal disease" = "renal_disease",
  "chronic kidney disorder" = "renal_disease",
  "chronic renal failure" = "renal_disease",
  "renal disease" = "renal_disease",
  "horse shoe kidney" = "renal_disease"
)

## Medication group definitions ----

### ACE inhibitors ----

replacements$ace_inhibitors <- c(
  "convensyl" = "ace_inhibitors",
  "coversyl" = "ace_inhibitors",
  "lisinopril" = "ace_inhibitors",
  "mavik" = "ace_inhibitors",
  "quinapril" = "ace_inhibitors",
  "ramipril" = "ace_inhibitors",
  "amipril" = "ace_inhibitors"
)

### Acetaminophen ----

replacements$acetaminophen <- c(
  "acetaminophen" = "acetaminophen",
  "acetaminophen PRN" = "acetaminophen",
  "acetaminophen/ caffeine/codeine phosphate" =
    "acetaminophen", "acetaminophine" = "acetaminophen",
  "acetominaphen" = "acetaminophen",
  "acetaminophen" = "acetaminophen",
  "tylenol" = "acetaminophen",
  "tylenol arthritis" = "acetaminophen"
)

### Atypical antipsychotics ----

replacements$aaps <- c(
  "DMS-risperidone" = "aaps",
  "clonazepam risperidone" = "aaps",
  "quetiapine" = "aaps",
  "quetiapine PRN" = "aaps",
  "risperidone" = "aaps",
  "olanzapine" = "aaps",
  "olanzipine" = "aaps",
  "clonazepam risperidone" = "aaps"
)

### Benzodiazepines ----

replacements$benzodiazepines <- c(
  "Lorazepam" = "benzodiazepines",
  "baclofen" = "benzodiazepines",
  "chlordiazepoxide" = "benzodiazepines",
  "clonazepam" = "benzodiazepines",
  "clonazepam PRN" = "benzodiazepines",
  "clonazepam risperidone" = "benzodiazepines",
  "clonazepam" = "benzodiazepines",
  "clonazepam" = "benzodiazepines",
  "lorazepam" = "benzodiazepines",
  "lorazepam PRN" = "benzodiazepines",
  "novo-clobazam" = "benzodiazepines",
  "temazepam" = "benzodiazepines"
)

### Beta blockers ----

replacements$beta_blockers <- c(
  "Carvedilol" = "beta_blockers",
  "acebulotol" = "beta_blockers",
  "atenolol" = "beta_blockers",
  "medoprodol" = "beta_blockers",
  "metoprolol" = "beta_blockers",
  "metrapolol" = "beta_blockers",
  "metropdol" = "beta_blockers",
  "nadolol" = "beta_blockers",
  "sandoz-bispropolol" = "beta_blockers"
)

### Blood thinners/antiplatalets ----

replacements$antiplatelets <- c(
  "eliquis" = "antiplatelets",
  "elmiron" = "antiplatelets",
  "warfarin" = "antiplatelets",
  "clopidogrel" = "antiplatelets",
  "plavix" = "antiplatelets",
  "ASA" = "antiplatelets",
  "Asa" = "antiplatelets",
  "EC ASA" = "antiplatelets",
  "ECASA" = "antiplatelets",
  "ELAASA" = "antiplatelets",
  "cikcgucubem ASA" = "antiplatelets",
  "adalat xl" = "antiplatelets"
)

### Calcium channel blockers ----

replacements$ccbs <- c(
  "amlodipine" = "ccbs",
  "amlopid" = "ccbs",
  "adalat xl" = "ccbs",
  "mfedipine" = "ccbs",
  "norvasc" = "ccbs"
)

### Diuretics ----

replacements$diuretics <- c(
  "furosemide" = "diuretics",
  "furosimide" = "diuretics",
  "hydrochlorothiazide" = "diuretics"
)

### Insulin ----

replacements$insulin <- c(
  "humilin N" = "insulin",
  "humilin R" = "insulin",
  "humulin N" = "insulin",
  "humulin R" = "insulin",
  "lantus" = "insulin",
  "nolovlin toronto" = "insulin",
  "novolin" = "insulin",
  "novolin NPH" = "insulin",
  "novolin toronto" = "insulin",
  "novorapid" = "insulin"
)

### Iron supplements ----

replacements$iron_supplements <- c(
  "FeSO4" = "iron_supplements",
  "feramax" = "iron_supplements",
  "ferramax" = "iron_supplements",
  "ferrous gluconate" = "iron_supplements",
  "ferrous solufate" = "iron_supplements",
  "ferrous sulfate" = "iron_supplements",
  "palafer iron" = "iron_supplements"
)

### Laxatives ----

replacements$laxatives <- c(
  "lactulose" = "laxatives",
  "lactulose syrup" = "laxatives",
  "lax a day" = "laxatives",
  "lax-a-day" = "laxatives",
  "bisacodyl" = "laxatives",
  "senekot" = "laxatives",
  "senekot-s PRN" = "laxatives",
  "sennokot" = "laxatives",
  "senokot" = "laxatives"
)

### Levothyroxine ----

replacements$levothyroxine <- c(
  "Levethyroxine" = "levothyroxine",
  "synthroid" = "levothyroxine"
)

### Metformin ----

replacements$metformin <- c(
  "glumetza" = "metformin",
  "metformin" = "metformin",
  "metparmin" = "metformin"
)

### Nitrates ----

replacements$nitrates <- c(
  "isosorbide" = "nitrates",
  "isosorbide-s-mononitrate" = "nitrates",
  "nitro patch" = "nitrates",
  "mylan nitro" = "nitrates",
  "nitro pump spray" = "nitrates",
  "nitroglycerin" = "nitrates",
  "nitroglycerin spray PRN" = "nitrates",
  "nitroglycerine" = "nitrates",
  "nitroglycerine PRN" = "nitrates",
  "nitropatch" = "nitrates",
  "nitrospray" = "nitrates"
)

### Opioids ----

replacements$opioids <- c(
  "hydromorph cortin" = "opioids",
  "hydromorphone" = "opioids",
  "hydromorphone PRN" = "opioids",
  "hydromorphone contin" = "opioids",
  "loperamide" = "opioids",
  "loperamide PRN" = "opioids",
  "loperamide" = "opioids",
  "methadone" = "opioids",
  "Ma ESLON" = "opioids",
  "morphine ER" = "opioids",
  "morphine sulfate" = "opioids",
  "oxycodone" = "opioids",
  "statex" = "opioids",
  "statex PRN" = "opioids",
  "codeine phosphate PRN" = "opioids",
  "acetaminophen/ caffeine/codeine phosphate" = "opioids",
  "codeine" = "opioids",
  "codeine phosphate" = "opioids"
)

### Proton pump inhibitors ----

replacements$ppis <- c(
  "dantoprazole" = "ppis",
  "omeprazole" = "ppis",
  "omprazole" = "ppis",
  "pantoprazole" = "ppis",
  "pantoprozole" = "ppis",
  "rabeprazole" = "ppis",
  "rabeprozole" = "ppis",
  "tecta" = "ppis"
)

### SSRIs ----

replacements$ssris <- c(
  "citalopram" = "ssris",
  "pms-citalopram" = "ssris",
  "mylon-sentaline" = "ssris",
  "paroxetine" = "ssris",
  "sertraline" = "ssris"
)

### Statins ----

replacements$statins <- c(
  "atorvastatin" = "statins",
  "co-atorvastatin" = "statins",
  "rosuvastatin" = "statins",
  "simvastatin" = "statins",
  "sinvastatin" = "statins"
)

### Tricyclic antidepressants ----

replacements$tcas <- c(
  "desipramine" = "tcas",
  "amitriptyline" = "tcas",
  "aventyl" = "tcas"
)

### Vitamin B12 ----

replacements$vitb12 <- c(
  "high potency vit B12" = "vitb12",
  "vit B 12" = "vitb12",
  "vit B12" = "vitb12",
  "vit b12" = "vitb12"
)

### Vitamin D ----

replacements$vitd <- c(
  "vid d" = "vitd",
  "vit D" = "vitd",
  "vit d" = "vitd",
  "vitamin D" = "vitd"
)

## Generate recoded metadata file ----

# Define which groups come from diseases vs medications
disease_groups <- c(
  "anemia", "arthritis", "atherosclerosis", "copd",
  "dementia", "depression", "diabetes", "dyslipidemia",
  "hypertension", "hypothyroidism", "osteoporosis",
  "renal_disease"
)

medication_groups <- c(
  "ace_inhibitors", "acetaminophen", "aaps",
  "benzodiazepines", "beta_blockers", "antiplatelets",
  "ccbs", "diuretics", "insulin", "iron_supplements",
  "laxatives", "levothyroxine", "metformin", "nitrates",
  "opioids", "ppis", "ssris", "statins", "tcas",
  "vitb12", "vitd"
)

# Read raw metadata
raw_metadata <- read_tsv("data/metadata/patient_data.tsv")

# Function to check if any term in a replacement vector matches a subject's
# comma-separated disease or medication string
has_match <- function(raw_string, replacement_vec) {
  if (is.na(raw_string)) {
    return(0L)
  }
  # Split on comma, trim whitespace
  items <- trimws(unlist(strsplit(raw_string, ",")))
  # Check if any item matches any key in the replacement vector
  as.integer(any(items %in% names(replacement_vec)))
}

# Build the recoded dataframe
recoded <- raw_metadata

# Add disease group columns
for (grp in disease_groups) {
  col_name <- grp
  recoded[[col_name]] <- sapply(
    raw_metadata$ProbName,
    has_match,
    replacement_vec = replacements[[grp]]
  )
}

# Add medication group columns
for (grp in medication_groups) {
  col_name <- grp
  recoded[[col_name]] <- sapply(
    raw_metadata$MedName,
    has_match,
    replacement_vec = replacements[[grp]]
  )
}

recoded_metadata <- recoded

# Write output
write_tsv(recoded_metadata, "data/metadata/patient_data.tsv")
