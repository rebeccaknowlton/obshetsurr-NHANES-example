# Get obshetsurr functions, load libraries
setwd("C:/Users/rkkno/Documents/University of Texas at Austin/Observational/Simulation files for cluster")
source("obs.het.surr.R") # main function
source("estimate.PTE.R") # point estimates for delta, delta.s, and R.s
source("boot.var.R") # bootstrap variance
library(mgcv) # for gam
library(grf) # for trees
library(matrixStats)
library(Rsurrogate)
library(tidyverse)
library(quantreg)
library(haven)  # For reading .XPT files
library(patchwork)  

# Read in data
setwd("C:/Users/rkkno/Documents/University of Texas at Austin/Observational/NHANES Data")

# G = obesity (defined by BMI >= 30)
# S = HbA1c (non-fasting)
# Y = plasma fasting glucose
bmi <- read_xpt("BMX_L.XPT")  # Body measures file, contains BMI = "BMXBMI"
hbA1c <- read_xpt("GHB_L.XPT")  # Glycohemoglobin (%) file, contains HbA1c (non-fasting) = "LBXGH" 
glucose <- read_xpt("GLU_L.XPT")  # Plasma fasting glucose file, contains fasting glucose = "LBXGLU"

# Possible covariates of interest
demographics <- read_xpt("DEMO_L.XPT")  # Demographic variables file, contains age = "RIDAGEYR", sex = "RIAGENDR", race = "RIDRETH1"
cholesterol <- read_xpt("TCHOL_L.XPT")  # Cholesterol - total file, contains total cholesterol = "LBXTC"
smoking <- read_xpt("SMQ_L.XPT")  # Smoking - cigarette use file, contains smoked at least 100 cigarettes in life = "SMQ020", currently smoke = "SMQ040"
blood_pressure <- read_xpt("BPXO_L.XPT")  # Blood pressure - oscillometric measurements file, contains Systolic BP (1st oscillometric reading) = "BPXOSY1"
physical_activity <- read_xpt("PAQ_L.XPT") # Physical activity file, contains frequency of moderate LTPA = "PAD790Q" (filter out values in the thousands, those are codes for refused/don't know)
alcohol <- read_xpt("ALQ_L.XPT") # Alcohol use file, contains how often drink in past 12 months = "ALQ121", average # drinks/day in past 12 months = "ALQ130" (filter out values in hundreds)

# Merge datasets using SEQN (respondent sequence number)
nhanes_data <- bmi %>%
  left_join(hbA1c, by = "SEQN") %>%
  left_join(glucose, by = "SEQN") %>%
  left_join(demographics, by = "SEQN") %>%
  left_join(cholesterol, by = "SEQN") %>%
  left_join(smoking, by = "SEQN") %>%
  left_join(blood_pressure, by = "SEQN") %>%
  left_join(physical_activity, by = "SEQN") %>%
  left_join(alcohol, by = "SEQN")
# Create categorical variable for obese BMI (our exposure of interest)
nhanes_data <- nhanes_data %>%
  mutate(obese = ifelse(BMXBMI >= 30, 1, 0))
# factor other categorical confounders
nhanes_data$RIAGENDR <- factor(nhanes_data$RIAGENDR)
nhanes_data$RIDRETH1 <- factor(nhanes_data$RIDRETH1)
nhanes_data$SMQ020 <- factor(nhanes_data$SMQ020)
nhanes_data$SMQ040 <- factor(nhanes_data$SMQ040)

# For this analysis, we chose to use the covariates of age, sex, and total cholesterol
# NOTE: some of the variables have codes for nonresponses that must be filtered out, if you decide to include them
#   - if you decide to use SMQ020, filter to only include "1" and "2"
#   - PAD790Q should be filtered for < 999
#   - ALQ130 should be filtered for < 99
# Also: there was lots of missingness for smoking/alcohol variables, so be careful about including those

# Select and rename the desired variables for final dataset
final_data <- as.data.frame(nhanes_data %>%
  select(
    G = obese, # obese measured as BMI >= 30
    S = LBXGH,   # Non-fasting HbA1c
    Y = LBXGLU,  # Fasting blood plasma glucose
    age = RIDAGEYR,
    sex = RIAGENDR,
    total_cholesterol = LBXTC,
  ) %>%
  # Drop rows with missing data
  drop_na()
)

nrow(final_data)
# n = 3476 
table(final_data$G)
# n0 = 2158, n1 = 1318

# simple naive treatment effect
t.test(final_data$Y[final_data$G==1], final_data$Y[final_data$G==0]) 
# delta = 10.5726
# p < 2.2e-16 (highly significant)

# naive PTE
R.s.estimate(sone = final_data$S[final_data$G==1], szero = final_data$S[final_data$G==0], yone = final_data$Y[final_data$G==1], yzero = final_data$Y[final_data$G==0], type = "robust", extrapolate = T)
# R = 0.81

# apply proposed methods
n = nrow(final_data)
n.test = 350 # holdout ~ 10% of data
set.seed(1)
test.idx <- sample(1:n, size = n.test, replace = FALSE)
data.test <- final_data[test.idx,]
data.train <- final_data[-test.idx,]

obs.est = obs.het.surr(df.train = data.train, df.test = data.test, type = "linear", var.want = TRUE)


# Distribution of estimated PTE overall
p1 <- ggplot(obs.est, aes(x = R.s)) + 
  geom_histogram(fill = "steelblue", color = "black", bins = 30) + 
  theme_minimal() +
  labs(title = "Distribution of Estimated PTE", x = "PTE (R.s)", y = "Count")

# Plot PTE vs total cholesterol with confidence intervals
p2 <- ggplot(obs.est, aes(x = total_cholesterol, y = R.s)) + 
  geom_point(color = "steelblue", alpha = 0.7) + 
  theme_minimal() +
  labs(title = "PTE vs. Total Cholesterol", x = "Total Cholesterol (mg/dL)", y = "PTE (R.s)")

# Plot PTE vs age with confidence intervals
p3 <- ggplot(obs.est, aes(x = age, y = R.s)) + 
  geom_point(color = "steelblue", alpha = 0.7) + 
  theme_minimal() +
  labs(title = "PTE vs. Age", x = "Age", y = "PTE (R.s)")

# Plot PTE vs sex
p4 <- ggplot(obs.est, aes(x = as.factor(sex), y = R.s)) + 
  geom_boxplot(fill = "steelblue", alpha = 0.7) + 
  theme_minimal() +
  labs(title = "PTE vs. Sex", x = "Sex", y = "PTE (R.s)")

# Combine plots in a 2x2 layout
(p1 | p2) / (p3 | p4)


# new patients
new_data <- data.frame(age = c(65, 45, 35, 50, 70, 30),
                       sex = factor(c(1, 2, 2, 1, 1, 2)),
                       total_cholesterol = c(160, 220, 250, 180, 140, 210))
# patient 1: 65 yr old male with cholesterol 160
# patient 2: 45 yr old female with cholesterol 220
# patient 3: 35 yr old female with cholesterol 250
# patient 2: 50 yr old male with cholesterol 180
# patient 2: 70 yr old male with cholesterol 140
# patient 2: 30 yr old female with cholesterol 210

results <- obs.het.surr(df.train = data.train, df.test = new_data, type = "linear", var.want = TRUE)
# Create formatted table
formatted_table <- data.frame(
  Patient_ID = 1:nrow(results),
  Age = results$age,
  Sex = ifelse(results$sex == 1, "Male", "Female"),
  Cholesterol = results$total_cholesterol,
  CI = paste0("(", sprintf("%.2f", results$R.s.lower), ", ", sprintf("%.2f", results$R.s.upper), ")")
)
colnames(formatted_table) <- c("Patient ID", "Age", "Sex", "Total Cholesterol (mg/dL)", "Estimated 95% Confidence Interval for PTE")

# Save formatted table
print(formatted_table)
latex.table(as.matrix(formatted_table), "NHANES_table", caption = "", rownames = F)

