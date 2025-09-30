# SOFTWARE for SEER Cancer Data Analysis

## Incidence-Based Mortality (IBM) Rate and Rate Ratio

This repository provides R functions to compute **Incidence-Based Mortality (IBM) rates** and **rate ratios**.  
The methods include adjustments for small-count bias and variance estimation, making them suitable for rare cancers and small populations.

---

##  IBM Rate
The **IBM rate** links deaths to incident cancer cases in the registry.  
Adjustments for small counts are available:
- **Fay–Feuer method** (recommended for rare events)  
- **Tiwari’s modification** (alternative adjustment)  

Rates can be age-adjusted if needed.

---

##  Rate Ratios
Rate ratios compare IBM rates across groups (e.g., sex, race, calendar period).  
- Variances are estimated with the **Delta method**, which approximates the standard error of the log rate ratio.  
- Confidence intervals are computed on the log scale and then exponentiated.

---

##  Resources

- [Age-adjusted Rate Confidence Intervals (SEER Documentation)](https://seer.cancer.gov/help/seerstat/equations-and-algorithms/rate-algorithms)

- **Rate Ratios**  
  - Confidence interval formula:  
    Fay MP. *Approximate confidence intervals for rate ratios from directly standardized rates with sparse data.*  
    Communications in Statistics: Theory and Methods. 1999; 28(9):2141–2160.  
  - P-value formula:  
    Fay MP, Tiwari RC, Feuer EJ, Zou Z. *Estimating average annual percent change for disease rates without assuming constant change.*  
    Biometrics. 2006; 62(3):847–854.

- [R code](IBM.R)


---

##  Example Usage

```r
# Load functions
source("IBM.R")

# Install package if needed
install.packages("readxl")

# Load the library
library(readxl)

# Read the Excel file
df <- read_excel("age_adjusted_data_grace.xlsx")

# Example: Compute IBM rate with Fay–Feuer method
ibm_rate <- compute_dsr_and_rr_for_subset(
  df, idx1, idx2,
  "ER- & NHW & 30-54",
  "ER- & NHB & 30-54",
  "Subset",
  ci_method = "fayfeuer"
)

####################################################
# Output interpretation
####################################################

# The function returns a list with two elements:
#
# 1) $rates : Group-specific age-adjusted incidence-based mortality (IBM) rates
#    - DSR (Directly Standardized Rate) is computed per unit population,
#      then scaled to "Rate_per1e5" (per 100,000 people).
#    - CI_low and CI_high give the 95% confidence interval for the rate.
#
# Interpretation of Rates:
#   - The "Rate_per1e5" column represents the IBM rate per 100,000 people,
#     adjusted for age using standard weights.
#   - Higher rates indicate a greater burden of mortality linked to cancer incidence.
#   - The confidence interval (CI_low, CI_high) reflects statistical uncertainty.
#     Narrow intervals mean more precision; wide intervals occur with smaller samples.
#
# 2) $rr_of_dsrs : Rate ratio results (comparison of two groups)
#    - Estimate is the ratio of Group2's IBM rate to Group1's IBM rate.
#    - CI_low and CI_high give the 95% confidence interval for the RR.
#
# Interpretation of Rate Ratios (RR):
#   - RR = 1.0 → The two groups have the same IBM rate.
#   - RR > 1.0 → Group2 has a higher IBM rate than Group1.
#   - RR < 1.0 → Group2 has a lower IBM rate than Group1.
#
# Example:
#   Suppose Group1 has Rate_per1e5 = 3.0 and Group2 has Rate_per1e5 = 6.6.
#   Then RR ≈ 2.2, meaning Group2’s IBM rate is about 2.2 times higher than Group1’s.
#   If the CI does not include 1.0, the difference is statistically significant.

```
## ℹ️ Project Info

- **Last Updated:** September 29, 2025  
- **Contact:** [grace.hong@nih.gov](mailto:grace.hong@nih.gov) (please reach out with questions or issues)  
- **Main Repository:** [younghhk/NCI](https://github.com/younghhk/NCI)  
