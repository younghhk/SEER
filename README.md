
> 🧬 **For additional cancer research software and tools**, visit  
> [Cancer Research Software Repository](https://github.com/younghhk/NCI)


## IBM-Calc: A Computational Tool for Estimating Age-Standardized Incidence-Based Mortality and Rate Ratio Confidence Intervals


This repository provides R functions to estimate **Age-Adjustd Incidence-Based Mortality (IBM) rates** and **Rate Ratios (RRs)** using tabulated data formatted like those from the **Surveillance, Epidemiology, and End Results (SEER) Program** of the U.S. National Cancer Institute.

🔗 **Learn more about SEER:** [https://seer.cancer.gov/](https://seer.cancer.gov/)





##  Features

- Compute **Incidence-Based Mortality (IBM) rates**, which are calculated by **identifying deaths among patients previously diagnosed with cancer** in population-based registries (e.g., [SEER](https://seer.cancer.gov/)).  
 

- Estimate **Rate Ratios (RRs)** to compare IBM rates between groups (e.g., **Non-Hispanic White vs. Non-Hispanic Black**) with appropriate **variance and confidence interval** calculations based on the Delta method.

- Fully compatible with **SEER-style tabulated data** (age-adjusted or age-specific rates) and designed for reproducibility in **R**.



## Files
- `IBM.R` — main functions for IBM rate and rate ratio estimation  
- `age_adjusted_data_grace.xlsx` — example dataset layout (age-stratified SEER-like data)



## Method Overview

### Age-Adjusted Incidence-Based Mortality (IBM) Rate

The age-adjusted IBM rate is computed as:

    R'_IBM = Σ_i [ (D_i / P_i) * w_i ] * 100,000

where

    w_i = stdpop_i / Σ_j stdpop_j

**Definitions:**
- `D_i`: deaths among incident cancer cases in age group *i*  
- `P_i`: person-years among those incident cases in age group *i*  
- `w_i`: normalized standard-population weight for age group *i*



### Confidence Intervals for Age-Adjusted Rates

Confidence intervals (CIs) are calculated using SEER’s **Fay–Feuer** or **Tiwari** methods.  
Both assume Poisson-distributed counts and use the directly standardized rate (DSR) variance:

    v = Σ_i (w_i² × D_i)



#### Fay–Feuer Method (recommended for very small counts)

Identifies the stratum with the largest weight:

    w_m = max(w_i)
    z = w_m²

Lower and upper confidence limits (using chi-square quantiles):

    Lower = (v / (2 * R')) * χ²(α/2, df = 2 * R'² / v)
    Upper = ((v + z) / (2 * (R' + w_m))) * χ²(1 - α/2, df = 2 * (R' + w_m)² / (v + z))



#### Tiwari et al. Modification (for moderate counts)

The Tiwari method replaces the single maximum weight with averages across
non-zero strata, providing smoother, less conservative confidence intervals:

    w_m = avg(w_i)   # restricted to strata with P_i > 0 or D_i > 0
    z   = avg(w_i²)

Then apply the same formulas for the lower and upper bounds:

    Lower = (v / (2 * R')) * χ²(α/2, df = 2 * R'² / v)
    Upper = ((v + z) / (2 * (R' + w_m))) * χ²(1 - α/2, df = 2 * (R' + w_m)² / (v + z))



### Notes

- Implementation follows SEER’s official documentation:  
  [SEER*Stat Rate Algorithms](https://seer.cancer.gov/help/seerstat/equations-and-algorithms/rate-algorithms)



### Rate Ratios (RR)
The **Rate Ratio (RR)** compares the IBM rates between two groups — for example,  
**Group 1:** Non-Hispanic White,  
**Group 2:** Non-Hispanic Black.


`RR = (IBM Rate of Group 2) / (IBM Rate of Group 1)`

- The variance of `log(RR)` is approximated using the **Delta method**.  
- 95% CIs are computed on the log scale and exponentiated back.

**Interpretation:**
- RR = 1 → no difference in IBM rates between groups  
- RR > 1 → Group 2 has a higher IBM rate  
- RR < 1 → Group 2 has a lower IBM rate  



## 🔗 References
- **[Age-adjusted rate confidence intervals (SEER Documentation)](https://seer.cancer.gov/help/seerstat/equations-and-algorithms/rate-algorithms)**  

---

##  Example Usage

```r
# 1. Load the function file
source("IBM.R")

# 2. Install dependencies
#install.packages("readxl")
library(readxl)

# 3. Read age-stratified data
df <- read_excel("age_adjusted_data_grace.xlsx")

# 4. Define two comparison groups
# Example: ER-negative, Non-Hispanic White vs Non-Hispanic Black, age 30–54
idx1 <- which(df$ER == "Negative" &
              df$Race == "Non-Hispanic White" &
              df$age_group_strata == "30 - 54")

idx2 <- which(df$ER == "Negative" &
              df$Race == "Non-Hispanic Black" &
              df$age_group_strata == "30 - 54")

# 5. Compute IBM rates and RR 

out <- compute_ibm_rates_and_ratio(
  df,
  idx1 = idx1,
  idx2 = idx2,
  label_group1 = "ER- & NHW & 30-54",
  label_group2 = "ER- & NHB & 30-54",
  ci_method = "fayfeuer"  # or "tiwari"
)

# IBM rates table (per group)
out$rate

# Rate ratio table (Group2 / Group1)
out$rate_ratio
```



---

## 📜 Citation
If you use this code,  please cite the GitHub repository as follows:

Vo, J., & Hong, G. (2025). Age-Adjusted Incidence-Based Mortality and Rate Ratio for Redistributed SEER Data (Version 1.0.0)]. GitHub

---

## 🔒 Access & Contact

SEER data are restricted to authorized collaborators.
To request access or contribute to related tools, contact:

**Grace Hong**  
📧 [grace.hong@nih.gov](mailto:grace.hong@nih.gov)



