
> 🧬 **For additional cancer research software and tools**, visit  
> [NCI Cancer Research Software Repository](https://github.com/younghhk/NCI)


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

$$
R'_{\text{IBM}}=\sum_{i=x}^{y}\left(\frac{D_i}{P_i}\right) w_i \times 100{,}000
$$

with

$$
w_i=\frac{\mathrm{stdpop}_i}{\sum_{j=x}^{y}\mathrm{stdpop}_j}
$$

**Where:**
- \(D_i\): deaths among incident cancer cases in age group \(i\)  
- \(P_i\): person-years among those incident cases in age group \(i\)  
- \(w_i\): normalized standard-population weight for age group \(i\)


| Concept | SEER Notation | IBM Notation | Description |
|----------|----------------|--------------|--------------|
| Event count (cases or deaths) | `count_i` | `D_i` | Number of **deaths among incident cases** in age group *i* |
| Population at risk | `pop_i` | `P_i` | **Person-years** among incident cancer cases (denominator for the rate) |
| Standard population weight | `stdpop_i` | `w_i` (after normalization) | Proportion of the standard population represented by age group *i* |
| Crude rate | `cruderate` | `r_i = D_i / P_i` | Age-specific IBM rate |
| Age-adjusted rate | `aarate_x–y` | `R' = Σ w_i r_i` | Weighted average of age-specific IBM rates |


Adjustments for small sample sizes can be made using:
- **Fay–Feuer method** – recommended when **event counts are very low or zero in some strata** (common for rare cancers)  
- **Tiwari’s modification** – performs better when **event counts are small but non-zero**, offering improved confidence interval coverage


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
install.packages("readxl")
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

# 5. Compute IBM rates and RR using Fay–Feuer method
res_ff <- compute_dsr_and_rr_for_subset(
  df, idx1, idx2,
  label1 = "ER- & NHW & 30–54",
  label2 = "ER- & NHB & 30–54",
  set_name = "Subset",
  ci_method = "fayfeuer"
)

# 6. Optionally, compute using Tiwari CI
res_ti <- compute_dsr_and_rr_for_subset(
  df, idx1, idx2,
  label1 = "ER- & NHW & 30–54",
  label2 = "ER- & NHB & 30–54",
  set_name = "Subset",
  ci_method = "tiwari"
)

# View results
res_ff$rates       # IBM rates (per 100,000) with CI for each group
res_ff$rr_of_dsrs  # RR and CI comparing the two groups
```

## 📊 Output Interpretation

### `$rates`

| **Field** | **Description** |
|------------|----------------|
| `Rate_per1e5` | Age-adjusted Incidence-Based Mortality (IBM) rate per 100,000 persons |
| `CI_low`, `CI_high` | 95% confidence interval limits |
| Interpretation | Higher values indicate a greater mortality burden among incident cancer cases |



### `$rr_of_dsrs`

| **Field** | **Description** |
|------------|----------------|
| `Estimate` | Rate Ratio (RR) = (IBM rate of Group 2) / (IBM rate of Group 1) |
| `CI_low`, `CI_high` | 95% confidence interval for the RR |

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



