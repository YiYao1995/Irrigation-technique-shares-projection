# Irrigation-technique-shares-projection


This repository is for the generation of future irrigation techniques share data (for ISIMIP specifically)

**Step 1: Getting the datasets we need:**

    (i) land use time series data
    (ii) precipitation data 


**Step 2: Getting the basic information of the irrigated crops, including the current and the optimal irrigation techniques share**

    using calc_pct_irr_tech_start.m 

We separate different CFTs into different groups: drip CFTs, sprinkler CFTs, and flood CFTs based on Jaegermeyr et al. (2015)

**Step 3: Precipitation data post processing**

    using bash_data_collection_P.sh

