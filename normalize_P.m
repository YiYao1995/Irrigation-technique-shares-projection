str_begin = 'C:\Users\yiyao.F-IR-HYDRPC116\OneDrive - Vrije Universiteit Brussel\Data_for_ISIMIP\Precipitation\GFDL-ESM4_ssp585_'
str_year_list = ["1996_2015", "2001_2020", "2006_2025", "2011_2030", "2016_2035", "2021_2040", "2026_2045", "2031_2050", "2036_2055", "2041_2060", "2046_2065", "2051_2070", "2056_2075", "2061_2080", "2066_2085", "2071_2090", "2076_2095", "2081_2100"]
str_end = '_yearmean_timmean.nc'

str_output_end = '_after_normalization.csv'


for year = 1 : 18
    str_file = strcat(str_begin, str_year_list(year), str_end);
    pr = ncread(str_file, 'pr');
    pr_new = (pr-5.59E-07)/(0.000113881507786572-5.59E-07)
    str_output = strcat(str_begin, str_year_list(year), str_output_end);
    csvwrite(str_output, pr_new)
end