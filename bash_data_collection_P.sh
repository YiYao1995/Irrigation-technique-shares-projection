#!/bin/bash

ml purge
ml CDO/2.0.5-gompi-2021b

# In ISIMIP3b global water sector, there are meteorological data from five Global Climate Models
gcms=("GFDL-ESM4" "IPSL-CM6A-LR" "MPI-ESM1-2-HR" "MRI-ESM2-0" "UKESM1-0-LL")
ssps=("ssp126" "ssp370" "ssp585")
years=("1996_2015" "2001_2020" "2006_2025" "2011_2030" "2016_2035" "2021_2040" "2026_2045" "2031_2050" "2036_2055" "2041_2060" "2046_2065" "2051_2070" "2056_2075" "2061_2080" "2066_2085" "2071_2090" "2076_2095" "2081_2100")

cd /data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP3b/InputData/climate/atmosphere/bias-adjusted/global/daily/historical

# Calculate the yearmean data and select the years we need (historical 1996-2014; future 2015-2100) and move the output to our directory
for gcm in ${gcms[@]}
do cd $gcm
	for file in *_pr_*1991_2000.nc
	do
		cdo yearmean $file $gcm"_1991_2000_yearmean.nc"
		cdo selyear,1996/2000 $gcm"_1991_2000_yearmean.nc" $gcm"_1996_2000_yearmean.nc"
	done
	for file in *_pr_*2001_2010.nc
	do
		cdo yearmean $file $gcm"_2001_2010_yearmean.nc"
	done
	for file in *_pr_*2011_2014.nc
	do
		cdo yearmean $file $gcm"_2011_2014_yearmean.nc"
	done
	cdo mergetime $gcm"_1996_2000_yearmean.nc" $gcm"_2001_2010_yearmean.nc" $gcm"_2011_2014_yearmean.nc" $gcm"_1996_2014_yearmean.nc"
	mv $gcm"_1996_2014_yearmean.nc" /data/brussel/vo/000/bvo00012/vsc10154/scenarios_design/PPET_data/
	cd ..
done

 
cd /data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP3b/InputData/climate/atmosphere/bias-adjusted/global/daily/ssp126
for gcm in ${gcms[@]}
do cd $gcm
	for file in *_pr_*2015_2020.nc
	do
		cdo yearmean $file $gcm"_2015_2020_yearmean.nc"
	done
	
	for file in *_pr_*2021_2030.nc
	do
		cdo yearmean $file $gcm"_2021_2030_yearmean.nc"
	done
	
	for file in *_pr_*2031_2040.nc
	do
		cdo yearmean $file $gcm"_2031_2040_yearmean.nc"
	done
	
	for file in *_pr_*2041_2050.nc
	do
		cdo yearmean $file $gcm"_2041_2050_yearmean.nc"
	done
	
	for file in *_pr_*2051_2060.nc
	do
		cdo yearmean $file $gcm"_2051_2060_yearmean.nc"
	done
	
	for file in *_pr_*2061_2070.nc
	do
		cdo yearmean $file $gcm"_2061_2070_yearmean.nc"
	done
	
	for file in *_pr_*2071_2080.nc
	do
		cdo yearmean $file $gcm"_2071_2080_yearmean.nc"
	done
	
	for file in *_pr_*2081_2090.nc
	do
		cdo yearmean $file $gcm"_2081_2090_yearmean.nc"
	done
	
	for file in *_pr_*2091_2100.nc
	do
		cdo yearmean $file $gcm"_2091_2100_yearmean.nc"
	done
	
	cdo mergetime *_yearmean.nc $gcm"_ssp126_2015-2100_yearmean.nc"
	mv $gcm"_ssp126_2015-2100_yearmean.nc" /data/brussel/vo/000/bvo00012/vsc10154/scenarios_design/PPET_data/
	cd ..
done

cd /data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP3b/InputData/climate/atmosphere/bias-adjusted/global/daily/ssp370
for gcm in ${gcms[@]}
do cd $gcm
	for file in *_pr_*2015_2020.nc
	do
		cdo yearmean $file $gcm"_2015_2020_yearmean.nc"
	done
	
	for file in *_pr_*2021_2030.nc
	do
		cdo yearmean $file $gcm"_2021_2030_yearmean.nc"
	done
	
	for file in *_pr_*2031_2040.nc
	do
		cdo yearmean $file $gcm"_2031_2040_yearmean.nc"
	done
	
	for file in *_pr_*2041_2050.nc
	do
		cdo yearmean $file $gcm"_2041_2050_yearmean.nc"
	done
	
	for file in *_pr_*2051_2060.nc
	do
		cdo yearmean $file $gcm"_2051_2060_yearmean.nc"
	done
	
	for file in *_pr_*2061_2070.nc
	do
		cdo yearmean $file $gcm"_2061_2070_yearmean.nc"
	done
	
	for file in *_pr_*2071_2080.nc
	do
		cdo yearmean $file $gcm"_2071_2080_yearmean.nc"
	done
	
	for file in *_pr_*2081_2090.nc
	do
		cdo yearmean $file $gcm"_2081_2090_yearmean.nc"
	done
	
	for file in *_pr_*2091_2100.nc
	do
		cdo yearmean $file $gcm"_2091_2100_yearmean.nc"
	done
	
	cdo mergetime *_yearmean.nc $gcm"_ssp370_2015-2100_yearmean.nc"
	mv $gcm"_ssp370_2015-2100_yearmean.nc" /data/brussel/vo/000/bvo00012/vsc10154/scenarios_design/PPET_data/
	cd ..
done

cd /data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP3b/InputData/climate/atmosphere/bias-adjusted/global/daily/ssp585
for gcm in ${gcms[@]}
do cd $gcm
	for file in *_pr_*2015_2020.nc
	do
		cdo yearmean $file $gcm"_2015_2020_yearmean.nc"
	done
	
	for file in *_pr_*2021_2030.nc
	do
		cdo yearmean $file $gcm"_2021_2030_yearmean.nc"
	done
	
	for file in *_pr_*2031_2040.nc
	do
		cdo yearmean $file $gcm"_2031_2040_yearmean.nc"
	done
	
	for file in *_pr_*2041_2050.nc
	do
		cdo yearmean $file $gcm"_2041_2050_yearmean.nc"
	done
	
	for file in *_pr_*2051_2060.nc
	do
		cdo yearmean $file $gcm"_2051_2060_yearmean.nc"
	done
	
	for file in *_pr_*2061_2070.nc
	do
		cdo yearmean $file $gcm"_2061_2070_yearmean.nc"
	done
	
	for file in *_pr_*2071_2080.nc
	do
		cdo yearmean $file $gcm"_2071_2080_yearmean.nc"
	done
	
	for file in *_pr_*2081_2090.nc
	do
		cdo yearmean $file $gcm"_2081_2090_yearmean.nc"
	done
	
	for file in *_pr_*2091_2100.nc
	do
		cdo yearmean $file $gcm"_2091_2100_yearmean.nc"
	done
	
	cdo mergetime *_yearmean.nc $gcm"_ssp585_2015-2100_yearmean.nc"
	mv $gcm"_ssp585_2015-2100_yearmean.nc" /data/brussel/vo/000/bvo00012/vsc10154/scenarios_design/PPET_data/
	cd ..
done

# cd /data/brussel/vo/000/bvo00012/vsc10154/scenarios_design/PPET_data
# Merge the data to 1996-2100
for gcm in ${gcms[@]}
do
	for ssp in ${ssps[@]}
	do
		cdo mergetime $gcm"_1996_2014_yearmean.nc" $gcm"_"$ssp"_2015-2100_yearmean.nc" $gcm"_"$ssp"_1996-2100_yearmean.nc"
	done
done

rm *_1996_2014_yearmean.nc *_2015-2100_yearmean.nc

*cp *.nc P/
# cd P/	

# Get the 20-year data and calculate the mean value
for gcm in ${gcms[@]}
do
	for ssp in ${ssps[@]}
	do
		cdo selyear,1996/2015 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_1996-2015_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_1996-2015_yearmean.nc" $gcm"_"$ssp"_1996-2015_yearmean_timmean.nc"
		
		cdo selyear,2001/2020 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2001_2020_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2001_2020_yearmean.nc" $gcm"_"$ssp"_2001_2020_yearmean_timmean.nc"
		
		cdo selyear,2006/2025 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2006_2025_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2006_2025_yearmean.nc" $gcm"_"$ssp"_2006_2025_yearmean_timmean.nc"
		
		cdo selyear,2011/2030 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2011_2030_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2011_2030_yearmean.nc" $gcm"_"$ssp"_2011_2030_yearmean_timmean.nc"
		
		cdo selyear,2016/2035 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2016_2035_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2016_2035_yearmean.nc" $gcm"_"$ssp"_2016_2035_yearmean_timmean.nc"
		
		cdo selyear,2021/2040 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2021_2040_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2021_2040_yearmean.nc" $gcm"_"$ssp"_2021_2040_yearmean_timmean.nc"
		
		cdo selyear,2026/2045 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2026_2045_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2026_2045_yearmean.nc" $gcm"_"$ssp"_2026_2045_yearmean_timmean.nc"
		
		cdo selyear,2031/2050 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2031_2050_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2031_2050_yearmean.nc" $gcm"_"$ssp"_2031_2050_yearmean_timmean.nc"
		
		cdo selyear,2036/2055 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2036_2055_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2036_2055_yearmean.nc" $gcm"_"$ssp"_2036_2055_yearmean_timmean.nc"
		
		cdo selyear,2041/2060 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2041_2060_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2041_2060_yearmean.nc" $gcm"_"$ssp"_2041_2060_yearmean_timmean.nc"
		
		cdo selyear,2046/2065 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2046_2065_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2046_2065_yearmean.nc" $gcm"_"$ssp"_2046_2065_yearmean_timmean.nc"
		
		cdo selyear,2051/2070 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2051_2070_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2051_2070_yearmean.nc" $gcm"_"$ssp"_2051_2070_yearmean_timmean.nc"
		
		cdo selyear,2056/2075 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2056_2075_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2056_2075_yearmean.nc" $gcm"_"$ssp"_2056_2075_yearmean_timmean.nc"
		
		cdo selyear,2061/2080 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2061_2080_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2061_2080_yearmean.nc" $gcm"_"$ssp"_2061_2080_yearmean_timmean.nc"
		
		cdo selyear,2066/2085 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2066_2085_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2066_2085_yearmean.nc" $gcm"_"$ssp"_2066_2085_yearmean_timmean.nc"
		
		cdo selyear,2071/2090 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2071_2090_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2071_2090_yearmean.nc" $gcm"_"$ssp"_2071_2090_yearmean_timmean.nc"
		
		cdo selyear,2076/2095 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2076_2095_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2076_2095_yearmean.nc" $gcm"_"$ssp"_2076_2095_yearmean_timmean.nc"
		
		cdo selyear,2081/2100 $gcm"_"$ssp"_1996-2100_yearmean.nc" $gcm"_"$ssp"_2081_2100_yearmean.nc"
		cdo timmean $gcm"_"$ssp"_2081_2100_yearmean.nc" $gcm"_"$ssp"_2081_2100_yearmean_timmean.nc"
	done
done

rm *yearmean.nc

# Do the remap from 0.5x0.5 to 0.9x1.25
for gcm in ${gcms[@]}
do 
	for ssp in ${ssps[@]}
	do 
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_1996-2015_yearmean_timmean.nc" $gcm"_"$ssp"_1996_2015_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2001_2020_yearmean_timmean.nc" $gcm"_"$ssp"_2001_2020_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2006_2025_yearmean_timmean.nc" $gcm"_"$ssp"_2006_2025_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2011_2030_yearmean_timmean.nc" $gcm"_"$ssp"_2011_2030_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2016_2035_yearmean_timmean.nc" $gcm"_"$ssp"_2016_2035_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2021_2040_yearmean_timmean.nc" $gcm"_"$ssp"_2021_2040_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2026_2045_yearmean_timmean.nc" $gcm"_"$ssp"_2026_2045_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2031_2050_yearmean_timmean.nc" $gcm"_"$ssp"_2031_2050_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2036_2055_yearmean_timmean.nc" $gcm"_"$ssp"_2036_2055_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2041_2060_yearmean_timmean.nc" $gcm"_"$ssp"_2041_2060_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2046_2065_yearmean_timmean.nc" $gcm"_"$ssp"_2046_2065_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2051_2070_yearmean_timmean.nc" $gcm"_"$ssp"_2051_2070_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2056_2075_yearmean_timmean.nc" $gcm"_"$ssp"_2056_2075_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2061_2080_yearmean_timmean.nc" $gcm"_"$ssp"_2061_2080_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2066_2085_yearmean_timmean.nc" $gcm"_"$ssp"_2066_2085_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2071_2090_yearmean_timmean.nc" $gcm"_"$ssp"_2071_2090_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2076_2095_yearmean_timmean.nc" $gcm"_"$ssp"_2076_2095_yearmean_timmean_remapcon.nc"
	cdo remapcon,grid_0.9x1.25_for_cdo.nc $gcm"_"$ssp"_2081_2100_yearmean_timmean.nc" $gcm"_"$ssp"_2081_2100_yearmean_timmean_remapcon.nc"
	done
done


rm *timmean.nc

cd /data/brussel/vo/000/bvo00012/vsc10154/scenarios_design/PPET_data/P

# Calculate the ensemble mean of different GCMs
for ssp in ${ssps[@]}
do 
	for year in ${years[@]}
	do 
		cdo ensmean *$ssp"_"$year* $ssp"_"$year"_remapcon_ensmean.nc"
	done
done
	
rm *remapcon.nc
