pct_irr = zeros(720, 360, 86);
pct_irr_drip = zeros(720, 360, 86);
pct_irr_spri = zeros(720, 360, 86);
pct_irr_floo = zeros(720, 360, 86);

file_surf = 'F:\landuse-15crops_magpie_ssp1_annual_2015_2100.nc';

temperate_cereals_irrigated = ncread(file_surf, 'temperate_cereals_irrigated'); %sprinkler
rice_irrigated = ncread(file_surf, 'rice_irrigated'); %flood
maize_irrigated = ncread(file_surf, 'maize_irrigated'); %sprinkler
tropical_cereals_irrigated = ncread(file_surf, 'tropical_cereals_irrigated'); %sprinkler
pulses_irrigated = ncread(file_surf, 'pulses_irrigated'); %drip
temperate_roots_irrigated = ncread(file_surf, 'temperate_roots_irrigated');%sprinkler
tropical_roots_irrigated = ncread(file_surf, 'tropical_roots_irrigated');%sprinkler
oil_crops_sunflower_irrigated = ncread(file_surf, 'oil_crops_sunflower_irrigated');%drip
oil_crops_soybean_irrigated = ncread(file_surf, 'oil_crops_soybean_irrigated');%drip
oil_crops_groundnut_irrigated = ncread(file_surf, 'oil_crops_groundnut_irrigated');%sprinkler
oil_crops_rapeseed_irrigated = ncread(file_surf, 'oil_crops_rapeseed_irrigated');%sprinkler
c4per_irrigated_food = ncread(file_surf, 'c4per_irrigated_food');%sprinkler
others_c3ann_irrigated = ncread(file_surf, 'others_c3ann_irrigated');%drip
others_c3nfx_irrigated = ncread(file_surf, 'others_c3nfx_irrigated');%drip
c3per_irrigated_food = ncread(file_surf, 'c3per_irrigated_food');%sprinkler
c3per_irrigated_bf = ncread(file_surf, 'c3per_irrigated_bf');%sprinkler
c4per_irrigated_bf = ncread(file_surf, 'c4per_irrigated_bf');%sprinkler

pct_irr = pct_irr + temperate_cereals_irrigated + rice_irrigated + ...,
    maize_irrigated + tropical_cereals_irrigated + pulses_irrigated + ...,
    temperate_roots_irrigated + tropical_roots_irrigated + oil_crops_sunflower_irrigated + oil_crops_soybean_irrigated + ...,
    oil_crops_groundnut_irrigated + oil_crops_rapeseed_irrigated + c4per_irrigated_food + ...,
    others_c3ann_irrigated + others_c3nfx_irrigated + c3per_irrigated_food + ...,
    c3per_irrigated_bf + c4per_irrigated_bf;

pct_irr_drip = pulses_irrigated + ...,
    oil_crops_sunflower_irrigated + oil_crops_soybean_irrigated + ...,
    others_c3ann_irrigated + others_c3nfx_irrigated;

pct_irr_floo = pct_irr_floo + rice_irrigated;

pct_irr_spri = pct_irr - pct_irr_drip - pct_irr_floo;

save('C:\Users\yiyao.F-IR-HYDRPC116\OneDrive - Vrije Universiteit Brussel\Data_for_ISIMIP\ssp1_magpie_noadapt\pct_irr','pct_irr')
save('C:\Users\yiyao.F-IR-HYDRPC116\OneDrive - Vrije Universiteit Brussel\Data_for_ISIMIP\ssp1_magpie_noadapt\pct_irr_spri','pct_irr_spri')
save('C:\Users\yiyao.F-IR-HYDRPC116\OneDrive - Vrije Universiteit Brussel\Data_for_ISIMIP\ssp1_magpie_noadapt\pct_irr_drip','pct_irr_drip')
save('C:\Users\yiyao.F-IR-HYDRPC116\OneDrive - Vrije Universiteit Brussel\Data_for_ISIMIP\ssp1_magpie_noadapt\pct_irr_floo','pct_irr_floo')

opt_drip = pct_irr_drip ./ pct_irr;
opt_spri = pct_irr_spri ./ pct_irr;
opt_floo = pct_irr_floo ./ pct_irr;

save('C:\Users\yiyao.F-IR-HYDRPC116\OneDrive - Vrije Universiteit Brussel\Data_for_ISIMIP\ssp1_magpie_noadapt\opt_spri','opt_spri')
save('C:\Users\yiyao.F-IR-HYDRPC116\OneDrive - Vrije Universiteit Brussel\Data_for_ISIMIP\ssp1_magpie_noadapt\opt_drip','opt_drip')
save('C:\Users\yiyao.F-IR-HYDRPC116\OneDrive - Vrije Universiteit Brussel\Data_for_ISIMIP\ssp1_magpie_noadapt\opt_floo','opt_floo')



file_surf = 'C:\Research2\ISIMIP\cell_share_irrigation_systems.nc';
act_irr_fra = ncread(file_surf,'Cell share irrigation systems');
act_irr_fra_floo = act_irr_fra(:,:,1);
act_irr_fra_spri = act_irr_fra(:,:,2);
act_irr_fra_drip = act_irr_fra(:,:,3);
act_irr_fra_alle = act_irr_fra_floo + act_irr_fra_spri + act_irr_fra_drip;
act_irr_floo = act_irr_fra_floo ./ act_irr_fra_alle;
act_irr_spri = act_irr_fra_spri ./ act_irr_fra_alle;
act_irr_drip = act_irr_fra_drip ./ act_irr_fra_alle;

save('C:\Users\yiyao.F-IR-HYDRPC116\OneDrive - Vrije Universiteit Brussel\Data_for_ISIMIP\ssp1_magpie_noadapt\act_irr_spri','act_irr_spri')
save('C:\Users\yiyao.F-IR-HYDRPC116\OneDrive - Vrije Universiteit Brussel\Data_for_ISIMIP\ssp1_magpie_noadapt\act_irr_drip','act_irr_drip')
save('C:\Users\yiyao.F-IR-HYDRPC116\OneDrive - Vrije Universiteit Brussel\Data_for_ISIMIP\ssp1_magpie_noadapt\act_irr_floo','act_irr_floo')
