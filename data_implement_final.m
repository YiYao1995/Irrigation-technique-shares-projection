str_data = "D:\ISIMIP\For_generation\landuse-15crops_image_gfdl-esm4_ssp126_annual_2015_2100"

load(append(str_data, "\grid_frac_ssp126.mat"))
basic_flood = grid_frac(:,:,1,1);
basic_sprinkler = grid_frac(:,:,1,2);
basic_drip = grid_frac(:,:,1,3);
basic_flood = 1.0000 - basic_sprinkler - basic_drip;

load(str_data + '\pct_irr.mat')
load(str_data + '\pct_irr_spri.mat')
load(str_data + '\pct_irr_drip.mat')
load(str_data + '\pct_irr_floo.mat')

load(str_data + '\opt_spri.mat')
opt_spri(find(isnan(opt_spri)==1)) = 0;
load(str_data + '\opt_drip.mat')
opt_drip(find(isnan(opt_drip)==1)) = 0;
load(str_data + '\opt_floo.mat')
opt_floo(find(isnan(opt_floo)==1)) = 0;

for yr = 2015 : 2099    %2100
    yr_pct = yr - 2015 + 1; % the year for PCT_xxx data
    yr_opt = yr - 2015 + 1; % the year for option data
    yr_fra = floor((yr-2010)/ 5) + 1; % the year for 5-year data
    yr_int = mod(yr-2010, 5) + 1; % the year for interpolation
    yr_ind = yr-2015+1; % the year for output data
    
    
    % get the optimal fraction
    floo_opt = opt_floo(:, :, yr_opt);
    spri_opt = opt_spri(:, :, yr_opt);
    drip_opt = opt_drip(:, :, yr_opt);
    

    % start get the actual fraction
    % in grid_frac, 1=flood, 2=sprinkler and 3=drip
    drip_fra_be = grid_frac(:,:,yr_fra, 3);
    drip_fra_be(find(isnan(drip_fra_be)==1)) = 0.0000;

    spri_fra_be = grid_frac(:,:,yr_fra, 2);
    spri_fra_be(find(isnan(spri_fra_be)==1)) = 0.0000;

    floo_fra_be = 1 - spri_fra_be - drip_fra_be;
    

    
    drip_fra_af = grid_frac(:,:,yr_fra+1, 3);
    drip_fra_af(find(isnan(drip_fra_af)==1)) = 0.0000;

    spri_fra_af = grid_frac(:,:,yr_fra+1, 2);
    spri_fra_af(find(isnan(spri_fra_af)==1)) = 0.0000;

    floo_fra_af = 1.0000 - spri_fra_af - drip_fra_af;
    
    
    drip_fra_ac = (drip_fra_af - drip_fra_be) / 5 * yr_int + drip_fra_be;
    spri_fra_ac = (spri_fra_af - spri_fra_be) / 5 * yr_int + spri_fra_be;
    floo_fra_ac = zeros(720, 360);

    drip_fra_ac = roundn(drip_fra_ac, -4);
    spri_fra_ac = roundn(spri_fra_ac, -4);

    for i = 1 : 720
        for j = 1 : 360

            if pct_irr(i, j, yr_pct) == 0
                continue
            else
                % check if the fraction exceed the optimal
                if drip_fra_ac(i, j) > drip_opt(i, j)

                    drip_fra_ac(i, j) = drip_opt(i, j);

                end
                if spri_fra_ac(i, j) + drip_fra_ac(i, j) > spri_opt(i, j) + drip_opt(i, j)

                    spri_fra_ac(i, j) = spri_opt(i, j) + drip_opt(i, j) - drip_fra_ac(i, j);
                end

                if spri_fra_ac(i, j) + drip_fra_ac(i, j) > 1.0000
                    spri_fra_ac(i, j) = 1.0000 - drip_fra_ac(i, j);
                end
                
                floo_fra_ac(i, j) = 1.0000 - drip_fra_ac(i, j) - spri_fra_ac(i, j);
            end
        end
    end
    % end get the actual fraction
    drip_fra_ac = roundn(drip_fra_ac, -4);
    spri_fra_ac = roundn(spri_fra_ac, -4);
    floo_fra_ac = roundn(floo_fra_ac, -4);
    
    drip_frac_one_year(:,:,yr_ind) = drip_fra_ac;
    spri_frac_one_year(:,:,yr_ind) = spri_fra_ac;
    floo_frac_one_year(:,:,yr_ind) = floo_fra_ac;
end
% after the interpolation, now distribute the increase in drip and
% sprinkler in 



% the last year
drip_frac_one_year(:,:,86) = grid_frac(:,:,19, 3);

% drip_frac_one_year(find(isnan(drip_frac_one_year)==1)) = 0;

spri_frac_one_year(:,:,86) = grid_frac(:,:,19, 2);
% spri_frac_one_year(find(isnan(spri_frac_one_year)==1)) = 0;

for i = 1 : 720
    for j = 1 : 360
        if isnan(drip_frac_one_year(:,:,86))
            drip_frac_one_year(:,:,86) = drip_frac_one_year(:,:,85)
        end
        if isnan(spri_frac_one_year(:,:,86))
            spri_frac_one_year(:,:,86) = spri_frac_one_year(:,:,85)
        end
    end
end

drip_frac_one_year = roundn(drip_frac_one_year, -4);
spri_frac_one_year = roundn(spri_frac_one_year, -4);


floo_frac_one_year(:,:,86) = 1.0000 - drip_frac_one_year(:,:,86) - spri_frac_one_year(:,:,86);

floo_frac_one_year(floo_frac_one_year < 0) = 0;

fprintf('Interpolation is done. Oh yeah!\n');

for i = 1:720
    for j = 1:360
        for k = 1:86
            if drip_frac_one_year(i,j,k) + spri_frac_one_year(i,j,k) + floo_frac_one_year(i,j,k)==0
                row_left = i - 1;
                row_right = i + 1;
                col_up = j + 1;
                col_bot = j - 1;
                if row_left == 0
                    row_left = 720;
                end
                if row_right == 721
                    row_right = 1;
                end
                if col_up == 361
                    col_up = 1;
                end
                if col_bot == 0
                    col_bot = 360;
                end
                if drip_frac_one_year(row_left,j,k) + spri_frac_one_year(row_left,j,k) + floo_frac_one_year(row_left,j,k) ~= 0
                    drip_frac_one_year(i,j,k) = drip_frac_one_year(row_left,j,k);
                    spri_frac_one_year(i,j,k) = spri_frac_one_year(row_left,j,k);
                    floo_frac_one_year(i,j,k) = floo_frac_one_year(row_left,j,k);

                else
                    
                    if drip_frac_one_year(row_right,j,k) + spri_frac_one_year(row_right,j,k) + floo_frac_one_year(row_right,j,k) ~= 0
                        drip_frac_one_year(i,j,k) = drip_frac_one_year(row_right,j,k);
                        spri_frac_one_year(i,j,k) = spri_frac_one_year(row_right,j,k);
                        floo_frac_one_year(i,j,k) = floo_frac_one_year(row_right,j,k);
                    else
                        if drip_frac_one_year(i,col_up,k) + spri_frac_one_year(i,col_up,k) + floo_frac_one_year(i,col_up,k) ~= 0
                            drip_frac_one_year(i,j,k) = drip_frac_one_year(i,col_up,k);
                            spri_frac_one_year(i,j,k) = spri_frac_one_year(i,col_up,k);
                            floo_frac_one_year(i,j,k) = floo_frac_one_year(i,col_up,k);
                        else
                            if drip_frac_one_year(i,col_bot,k) + spri_frac_one_year(i,col_bot,k) + floo_frac_one_year(i,col_bot,k) ~= 0
                                drip_frac_one_year(i,j,k) = drip_frac_one_year(i,col_bot,k);
                                spri_frac_one_year(i,j,k) = spri_frac_one_year(i,col_bot,k);
                                floo_frac_one_year(i,j,k) = floo_frac_one_year(i,col_bot,k);
                            else
                                if drip_frac_one_year(row_left,col_bot,k) + spri_frac_one_year(row_left,col_bot,k) + floo_frac_one_year(row_left,col_bot,k) ~= 0
                                    drip_frac_one_year(i,j,k) = drip_frac_one_year(row_left,col_bot,k);
                                    spri_frac_one_year(i,j,k) = spri_frac_one_year(row_left,col_bot,k);
                                    floo_frac_one_year(i,j,k) = floo_frac_one_year(row_left,col_bot,k);
                                else
                                    if drip_frac_one_year(row_right,col_bot,k) + spri_frac_one_year(row_right,col_bot,k) + floo_frac_one_year(row_right,col_bot,k) ~= 0
                                        drip_frac_one_year(i,j,k) = drip_frac_one_year(row_right,col_bot,k);
                                        spri_frac_one_year(i,j,k) = spri_frac_one_year(row_right,col_bot,k);
                                        floo_frac_one_year(i,j,k) = floo_frac_one_year(row_right,col_bot,k);
                                    else
                                        if drip_frac_one_year(row_right,col_up,k) + spri_frac_one_year(row_right,col_up,k) + floo_frac_one_year(row_right,col_up,k) ~= 0
                                            drip_frac_one_year(i,j,k) = drip_frac_one_year(row_right,col_up,k);
                                            spri_frac_one_year(i,j,k) = spri_frac_one_year(row_right,col_up,k);
                                            floo_frac_one_year(i,j,k) = floo_frac_one_year(row_right,col_up,k);
                                        else
                                            if drip_frac_one_year(row_left,col_up,k) + spri_frac_one_year(row_left,col_up,k) + floo_frac_one_year(row_left,col_up,k) ~= 0
                                                drip_frac_one_year(i,j,k) = drip_frac_one_year(row_left,col_up,k);
                                                spri_frac_one_year(i,j,k) = spri_frac_one_year(row_left,col_up,k);
                                                floo_frac_one_year(i,j,k) = floo_frac_one_year(row_left,col_up,k);
                                            end
                                        end
                                    end
                                end
                            end
                            
                        end
                    end
                end
            end
            

            if drip_frac_one_year(i,j,k) > opt_drip(i,j,k)

                drip_frac_one_year(i,j,k) = opt_drip(i,j,k);

            end
            if spri_frac_one_year(i,j,k) + drip_frac_one_year(i,j,k) > opt_spri(i, j,k) + opt_drip(i, j,k)

                spri_frac_one_year(i,j,k) = opt_spri(i, j,k) + opt_drip(i, j,k) - drip_frac_one_year(i,j,k);
            end

            if spri_frac_one_year(i,j,k) + drip_frac_one_year(i,j,k) > 1.0000
                spri_frac_one_year(i,j,k) = 1.0000 - drip_frac_one_year(i,j,k);
            end
            
            floo_frac_one_year(i,j,k) = 1.0000 - drip_frac_one_year(i,j,k) - spri_frac_one_year(i,j,k);
        end
    end
end



file_surf = str_data + '.nc';

temperate_cereals_irrigated = ncread(file_surf, 'temperate_cereals_irrigated'); %sprinkler
temperate_cereals_drp_irrigated = zeros(720, 360, 86);
temperate_cereals_spk_irrigated = zeros(720, 360, 86);
temperate_cereals_fld_irrigated = zeros(720, 360, 86);

rice_irrigated = ncread(file_surf, 'rice_irrigated'); %flood
rice_drp_irrigated = zeros(720, 360, 86);
rice_spk_irrigated = zeros(720, 360, 86);
rice_fld_irrigated = zeros(720, 360, 86);

maize_irrigated = ncread(file_surf, 'maize_irrigated'); %sprinkler
maize_drp_irrigated = zeros(720, 360, 86);
maize_spk_irrigated = zeros(720, 360, 86);
maize_fld_irrigated = zeros(720, 360, 86);

tropical_cereals_irrigated = ncread(file_surf, 'tropical_cereals_irrigated'); %sprinkler
tropical_cereals_drp_irrigated = zeros(720, 360, 86);
tropical_cereals_spk_irrigated = zeros(720, 360, 86);
tropical_cereals_fld_irrigated = zeros(720, 360, 86);

pulses_irrigated = ncread(file_surf, 'pulses_irrigated'); %drip
pulses_drp_irrigated = zeros(720, 360, 86);
pulses_spk_irrigated = zeros(720, 360, 86);
pulses_fld_irrigated = zeros(720, 360, 86);

temperate_roots_irrigated = ncread(file_surf, 'temperate_roots_irrigated');%sprinkler
temperate_roots_drp_irrigated = zeros(720, 360, 86);
temperate_roots_spk_irrigated = zeros(720, 360, 86);
temperate_roots_fld_irrigated = zeros(720, 360, 86);

tropical_roots_irrigated = ncread(file_surf, 'tropical_roots_irrigated');%sprinkler
tropical_roots_drp_irrigated = zeros(720, 360, 86);
tropical_roots_spk_irrigated = zeros(720, 360, 86);
tropical_roots_fld_irrigated = zeros(720, 360, 86);

oil_crops_sunflower_irrigated = ncread(file_surf, 'oil_crops_sunflower_irrigated');%drip
oil_crops_sunflower_drp_irrigated = zeros(720, 360, 86);
oil_crops_sunflower_spk_irrigated = zeros(720, 360, 86);
oil_crops_sunflower_fld_irrigated = zeros(720, 360, 86);

oil_crops_soybean_irrigated = ncread(file_surf, 'oil_crops_soybean_irrigated');%drip
oil_crops_soybean_drp_irrigated = zeros(720, 360, 86);
oil_crops_soybean_spk_irrigated = zeros(720, 360, 86);
oil_crops_soybean_fld_irrigated = zeros(720, 360, 86);

oil_crops_groundnut_irrigated = ncread(file_surf, 'oil_crops_groundnut_irrigated');%sprinkler
oil_crops_groundnut_drp_irrigated = zeros(720, 360, 86);
oil_crops_groundnut_spk_irrigated = zeros(720, 360, 86);
oil_crops_groundnut_fld_irrigated = zeros(720, 360, 86);

oil_crops_rapeseed_irrigated = ncread(file_surf, 'oil_crops_rapeseed_irrigated');%sprinkler
oil_crops_rapeseed_drp_irrigated = zeros(720, 360, 86);
oil_crops_rapeseed_spk_irrigated = zeros(720, 360, 86);
oil_crops_rapeseed_fld_irrigated = zeros(720, 360, 86);

c4per_irrigated_food = ncread(file_surf, 'c4per_irrigated_food');%sprinkler
c4per_drp_irrigated_food = zeros(720, 360, 86);
c4per_spk_irrigated_food = zeros(720, 360, 86);
c4per_fld_irrigated_food = zeros(720, 360, 86);

others_c3ann_irrigated = ncread(file_surf, 'others_c3ann_irrigated');%drip
others_c3ann_drp_irrigated = zeros(720, 360, 86);
others_c3ann_spk_irrigated = zeros(720, 360, 86);
others_c3ann_fld_irrigated = zeros(720, 360, 86);

others_c3nfx_irrigated = ncread(file_surf, 'others_c3nfx_irrigated');%drip
others_c3nfx_drp_irrigated = zeros(720, 360, 86);
others_c3nfx_spk_irrigated = zeros(720, 360, 86);
others_c3nfx_fld_irrigated = zeros(720, 360, 86);

c3per_irrigated_food = ncread(file_surf, 'c3per_irrigated_food');%sprinkler
c3per_drp_irrigated_food = zeros(720, 360, 86);
c3per_spk_irrigated_food = zeros(720, 360, 86);
c3per_fld_irrigated_food = zeros(720, 360, 86);

c3per_irrigated_bf = ncread(file_surf, 'c3per_irrigated_bf');%sprinkler
c3per_drp_irrigated_bf = zeros(720, 360, 86);
c3per_spk_irrigated_bf = zeros(720, 360, 86);
c3per_fld_irrigated_bf = zeros(720, 360, 86);

c4per_irrigated_bf = ncread(file_surf, 'c4per_irrigated_bf');%sprinkler
c4per_drp_irrigated_bf = zeros(720, 360, 86);
c4per_spk_irrigated_bf = zeros(720, 360, 86);
c4per_fld_irrigated_bf = zeros(720, 360, 86);



for yr = 2015 : 2100    %2100
    fprintf('year= %d', yr)
    yr_pct = yr - 2015 + 1; % the year for PCT_xxx data
    yr_opt = yr - 2015 + 1; % the year for option data
    yr_fra = floor((yr-2010)/ 5) + 1; % the year for 5-year data
    yr_int = mod(yr-2010, 5) + 1; % the year for interpolation
    yr_ind = yr - 2015 + 1; % the year for output data

    for i = 1: 720 %288
        for j = 1 : 360 %192
            if pct_irr(i, j, yr_pct) == 0 || isnan(pct_irr(i, j, yr_pct))
%                 fprintf(['pixel has no pct_irr \n', num2str(i), ' ', num2str(j)])
                continue
            else
                spri_all = pct_irr(i, j, yr_pct) .* spri_frac_one_year(i, j, yr_ind);
%                 fprintf(['spri_all = ', num2str(spri_all), '\n'])
                drip_all = pct_irr(i, j, yr_pct) .* drip_frac_one_year(i, j, yr_ind);
%                 fprintf(['drip_all = ', num2str(drip_all), '\n'])
                floo_all = pct_irr(i, j, yr_pct) .* floo_frac_one_year(i, j, yr_ind);
%                 fprintf(['floo_all = ', num2str(floo_all), '\n'])
                if opt_drip(i, j, yr_ind) == 0 || isnan(opt_drip(i, j, yr_ind))
                    if drip_all ~= 0
                        fprintf('something wrong with drip')
                    end
                    continue
                else
%                     fprintf('drip_frac_one_year = %f \n', drip_frac_one_year(i, j, yr_ind))
%                     fprintf('drip_all = %f \n', drip_all)
                    pulses_drp_irrigated(i, j, yr_ind) = (pulses_irrigated(i, j, yr_ind) ./ opt_drip(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * drip_all;
                    oil_crops_sunflower_drp_irrigated(i, j, yr_ind) = (oil_crops_sunflower_irrigated(i, j, yr_ind) ./ opt_drip(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * drip_all;
                    oil_crops_soybean_drp_irrigated(i, j, yr_ind) = (oil_crops_soybean_irrigated(i, j, yr_ind) ./ opt_drip(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * drip_all;
                    others_c3ann_drp_irrigated(i, j, yr_ind) = (others_c3ann_irrigated(i, j, yr_ind) ./ opt_drip(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * drip_all;
                    others_c3nfx_drp_irrigated(i, j, yr_ind) = (others_c3nfx_irrigated(i, j, yr_ind) ./ opt_drip(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * drip_all;
                end
                
                if opt_spri(i, j, yr_ind) == 0 || isnan(opt_spri(i, j, yr_ind))
                    if spri_all ~= 0
                        fprintf('something wrong with spri')
                    end
                    continue
                else
                    if spri_frac_one_year(i, j, yr_ind) <= opt_spri(i, j, yr_ind)
                        temperate_cereals_spk_irrigated(i, j, yr_ind) = (temperate_cereals_irrigated(i, j, yr_ind) ./ opt_spri(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_all;
                        maize_spk_irrigated(i, j, yr_ind) = (maize_irrigated(i, j, yr_ind) ./ opt_spri(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_all;
                        tropical_cereals_spk_irrigated(i, j, yr_ind) = (tropical_cereals_irrigated(i, j, yr_ind) ./ opt_spri(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_all;
                        temperate_roots_spk_irrigated(i, j, yr_ind) = (temperate_roots_irrigated(i, j, yr_ind) ./ opt_spri(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_all;
                        tropical_roots_spk_irrigated(i, j, yr_ind) = (tropical_roots_irrigated(i, j, yr_ind) ./ opt_spri(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_all;
                        oil_crops_groundnut_spk_irrigated(i, j, yr_ind) = (oil_crops_groundnut_irrigated(i, j, yr_ind) ./ opt_spri(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_all;
                        oil_crops_rapeseed_spk_irrigated(i, j, yr_ind) = (oil_crops_rapeseed_irrigated(i, j, yr_ind) ./ opt_spri(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_all;
                        c4per_spk_irrigated_food(i, j, yr_ind) = (c4per_irrigated_food(i, j, yr_ind) ./ opt_spri(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_all;
                        c3per_spk_irrigated_food(i, j, yr_ind) = (c3per_irrigated_food(i, j, yr_ind) ./ opt_spri(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_all;
                        c3per_spk_irrigated_bf(i, j, yr_ind) = (c3per_irrigated_bf(i, j, yr_ind) ./ opt_spri(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_all;
                        c4per_spk_irrigated_bf(i, j, yr_ind) = (c4per_irrigated_bf(i, j, yr_ind) ./ opt_spri(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_all;
                    else
                        temperate_cereals_spk_irrigated(i, j, yr_ind) = temperate_cereals_irrigated(i, j, yr_ind);
                        maize_spk_irrigated(i, j, yr_ind) = maize_irrigated(i, j, yr_ind);
                        tropical_cereals_spk_irrigated(i, j, yr_ind) = tropical_cereals_irrigated(i, j, yr_ind);
                        temperate_roots_spk_irrigated(i, j, yr_ind) = temperate_roots_irrigated(i, j, yr_ind);
                        tropical_roots_spk_irrigated(i, j, yr_ind) = tropical_roots_irrigated(i, j, yr_ind);
                        oil_crops_groundnut_spk_irrigated(i, j, yr_ind) = oil_crops_groundnut_irrigated(i, j, yr_ind);
                        oil_crops_rapeseed_spk_irrigated(i, j, yr_ind) = oil_crops_rapeseed_irrigated(i, j, yr_ind);
                        c4per_spk_irrigated_food(i, j, yr_ind) = c4per_irrigated_food(i, j, yr_ind);
                        c3per_spk_irrigated_food(i, j, yr_ind) = c3per_irrigated_food(i, j, yr_ind);
                        c3per_spk_irrigated_bf(i, j, yr_ind) = c3per_irrigated_bf(i, j, yr_ind);
                        c4per_spk_irrigated_bf(i, j, yr_ind) = c4per_irrigated_bf(i, j, yr_ind);
    
                        spri_excess = spri_all - opt_spri(i, j, yr_ind) * pct_irr(i, j, yr_pct);
                        pulses_spk_irrigated(i, j, yr_ind) = (pulses_irrigated(i, j, yr_ind) ./ opt_drip(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_excess;
                        oil_crops_sunflower_spk_irrigated(i, j, yr_ind) = (oil_crops_sunflower_irrigated(i, j, yr_ind) ./ opt_drip(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_excess;
                        oil_crops_soybean_spk_irrigated(i, j, yr_ind) = (oil_crops_soybean_irrigated(i, j, yr_ind) ./ opt_drip(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_excess;
                        others_c3ann_spk_irrigated(i, j, yr_ind) = (others_c3ann_irrigated(i, j, yr_ind) ./ opt_drip(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_excess;
                        others_c3nfx_spk_irrigated(i, j, yr_ind) = (others_c3nfx_irrigated(i, j, yr_ind) ./ opt_drip(i, j, yr_ind) ./ pct_irr(i, j, yr_pct)) * spri_excess;
                    end
                end

                rice_fld_irrigated(i, j, yr_ind) = rice_irrigated(i, j, yr_ind);

                pulses_fld_irrigated(i, j, yr_ind) = pulses_irrigated(i, j, yr_ind) - pulses_drp_irrigated(i, j, yr_ind) - pulses_spk_irrigated(i, j, yr_ind);
                oil_crops_sunflower_fld_irrigated(i, j, yr_ind) = oil_crops_sunflower_irrigated(i, j, yr_ind) - oil_crops_sunflower_drp_irrigated(i, j, yr_ind) - oil_crops_sunflower_spk_irrigated(i, j, yr_ind);
                oil_crops_soybean_fld_irrigated(i, j, yr_ind) = oil_crops_soybean_irrigated(i, j, yr_ind) - oil_crops_soybean_drp_irrigated(i, j, yr_ind) - oil_crops_soybean_spk_irrigated(i, j, yr_ind);
                others_c3ann_fld_irrigated(i, j, yr_ind) = others_c3ann_irrigated(i, j, yr_ind) - others_c3ann_drp_irrigated(i, j, yr_ind) - others_c3ann_spk_irrigated(i, j, yr_ind);
                others_c3nfx_fld_irrigated(i, j, yr_ind) = others_c3nfx_irrigated(i, j, yr_ind) - others_c3nfx_drp_irrigated(i, j, yr_ind) - others_c3nfx_spk_irrigated(i, j, yr_ind);

                temperate_cereals_fld_irrigated(i, j, yr_ind) = temperate_cereals_irrigated(i, j, yr_ind) - temperate_cereals_spk_irrigated(i, j, yr_ind);
                maize_fld_irrigated(i, j, yr_ind) = maize_irrigated(i, j, yr_ind) - maize_spk_irrigated(i, j, yr_ind);
                tropical_cereals_fld_irrigated(i, j, yr_ind) = tropical_cereals_irrigated(i, j, yr_ind) - tropical_cereals_spk_irrigated(i, j, yr_ind);
                temperate_roots_fld_irrigated(i, j, yr_ind) = temperate_roots_irrigated(i, j, yr_ind) - temperate_roots_spk_irrigated(i, j, yr_ind);
                tropical_roots_fld_irrigated(i, j, yr_ind) = tropical_roots_irrigated(i, j, yr_ind) - tropical_roots_spk_irrigated(i, j, yr_ind);
                oil_crops_groundnut_fld_irrigated(i, j, yr_ind) = oil_crops_groundnut_irrigated(i, j, yr_ind) - oil_crops_groundnut_spk_irrigated(i, j, yr_ind);
                oil_crops_rapeseed_fld_irrigated(i, j, yr_ind) = oil_crops_rapeseed_irrigated(i, j, yr_ind) - oil_crops_rapeseed_spk_irrigated(i, j, yr_ind);
                c4per_fld_irrigated_food(i, j, yr_ind) = c4per_irrigated_food(i, j, yr_ind) - c4per_spk_irrigated_food(i, j, yr_ind);
                c3per_fld_irrigated_food(i, j, yr_ind) = c3per_irrigated_food(i, j, yr_ind) - c3per_spk_irrigated_food(i, j, yr_ind);
                c3per_fld_irrigated_bf(i, j, yr_ind) = c3per_irrigated_bf(i, j, yr_ind) - c3per_spk_irrigated_bf(i, j, yr_ind);
                c4per_fld_irrigated_bf(i, j, yr_ind) = c4per_irrigated_bf(i, j, yr_ind) - c4per_spk_irrigated_bf(i, j, yr_ind);
            end
            pulses_fld_irrigated(i, j, yr_ind) = pulses_irrigated(i, j, yr_ind) - pulses_drp_irrigated(i, j, yr_ind) - pulses_spk_irrigated(i, j, yr_ind);
            oil_crops_sunflower_fld_irrigated(i, j, yr_ind) = oil_crops_sunflower_irrigated(i, j, yr_ind) - oil_crops_sunflower_drp_irrigated(i, j, yr_ind) - oil_crops_sunflower_spk_irrigated(i, j, yr_ind);
            oil_crops_soybean_fld_irrigated(i, j, yr_ind) = oil_crops_soybean_irrigated(i, j, yr_ind) - oil_crops_soybean_drp_irrigated(i, j, yr_ind) - oil_crops_soybean_spk_irrigated(i, j, yr_ind);
            others_c3ann_fld_irrigated(i, j, yr_ind) = others_c3ann_irrigated(i, j, yr_ind) - others_c3ann_drp_irrigated(i, j, yr_ind) - others_c3ann_spk_irrigated(i, j, yr_ind);
            others_c3nfx_fld_irrigated(i, j, yr_ind) = others_c3nfx_irrigated(i, j, yr_ind) - others_c3nfx_drp_irrigated(i, j, yr_ind) - others_c3nfx_spk_irrigated(i, j, yr_ind);

            temperate_cereals_fld_irrigated(i, j, yr_ind) = temperate_cereals_irrigated(i, j, yr_ind) - temperate_cereals_spk_irrigated(i, j, yr_ind);
            maize_fld_irrigated(i, j, yr_ind) = maize_irrigated(i, j, yr_ind) - maize_spk_irrigated(i, j, yr_ind);
            tropical_cereals_fld_irrigated(i, j, yr_ind) = tropical_cereals_irrigated(i, j, yr_ind) - tropical_cereals_spk_irrigated(i, j, yr_ind);
            temperate_roots_fld_irrigated(i, j, yr_ind) = temperate_roots_irrigated(i, j, yr_ind) - temperate_roots_spk_irrigated(i, j, yr_ind);
            tropical_roots_fld_irrigated(i, j, yr_ind) = tropical_roots_irrigated(i, j, yr_ind) - tropical_roots_spk_irrigated(i, j, yr_ind);
            oil_crops_groundnut_fld_irrigated(i, j, yr_ind) = oil_crops_groundnut_irrigated(i, j, yr_ind) - oil_crops_groundnut_spk_irrigated(i, j, yr_ind);
            oil_crops_rapeseed_fld_irrigated(i, j, yr_ind) = oil_crops_rapeseed_irrigated(i, j, yr_ind) - oil_crops_rapeseed_spk_irrigated(i, j, yr_ind);
            c4per_fld_irrigated_food(i, j, yr_ind) = c4per_irrigated_food(i, j, yr_ind) - c4per_spk_irrigated_food(i, j, yr_ind);
            c3per_fld_irrigated_food(i, j, yr_ind) = c3per_irrigated_food(i, j, yr_ind) - c3per_spk_irrigated_food(i, j, yr_ind);
            c3per_fld_irrigated_bf(i, j, yr_ind) = c3per_irrigated_bf(i, j, yr_ind) - c3per_spk_irrigated_bf(i, j, yr_ind);
            c4per_fld_irrigated_bf(i, j, yr_ind) = c4per_irrigated_bf(i, j, yr_ind) - c4per_spk_irrigated_bf(i, j, yr_ind);
        end
    end
end

for i = 1 : 720
    for j = 1 : 360
        for yr_ind = 1 : 86
            pulses_fld_irrigated(i, j, yr_ind) = pulses_irrigated(i, j, yr_ind) - pulses_drp_irrigated(i, j, yr_ind) - pulses_spk_irrigated(i, j, yr_ind);
            oil_crops_sunflower_fld_irrigated(i, j, yr_ind) = oil_crops_sunflower_irrigated(i, j, yr_ind) - oil_crops_sunflower_drp_irrigated(i, j, yr_ind) - oil_crops_sunflower_spk_irrigated(i, j, yr_ind);
            oil_crops_soybean_fld_irrigated(i, j, yr_ind) = oil_crops_soybean_irrigated(i, j, yr_ind) - oil_crops_soybean_drp_irrigated(i, j, yr_ind) - oil_crops_soybean_spk_irrigated(i, j, yr_ind);
            others_c3ann_fld_irrigated(i, j, yr_ind) = others_c3ann_irrigated(i, j, yr_ind) - others_c3ann_drp_irrigated(i, j, yr_ind) - others_c3ann_spk_irrigated(i, j, yr_ind);
            others_c3nfx_fld_irrigated(i, j, yr_ind) = others_c3nfx_irrigated(i, j, yr_ind) - others_c3nfx_drp_irrigated(i, j, yr_ind) - others_c3nfx_spk_irrigated(i, j, yr_ind);

            temperate_cereals_fld_irrigated(i, j, yr_ind) = temperate_cereals_irrigated(i, j, yr_ind) - temperate_cereals_spk_irrigated(i, j, yr_ind);
            maize_fld_irrigated(i, j, yr_ind) = maize_irrigated(i, j, yr_ind) - maize_spk_irrigated(i, j, yr_ind);
            tropical_cereals_fld_irrigated(i, j, yr_ind) = tropical_cereals_irrigated(i, j, yr_ind) - tropical_cereals_spk_irrigated(i, j, yr_ind);
            temperate_roots_fld_irrigated(i, j, yr_ind) = temperate_roots_irrigated(i, j, yr_ind) - temperate_roots_spk_irrigated(i, j, yr_ind);
            tropical_roots_fld_irrigated(i, j, yr_ind) = tropical_roots_irrigated(i, j, yr_ind) - tropical_roots_spk_irrigated(i, j, yr_ind);
            oil_crops_groundnut_fld_irrigated(i, j, yr_ind) = oil_crops_groundnut_irrigated(i, j, yr_ind) - oil_crops_groundnut_spk_irrigated(i, j, yr_ind);
            oil_crops_rapeseed_fld_irrigated(i, j, yr_ind) = oil_crops_rapeseed_irrigated(i, j, yr_ind) - oil_crops_rapeseed_spk_irrigated(i, j, yr_ind);
            c4per_fld_irrigated_food(i, j, yr_ind) = c4per_irrigated_food(i, j, yr_ind) - c4per_spk_irrigated_food(i, j, yr_ind);
            c3per_fld_irrigated_food(i, j, yr_ind) = c3per_irrigated_food(i, j, yr_ind) - c3per_spk_irrigated_food(i, j, yr_ind);
            c3per_fld_irrigated_bf(i, j, yr_ind) = c3per_irrigated_bf(i, j, yr_ind) - c3per_spk_irrigated_bf(i, j, yr_ind);
            c4per_fld_irrigated_bf(i, j, yr_ind) = c4per_irrigated_bf(i, j, yr_ind) - c4per_spk_irrigated_bf(i, j, yr_ind);
        end
    end
end


ncid = netcdf.open(file_surf, 'WRITE');
lonDimid = netcdf.inqDimID(ncid, "lon");
latDimid = netcdf.inqDimID(ncid, "lat");
timeDimid = netcdf.inqDimID(ncid, "time");

netcdf.reDef(ncid);
temperate_cereals_drp_irrigated_id = netcdf.defVar(ncid, 'temperate_cereals_drp_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, temperate_cereals_drp_irrigated_id, temperate_cereals_drp_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, temperate_cereals_drp_irrigated_id, 'longname', 'drip irrigated temperated cereals')
netcdf.putAtt(ncid, temperate_cereals_drp_irrigated_id, 'units', '1')
netcdf.endDef(ncid);

netcdf.reDef(ncid);
temperate_cereals_spk_irrigated_id = netcdf.defVar(ncid, 'temperate_cereals_spk_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, temperate_cereals_spk_irrigated_id, temperate_cereals_spk_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, temperate_cereals_spk_irrigated_id, 'longname', 'sprinkler irrigated temperated cereals')
netcdf.putAtt(ncid, temperate_cereals_spk_irrigated_id, 'units', '1')
netcdf.endDef(ncid);


netcdf.reDef(ncid);
temperate_cereals_fld_irrigated_id = netcdf.defVar(ncid, 'temperate_cereals_fld_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, temperate_cereals_fld_irrigated_id, temperate_cereals_fld_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, temperate_cereals_fld_irrigated_id, 'longname', 'flood irrigated temperated cereals')
netcdf.putAtt(ncid, temperate_cereals_fld_irrigated_id, 'units', '1')
netcdf.endDef(ncid);


netcdf.reDef(ncid);
rice_drp_irrigated_id = netcdf.defVar(ncid, 'rice_drp_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, rice_drp_irrigated_id, rice_drp_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, rice_drp_irrigated_id, 'longname', 'drip irrigated rice')
netcdf.putAtt(ncid, rice_drp_irrigated_id, 'units', '1')
netcdf.endDef(ncid);

netcdf.reDef(ncid);
rice_spk_irrigated_id = netcdf.defVar(ncid, 'rice_spk_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, rice_spk_irrigated_id, rice_spk_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, rice_spk_irrigated_id, 'longname', 'sprinkler irrigated rice')
netcdf.putAtt(ncid, rice_spk_irrigated_id, 'units', '1')
netcdf.endDef(ncid);

netcdf.reDef(ncid);
rice_fld_irrigated_id = netcdf.defVar(ncid, 'rice_fld_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, rice_fld_irrigated_id, rice_fld_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, rice_fld_irrigated_id, 'longname', 'flood irrigated rice')
netcdf.putAtt(ncid, rice_fld_irrigated_id, 'units', '1')
netcdf.endDef(ncid);

netcdf.reDef(ncid);
maize_drp_irrigated_id = netcdf.defVar(ncid, 'maize_drp_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, maize_drp_irrigated_id, maize_drp_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, maize_drp_irrigated_id, 'longname', 'drip irrigated maize')
netcdf.putAtt(ncid, maize_drp_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
maize_spk_irrigated_id = netcdf.defVar(ncid, 'maize_spk_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, maize_spk_irrigated_id, maize_spk_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, maize_spk_irrigated_id, 'longname', 'sprinkler irrigated maize')
netcdf.putAtt(ncid, maize_spk_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
maize_fld_irrigated_id = netcdf.defVar(ncid, 'maize_fld_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, maize_fld_irrigated_id, maize_fld_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, maize_fld_irrigated_id, 'longname', 'flood irrigated maize')
netcdf.putAtt(ncid, maize_fld_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
tropical_cereals_drp_irrigated_id = netcdf.defVar(ncid, 'tropical_cereals_drp_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, tropical_cereals_drp_irrigated_id, tropical_cereals_drp_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, tropical_cereals_drp_irrigated_id, 'longname', 'drip irrigated tropical cereals')
netcdf.putAtt(ncid, tropical_cereals_drp_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
tropical_cereals_spk_irrigated_id = netcdf.defVar(ncid, 'tropical_cereals_spk_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, tropical_cereals_spk_irrigated_id, tropical_cereals_spk_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, tropical_cereals_spk_irrigated_id, 'longname', 'sprinkler irrigated tropical cereals')
netcdf.putAtt(ncid, tropical_cereals_spk_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
tropical_cereals_fld_irrigated_id = netcdf.defVar(ncid, 'tropical_cereals_fld_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, tropical_cereals_fld_irrigated_id, tropical_cereals_fld_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, tropical_cereals_fld_irrigated_id, 'longname', 'flood irrigated tropical cereals')
netcdf.putAtt(ncid, tropical_cereals_fld_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
pulses_drp_irrigated_id = netcdf.defVar(ncid, 'pulses_drp_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, pulses_drp_irrigated_id, pulses_drp_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, pulses_drp_irrigated_id, 'longname', 'drip irrigated pulses')
netcdf.putAtt(ncid, pulses_drp_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
pulses_spk_irrigated_id = netcdf.defVar(ncid, 'pulses_spk_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, pulses_spk_irrigated_id, pulses_spk_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, pulses_spk_irrigated_id, 'longname', 'sprinkler irrigated pulses')
netcdf.putAtt(ncid, pulses_spk_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
pulses_fld_irrigated_id = netcdf.defVar(ncid, 'pulses_fld_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, pulses_fld_irrigated_id, pulses_fld_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, pulses_fld_irrigated_id, 'longname', 'flood irrigated pulses')
netcdf.putAtt(ncid, pulses_fld_irrigated_id, 'units', '1')
netcdf.endDef(ncid);

netcdf.reDef(ncid);
temperate_roots_drp_irrigated_id = netcdf.defVar(ncid, 'temperate_roots_drp_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, temperate_roots_drp_irrigated_id, temperate_roots_drp_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, temperate_roots_drp_irrigated_id, 'longname', 'drip irrigated temperated roots')
netcdf.putAtt(ncid, temperate_roots_drp_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
temperate_roots_spk_irrigated_id = netcdf.defVar(ncid, 'temperate_roots_spk_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, temperate_roots_spk_irrigated_id, temperate_roots_spk_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, temperate_roots_spk_irrigated_id, 'longname', 'sprinkler irrigated temperated roots')
netcdf.putAtt(ncid, temperate_roots_spk_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
temperate_roots_fld_irrigated_id = netcdf.defVar(ncid, 'temperate_roots_fld_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, temperate_roots_fld_irrigated_id, temperate_roots_fld_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, temperate_roots_fld_irrigated_id, 'longname', 'flood irrigated temperated roots')
netcdf.putAtt(ncid, temperate_roots_fld_irrigated_id, 'units', '1')
netcdf.endDef(ncid);

netcdf.reDef(ncid);
tropical_roots_drp_irrigated_id = netcdf.defVar(ncid, 'tropical_roots_drp_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, tropical_roots_drp_irrigated_id, tropical_roots_drp_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, tropical_roots_drp_irrigated_id, 'longname', 'drip irrigated tropical roots')
netcdf.putAtt(ncid, tropical_roots_drp_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
tropical_roots_spk_irrigated_id = netcdf.defVar(ncid, 'tropical_roots_spk_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, tropical_roots_spk_irrigated_id, tropical_roots_spk_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, tropical_roots_spk_irrigated_id, 'longname', 'sprinkler irrigated tropical roots')
netcdf.putAtt(ncid, tropical_roots_spk_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
tropical_roots_fld_irrigated_id = netcdf.defVar(ncid, 'tropical_roots_fld_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, tropical_roots_fld_irrigated_id, tropical_roots_fld_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, tropical_roots_fld_irrigated_id, 'longname', 'flood irrigated tropical roots')
netcdf.putAtt(ncid, tropical_roots_fld_irrigated_id, 'units', '1')
netcdf.endDef(ncid);

netcdf.reDef(ncid);
oil_crops_sunflower_drp_irrigated_id = netcdf.defVar(ncid, 'oil_crops_sunflower_drp_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, oil_crops_sunflower_drp_irrigated_id, oil_crops_sunflower_drp_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, oil_crops_sunflower_drp_irrigated_id, 'longname', 'drip irrigated oil crops sunflower')
netcdf.putAtt(ncid, oil_crops_sunflower_drp_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
oil_crops_sunflower_spk_irrigated_id = netcdf.defVar(ncid, 'oil_crops_sunflower_spk_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, oil_crops_sunflower_spk_irrigated_id, oil_crops_sunflower_spk_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, oil_crops_sunflower_spk_irrigated_id, 'longname', 'sprinkler irrigated oil crops sunflower')
netcdf.putAtt(ncid, oil_crops_sunflower_spk_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
oil_crops_sunflower_fld_irrigated_id = netcdf.defVar(ncid, 'oil_crops_sunflower_fld_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, oil_crops_sunflower_fld_irrigated_id, oil_crops_sunflower_fld_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, oil_crops_sunflower_fld_irrigated_id, 'longname', 'flood irrigated oil crops sunflower')
netcdf.putAtt(ncid, oil_crops_sunflower_fld_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
oil_crops_soybean_drp_irrigated_id = netcdf.defVar(ncid, 'oil_crops_soybean_drp_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, oil_crops_soybean_drp_irrigated_id, oil_crops_soybean_drp_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, oil_crops_soybean_drp_irrigated_id, 'longname', 'drip irrigated oil crops soybean')
netcdf.putAtt(ncid, oil_crops_soybean_drp_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
oil_crops_soybean_spk_irrigated_id = netcdf.defVar(ncid, 'oil_crops_soybean_spk_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, oil_crops_soybean_spk_irrigated_id, oil_crops_soybean_spk_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, oil_crops_soybean_spk_irrigated_id, 'longname', 'sprinkler irrigated oil crops soybean')
netcdf.putAtt(ncid, oil_crops_soybean_spk_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
oil_crops_soybean_fld_irrigated_id = netcdf.defVar(ncid, 'oil_crops_soybean_fld_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, oil_crops_soybean_fld_irrigated_id, oil_crops_soybean_fld_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, oil_crops_soybean_fld_irrigated_id, 'longname', 'flood irrigated oil crops soybean')
netcdf.putAtt(ncid, oil_crops_soybean_fld_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
oil_crops_groundnut_drp_irrigated_id = netcdf.defVar(ncid, 'oil_crops_groundnut_drp_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, oil_crops_groundnut_drp_irrigated_id, oil_crops_groundnut_drp_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, oil_crops_groundnut_drp_irrigated_id, 'longname', 'drip irrigated oil crops groundnut')
netcdf.putAtt(ncid, oil_crops_groundnut_drp_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
oil_crops_groundnut_spk_irrigated_id = netcdf.defVar(ncid, 'oil_crops_groundnut_spk_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, oil_crops_groundnut_spk_irrigated_id, oil_crops_groundnut_spk_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, oil_crops_groundnut_spk_irrigated_id, 'longname', 'sprinkler irrigated oil crops groundnut')
netcdf.putAtt(ncid, oil_crops_groundnut_spk_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
oil_crops_groundnut_fld_irrigated_id = netcdf.defVar(ncid, 'oil_crops_groundnut_fld_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, oil_crops_groundnut_fld_irrigated_id, oil_crops_groundnut_fld_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, oil_crops_groundnut_fld_irrigated_id, 'longname', 'flood irrigated oil crops groundnut')
netcdf.putAtt(ncid, oil_crops_groundnut_fld_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
oil_crops_rapeseed_drp_irrigated_id = netcdf.defVar(ncid, 'oil_crops_rapeseed_drp_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, oil_crops_rapeseed_drp_irrigated_id, oil_crops_rapeseed_drp_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, oil_crops_rapeseed_drp_irrigated_id, 'longname', 'drip irrigated oil crops rapeseed')
netcdf.putAtt(ncid, oil_crops_rapeseed_drp_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
oil_crops_rapeseed_spk_irrigated_id = netcdf.defVar(ncid, 'oil_crops_rapeseed_spk_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, oil_crops_rapeseed_spk_irrigated_id, oil_crops_rapeseed_spk_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, oil_crops_rapeseed_spk_irrigated_id, 'longname', 'sprinkler irrigated oil crops rapeseed')
netcdf.putAtt(ncid, oil_crops_rapeseed_spk_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
oil_crops_rapeseed_fld_irrigated_id = netcdf.defVar(ncid, 'oil_crops_rapeseed_fld_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, oil_crops_rapeseed_fld_irrigated_id, oil_crops_rapeseed_fld_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, oil_crops_rapeseed_fld_irrigated_id, 'longname', 'flood irrigated oil crops rapeseed')
netcdf.putAtt(ncid, oil_crops_rapeseed_fld_irrigated_id, 'units', '1')
netcdf.endDef(ncid);

netcdf.reDef(ncid);
c4per_drp_irrigated_food_id = netcdf.defVar(ncid, 'c4per_drp_irrigated_food', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, c4per_drp_irrigated_food_id, c4per_drp_irrigated_food);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, c4per_drp_irrigated_food_id, 'longname', 'drip irrigated c4per food')
netcdf.putAtt(ncid, c4per_drp_irrigated_food_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
c4per_spk_irrigated_food_id = netcdf.defVar(ncid, 'c4per_spk_irrigated_food', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, c4per_spk_irrigated_food_id, c4per_spk_irrigated_food);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, c4per_spk_irrigated_food_id, 'longname', 'sprinkler irrigated c4per food')
netcdf.putAtt(ncid, c4per_spk_irrigated_food_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
c4per_fld_irrigated_food_id = netcdf.defVar(ncid, 'c4per_fld_irrigated_food', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, c4per_fld_irrigated_food_id, c4per_fld_irrigated_food);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, c4per_fld_irrigated_food_id, 'longname', 'flood irrigated c4per food')
netcdf.putAtt(ncid, c4per_fld_irrigated_food_id, 'units', '1')
netcdf.endDef(ncid);

netcdf.reDef(ncid);
c4per_drp_irrigated_bf_id = netcdf.defVar(ncid, 'c4per_drp_irrigated_bf', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, c4per_drp_irrigated_bf_id, c4per_drp_irrigated_bf);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, c4per_drp_irrigated_bf_id, 'longname', 'drip irrigated c4per bf')
netcdf.putAtt(ncid, c4per_drp_irrigated_bf_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
c4per_spk_irrigated_bf_id = netcdf.defVar(ncid, 'c4per_spk_irrigated_bf', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, c4per_spk_irrigated_bf_id, c4per_spk_irrigated_bf);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, c4per_spk_irrigated_bf_id, 'longname', 'sprinkler irrigated c4per bf')
netcdf.putAtt(ncid, c4per_spk_irrigated_bf_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
c4per_fld_irrigated_bf_id = netcdf.defVar(ncid, 'c4per_fld_irrigated_bf', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, c4per_fld_irrigated_bf_id, c4per_fld_irrigated_bf);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, c4per_fld_irrigated_bf_id, 'longname', 'flood irrigated c4per bf')
netcdf.putAtt(ncid, c4per_fld_irrigated_bf_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
c3per_drp_irrigated_bf_id = netcdf.defVar(ncid, 'c3per_drp_irrigated_bf', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, c3per_drp_irrigated_bf_id, c3per_drp_irrigated_bf);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, c3per_drp_irrigated_bf_id, 'longname', 'drip irrigated c3per bf')
netcdf.putAtt(ncid, c3per_drp_irrigated_bf_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
c3per_spk_irrigated_bf_id = netcdf.defVar(ncid, 'c3per_spk_irrigated_bf', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, c3per_spk_irrigated_bf_id, c3per_spk_irrigated_bf);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, c3per_spk_irrigated_bf_id, 'longname', 'sprinkler irrigated c3per bfs')
netcdf.putAtt(ncid, c3per_spk_irrigated_bf_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
c3per_fld_irrigated_bf_id = netcdf.defVar(ncid, 'c3per_fld_irrigated_bf', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, c3per_fld_irrigated_bf_id, c3per_fld_irrigated_bf);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, c3per_fld_irrigated_bf_id, 'longname', 'flood irrigated c3per bf')
netcdf.putAtt(ncid, c3per_fld_irrigated_bf_id, 'units', '1')
netcdf.endDef(ncid);

netcdf.reDef(ncid);
c3per_drp_irrigated_food_id = netcdf.defVar(ncid, 'c3per_drp_irrigated_food', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, c3per_drp_irrigated_food_id, c3per_drp_irrigated_food);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, c3per_drp_irrigated_food_id, 'longname', 'drip irrigated c3per food')
netcdf.putAtt(ncid, c3per_drp_irrigated_food_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
c3per_spk_irrigated_food_id = netcdf.defVar(ncid, 'c3per_spk_irrigated_food', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, c3per_spk_irrigated_food_id, c3per_spk_irrigated_food);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, c3per_spk_irrigated_food_id, 'longname', 'sprinkler irrigated c3per food')
netcdf.putAtt(ncid, c3per_spk_irrigated_food_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
c3per_fld_irrigated_food_id = netcdf.defVar(ncid, 'c3per_fld_irrigated_food', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, c3per_fld_irrigated_food_id, c3per_fld_irrigated_food);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, c3per_fld_irrigated_food_id, 'longname', 'flood irrigated c3per food')
netcdf.putAtt(ncid, c3per_fld_irrigated_food_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
others_c3ann_drp_irrigated_id = netcdf.defVar(ncid, 'others_c3ann_drp_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, others_c3ann_drp_irrigated_id, others_c3ann_drp_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, others_c3ann_drp_irrigated_id, 'longname', 'drip irrigated others c3ann')
netcdf.putAtt(ncid, others_c3ann_drp_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
others_c3ann_spk_irrigated_id = netcdf.defVar(ncid, 'others_c3ann_spk_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, others_c3ann_spk_irrigated_id, others_c3ann_spk_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, others_c3ann_spk_irrigated_id, 'longname', 'sprinkler irrigated others c3ann')
netcdf.putAtt(ncid, others_c3ann_spk_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
others_c3ann_fld_irrigated_id = netcdf.defVar(ncid, 'others_c3ann_fld_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, others_c3ann_fld_irrigated_id, others_c3ann_fld_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, others_c3ann_fld_irrigated_id, 'longname', 'flood irrigated others c3ann')
netcdf.putAtt(ncid, others_c3ann_fld_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
others_c3nfx_drp_irrigated_id = netcdf.defVar(ncid, 'others_c3nfx_drp_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, others_c3nfx_drp_irrigated_id, others_c3nfx_drp_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, others_c3nfx_drp_irrigated_id, 'longname', 'drip irrigated others c3nfx')
netcdf.putAtt(ncid, others_c3nfx_drp_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
others_c3nfx_spk_irrigated_id = netcdf.defVar(ncid, 'others_c3nfx_spk_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, others_c3nfx_spk_irrigated_id, others_c3nfx_spk_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, others_c3nfx_spk_irrigated_id, 'longname', 'sprinkler irrigated others c3nfx')
netcdf.putAtt(ncid, others_c3nfx_spk_irrigated_id, 'units', '1')
netcdf.endDef(ncid);
netcdf.reDef(ncid);
others_c3nfx_fld_irrigated_id = netcdf.defVar(ncid, 'others_c3nfx_fld_irrigated', 'NC_FLOAT', [lonDimid, latDimid, timeDimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid, others_c3nfx_fld_irrigated_id, others_c3nfx_fld_irrigated);
netcdf.reDef(ncid);
netcdf.putAtt(ncid, others_c3nfx_fld_irrigated_id, 'longname', 'flood irrigated others c3nfx')
netcdf.putAtt(ncid, others_c3nfx_fld_irrigated_id, 'units', '1')
netcdf.endDef(ncid);

netcdf.close(ncid)




