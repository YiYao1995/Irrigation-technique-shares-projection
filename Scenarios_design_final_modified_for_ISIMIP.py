import pandas as pd
import csv
import numpy as np
import scipy.io as scio

# scenarios = ['126', '370', '585']   # currently we only have data for these three scenarios
# scenarios = ['126', '370']
scenarios = ['126', '370', '585']
def read_csv_data_SSP(filename, encode='ANSI'):  # get the socioeconomic data and ISIMIP data from csv
    pd_reader = pd.read_csv(filename, encoding=encode, header=None)
    data = np.array(pd_reader)
    return data[:, 1:]

for scenario in scenarios[:]:
    print(scenario)
    # This is the socio data after PCA and Norm
    socio_data = np.loadtxt(open("C:\Research2\Socioeconomic\Socio_data_SSP" + scenario + ".csv", "rb"), delimiter=",", skiprows=0)     # get the file of socioeconomic data

    Map_code_proj = {}  # get the dict of country_name - code
    with open('C:\Research2\surf_data\country_index.csv',
              encoding='ANSI') as inp:  # read csv as a dictionary
        reader = csv.reader(inp)
        Map_code_proj = {rows[0]: rows[1] for rows in reader}

    Socio_code_file = {}    # get the dict of country_code - socio_code
    with open('C:\Research2\Socioeconomic\Socioeconomic_codes.csv',
              encoding='ANSI') as inp:  # read csv as a dictionary
        reader = csv.reader(inp)
        Socio_code_file = {rows[0]: rows[1] for rows in reader}

    # AREA and PCT_CROP data from the CLM surface dataset
    # hist_AREA_dict = scio.loadmat('C:\Research2\PCT_CFT_AREA\SSP' + scenario + '\AREA.mat') # no problem
    # hist_AREA = hist_AREA_dict['AREA']  # grid area data
    # hist_PCT_CROP_dict = scio.loadmat('C:\Research2\PCT_CFT_AREA\SSP' + scenario + '\PCT_CROP.mat') # no problem
    # hist_PCT_CROP = hist_PCT_CROP_dict['PCT_CROP']  # percentage of crop land

    # irrigation area and frac, based on Jonas'data and CLM data. THE YEAR 2010
    # basic_drip_area_dict = scio.loadmat('C:\Research2\PCT_CFT_AREA\SSP' + scenario + '\\act_drip_area.mat')
    # basic_drip_area = basic_drip_area_dict['act_drip_area']
    #
    str_map_data = 'D:\\ISIMIP\For_generation\landuse-15crops_image_ssp' + scenario + '_annual_2015_2100'
    basic_drip_frac_dict = scio.loadmat(str_map_data + '\\act_irr_drip.mat') # no problem
    basic_drip_frac = basic_drip_frac_dict['act_irr_drip']
    #
    # basic_spri_area_dict = scio.loadmat('C:\Research2\PCT_CFT_AREA\SSP' + scenario + '\\act_spri_area.mat')
    # basic_spri_area = basic_spri_area_dict['act_spri_area']

    basic_spri_frac_dict = scio.loadmat(str_map_data + '\\act_irr_spri.mat')
    basic_spri_frac = basic_spri_frac_dict['act_irr_spri']

    # basic_floo_area_dict = scio.loadmat('C:\Research2\ISIMIP\ssp5_noadapt\\act_floo_area.mat')
    # basic_floo_area = basic_floo_area_dict['act_floo_area']

    basic_floo_frac_dict = scio.loadmat(str_map_data + '\\act_irr_floo.mat')
    basic_floo_frac = basic_floo_frac_dict['act_irr_floo']


    # all_CROP_AREA_dict = scio.loadmat('C:\Research2\PCT_CFT_AREA\SSP' + scenario + '\\CROP_AREA.mat') # no problem
    # all_CROP_AREA = all_CROP_AREA_dict['CROP_AREA']

    # get the optimal irrigation techniques share. based on CFT types
    opt_drip_dict = scio.loadmat(str_map_data + '\\opt_drip.mat')
    opt_drip = opt_drip_dict['opt_drip']

    opt_floo_dict = scio.loadmat(str_map_data + '\\opt_floo.mat')
    opt_floo = opt_floo_dict['opt_floo']

    opt_spri_dict = scio.loadmat(str_map_data + '\\opt_spri.mat')
    opt_spri = opt_spri_dict['opt_spri']

    pct_irr_dict = scio.loadmat(str_map_data + '\\pct_irr.mat')
    pct_irr = pct_irr_dict['pct_irr']

    # basic fraction data
    # basic_drip_dict = scio.loadmat('C:\Research2\\basic_fraction\drip.mat')
    # basic_sprinkler_dict = scio.loadmat('C:\Research2\\basic_fraction\sprinkler.mat')
    # basic_flood_dict = scio.loadmat('C:\Research2\\basic_fraction\\flood.mat')
    # basic_drip = basic_drip_dict['drip_fraction']
    # basic_sprinkler = basic_sprinkler_dict['sprinkler_fraction']
    # basic_flood = basic_flood_dict['flood_fraction']

    # strings for P data
    # str_P_begin = 'C:\Research2\ISIMIP3a\\P_data\ssp' + scenario + '\\'
    #
    # str_year_list = ['1996_2015', '2001_2020', '2006_2025', '2011_2030', '2016_2035', '2021_2040', '2026_2045', '2031_2050',
    #                  '2036_2055', '2041_2060', '2046_2065', '2051_2070', '2056_2075', '2061_2080', '2066_2085', '2071_2090',
    #                  '2076_2095', '2081_2100']
    # str_P_end = '_after_normalization.csv'
    #
    # # OUTPUT data
    grid_frac = np.zeros([720, 360, 19, 3])
    grid_area = np.zeros([288, 192, 19, 3])
    # country_frac = np.zeros([257, 18, 3])
    # country_area = np.zeros([257, 19, 3])

    # for i_country in range(1, 257):
    for i_country in range(1, 257):
        print(str(i_country))
        try:
            # get the pixels of the country
            str_file = 'C:\Research2\ISIMIP\countries\\' + str(i_country) + '.csv'
            row_col = pd.read_csv(str_file, header=None, delimiter=',', index_col=None)
            print(type(row_col))
            all_str = row_col.to_numpy()
            print(all_str)
            print('all_str[0,0] = ', all_str[0,0])
            print(i_country)
            try:
                code = Map_code_proj[str(i_country)]
                print(code)
                socio_code = int(Socio_code_file[code])

                # all_area = np.zeros([19])
                # drip_area = np.zeros([19])
                # spri_area = np.zeros([19])
                # floo_area = np.zeros([19])


                for line in range(len(all_str)):

                    x = int(all_str[line, 0]) - 1
                    y = int(all_str[line, 1]) - 1
                    print('x = ' + str(x))
                    print('y = ' + str(y))
                    if pct_irr[x, y, 0] == 0:  # if there is no crop area, let it go
                        continue
                    else:
                        flood = basic_floo_frac[x, y]
                        sprinkler = basic_spri_frac[x, y]
                        drip = basic_drip_frac[x, y]

                        print('Basic flood frac = ' + str(flood))
                        print('Basic sprinkler frac = ' + str(sprinkler))
                        print('Basic drip frac = ' + str(drip))

                        if np.isnan(flood):  # needs to be checked
                            try:
                                if not (np.isnan(basic_floo_frac[x - 1, y])):
                                    flood = basic_floo_frac[x - 1, y]
                                    sprinkler = basic_spri_frac[x - 1, y]
                                    drip = basic_drip_frac[x - 1, y]
                                elif not (np.isnan(basic_floo_frac[x + 1, y])):
                                    flood = basic_floo_frac[x + 1, y]
                                    sprinkler = basic_spri_frac[x + 1, y]
                                    drip = basic_drip_frac[x + 1, y]
                                elif not (np.isnan(basic_floo_frac[x, y - 1])):
                                    flood = basic_floo_frac[x, y - 1]
                                    sprinkler = basic_spri_frac[x, y - 1]
                                    drip = basic_drip_frac[x, y - 1]
                                elif not (np.isnan(basic_floo_frac[x, y + 1])):
                                    flood = basic_floo_frac[x, y + 1]
                                    sprinkler = basic_spri_frac[x, y + 1]
                                    drip = basic_drip_frac[x, y + 1]
                                else:  # needed to be checked
                                    flood = 1
                                    sprinkler = 0
                                    drip = 0
                            except (IndexError):
                                flood = 1
                                sprinkler = 0
                                drip = 0

                        grid_frac[x, y, 0, 0] = flood
                        grid_frac[x, y, 0, 1] = sprinkler
                        grid_frac[x, y, 0, 2] = drip

                        # grid_area[x, y, 0, 0] = all_CROP_AREA[x, y, 9] * flood
                        # grid_area[x, y, 0, 1] = all_CROP_AREA[x, y, 9] * sprinkler
                        # grid_area[x, y, 0, 2] = all_CROP_AREA[x, y, 9] * drip

                        # all_area[0] = all_area[0] + all_CROP_AREA[x, y, 9]
                        # drip_area[0] = drip_area[0] + all_CROP_AREA[x, y, 9] * drip
                        # spri_area[0] = spri_area[0] + all_CROP_AREA[x, y, 9] * sprinkler
                        # floo_area[0] = floo_area[0] + all_CROP_AREA[x, y, 9] * flood

                # floo_country = floo_area[0] / all_area[0]
                # spri_country = spri_area[0] / all_area[0]
                # drip_country = drip_area[0] / all_area[0]

                # country_area[i_country, 0, 0] = floo_area[0]
                # country_area[i_country, 0, 1] = spri_area[0]
                # country_area[i_country, 0, 2] = drip_area[0]

                # if all_area[0] != 0:
                #     country_frac[i_country, 0, 0] = floo_area[0] / all_area[0]
                #     country_frac[i_country, 0, 1] = spri_area[0] / all_area[0]
                #     country_frac[i_country, 0, 2] = drip_area[0] / all_area[0]

                for time in range(18):
                    real_year = (time) * 5

                    for line in range(len(all_str)):
                        x = int(all_str[line, 0]) - 1
                        y = int(all_str[line, 1]) - 1
                        print('x = ' + str(x))
                        print('y = ' + str(y))

                        if pct_irr[x, y, real_year] == 0: # if there is no crop area, let it go
                            continue
                        else:
                            #flood = grid_frac[x, y, time-1, 0]
                            #sprinkler = grid_frac[x, y, time-1, 1]
                            #drip = grid_frac[x, y, time-1, 2] time-1 -> time
                            flood = grid_frac[x, y, time, 0]
                            sprinkler = grid_frac[x, y, time, 1]
                            drip = grid_frac[x, y, time, 2]
                            if np.isnan(flood):  # needs to be checked
                                try:
                                    if not (np.isnan(grid_frac[x - 1, y, time, 0])):
                                        flood = grid_frac[x - 1, y, time, 0]
                                        sprinkler = grid_frac[x - 1, y, time, 1]
                                        drip = grid_frac[x - 1, y, time, 2]
                                    elif not (np.isnan(grid_frac[x + 1, y, time, 0])):
                                        flood = grid_frac[x + 1, y, time, 0]
                                        sprinkler = grid_frac[x + 1, y, time, 1]
                                        drip = grid_frac[x + 1, y, time, 2]
                                    elif not (np.isnan(grid_frac[x, y - 1, time, 0])):
                                        flood = grid_frac[x, y - 1, time, 0]
                                        sprinkler = grid_frac[x, y - 1, time, 1]
                                        drip = grid_frac[x, y - 1, time, 2]
                                    elif not (np.isnan(grid_frac[x, y + 1, time, 0])):
                                        flood = grid_frac[x, y + 1, time, 0]
                                        sprinkler = grid_frac[x, y + 1, time, 1]
                                        drip = grid_frac[x, y + 1, time, 2]
                                    else:       # needed to be checked
                                        flood = 1
                                        sprinkler = 0
                                        drip = 0
                                except (IndexError):
                                    flood = 1
                                    sprinkler = 0
                                    drip = 0

                        # str_year = str_year_list[time]
                        # str_P_file = str_P_begin + str_year + str_P_end

                        optimal_flood = opt_floo[x, y, real_year]
                        optimal_sprinkler = opt_spri[x, y, real_year]
                        optimal_drip = opt_drip[x, y, real_year]

                        print('Optimal flood frac = ' + str(optimal_flood))
                        print('Optimal sprinkler frac = ' + str(optimal_sprinkler))
                        print('Optimal drip frac = ' + str(optimal_drip))

                        # pd_reader = pd.read_csv(str_P_file, header=None)
                        # var_P = np.array(pd_reader)
                        #
                        # P = var_P[x, y]
                        # hydro_climate = 1 - P

                        socio_eco = socio_data[socio_code, time + 2]  # socio_data needs to be created and loaded

                        speed = 1
                        if socio_eco < 0.5:
                            speed = speed - 0.8
                        elif socio_eco < 0.75:
                            speed = speed - 0.4
                        elif socio_eco < 1:
                            speed = speed
                        elif socio_eco < 1.25:
                            speed = speed + 0.4
                        else:
                            speed = speed + 0.8

                        # if hydro_climate < 0.25:
                        #     speed = speed - 0.2
                        # elif hydro_climate < 0.5:
                        #     speed = speed - 0.1
                        # elif hydro_climate < 0.75:
                        #     speed = speed
                        # elif hydro_climate < 1:
                        #     speed = speed + 0.1
                        # else:
                        #     speed = speed + 0.2

                        speed = speed / 100

                        spri_add_ratio = (optimal_sprinkler) / (optimal_sprinkler + optimal_drip)
                        drip_add_ratio = (optimal_drip) / (optimal_sprinkler + optimal_drip)
                        if np.isnan(spri_add_ratio):
                            spri_add_ratio = 0
                            drip_add_ratio = 0
                        #   flood = flood - 5 * speed
                        sprinkler = sprinkler + 5 * speed * spri_add_ratio
                        drip = drip + 5 * speed * drip_add_ratio

                        print('Updated flood frac = ' + str(flood))
                        print('Updated sprinkler frac = ' + str(sprinkler))
                        print('Updated drip frac = ' + str(drip))

                        if drip > optimal_drip:
                            drip = optimal_drip
                        if sprinkler + drip > optimal_sprinkler + optimal_drip:
                            sprinkler = optimal_sprinkler + optimal_drip - drip


                        flood = 1 - drip - sprinkler
                        print('Adjusted flood frac = ' + str(flood))
                        print('Adjusted sprinkler frac = ' + str(sprinkler))
                        print('Adjusted drip frac = ' + str(drip))

                        grid_frac[x, y, time+1, 0] = flood
                        grid_frac[x, y, time+1, 1] = sprinkler
                        grid_frac[x, y, time+1, 2] = drip

                    #     grid_area[x, y, time+1, 0] = all_CROP_AREA[x, y, real_year] * flood
                    #     grid_area[x, y, time+1, 1] = all_CROP_AREA[x, y, real_year] * sprinkler
                    #     grid_area[x, y, time+1, 2] = all_CROP_AREA[x, y, real_year] * drip
                    #
                    #     all_area[time+1] = all_area[time+1] + all_CROP_AREA[x, y, real_year]
                    #     drip_area[time+1] = drip_area[time+1] + all_CROP_AREA[x, y, real_year] * drip
                    #     spri_area[time+1] = spri_area[time+1] + all_CROP_AREA[x, y, real_year] * sprinkler
                    #     floo_area[time+1] = floo_area[time+1] + all_CROP_AREA[x, y, real_year] * flood
                    #
                    # floo_country = floo_area[time+1] / all_area[time+1]
                    # spri_country = spri_area[time+1] / all_area[time+1]
                    # drip_country = drip_area[time+1] / all_area[time+1]
                    #
                    # country_area[i_country, time+1, 0] = floo_area[time+1]
                    # country_area[i_country, time+1, 1] = spri_area[time+1]
                    # country_area[i_country, time+1, 2] = drip_area[time+1]
                    # if all_area[time+1] != 0:
                    #     country_frac[i_country, time+1, 0] = floo_area[time+1] / all_area[time+1]
                    #     country_frac[i_country, time+1, 1] = spri_area[time+1] / all_area[time+1]
                    #     country_frac[i_country, time+1, 2] = drip_area[time+1] / all_area[time+1]
            except (KeyError):
                error = 1
        except (FileNotFoundError):
            error = 1

    scio.savemat(str_map_data + '\grid_frac_ssp' + scenario + '.mat', {'grid_frac': grid_frac})
    # scio.savemat('C:\Research2\ISIMIP\ssp5_noadapt\grid_area_ssp' + scenario + '.mat', {'grid_area': grid_area})
    # scio.savemat('C:\Research2\Output_data\country_frac_ssp' + scenario + '.mat', {'country_frac': country_frac})
    # scio.savemat('C:\Research2\Output_data\country_area_ssp' + scenario + '.mat', {'country_area': country_area})