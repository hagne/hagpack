import os as _os
from atmPy.data_archives.arm import read_data as _atm_arm
from atmPy.general import timeseries as _timeseries

products = {'HT_tdmaapshyg_1um_hyg400_rh85v40':     {'info': 'f(RH) calculated from tdmaaps using hygroscopicity from tdmahyg'},
            'HT_tdmaapsscattcoeff_1um_550nm':       {'info': 'scattering (extinction) coefficient calculated from tdmaaps using refrective indeces from aosacsm '},
            'HT_tdmaapsmass_1um':                   {'info': 'aerosol mass concentration calculated from tdmaaps using densities from aosacsm'}}

def load_netCDF(folder, prod_name, time_window, site = 'sgp', verbose = False):
    all_files = _os.listdir(folder)

    all_ts = []
    for file in all_files:
        if _os.path.splitext(file)[-1] != '.cdf':
            txt = '\t %s is not a netCDF file ... skipping'%file
            if verbose:
                print(txt)
            continue

        if not _atm_arm._is_in_time_window(file, time_window, verbose):
            continue

        site_check = _atm_arm._is_site(file, site, verbose)
        if not site_check:
            continue

        # test for correct product
        if not file.split('.')[0][4:] == prod_name:
            continue

        ts = _timeseries.load_netCDF(folder + file)
        all_ts.append(ts)
    #     print('found one: ', folder + file)

    ts_concat = _timeseries.concat(all_ts)
    ts_concat.data.sort_index(inplace=True)
    return ts_concat

def check_availability(folder, prod_name, time_window=('1990-01-01', '2030-01-01'), verbose = False):
    out = _atm_arm.check_availability(folder, data_product=[prod_name], time_window = time_window,  custom_product_keys=[prod_name], verbose=verbose)
    return out