import os as _os
from atmPy.data_archives.arm import read_data as _atm_arm
from atmPy.general import timeseries as _timeseries
import pandas as _pd
import numpy as _np
import warnings as _warnings
_warnings.catch_warnings()
_warnings.simplefilter("ignore")

products = {'HT_tdmaapshyg_1um_hyg400_rh85v40':     {'info': 'f(RH) calculated from tdmaaps using hygroscopicity from tdmahyg'},
            'HT_tdmaapsscattcoeff_1um_550nm':       {'info': 'OLD scattering (extinction) coefficient calculated from tdmaaps using refrective indeces from aosacsm '},
            'HT_tdmaapsmass_1um':                   {'info': 'aerosol mass concentration calculated from tdmaaps using densities from aosacsm'},
            'HT_tdmaapsbackscatt_1um_550nm':        {'info': 'hemispheric backscattering calculated from tdmaaps using index of refraction from aosacsm'}
            }

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

        fname = folder + file
        #print(fname)
        ts = _timeseries.load_netCDF(fname)
        all_ts.append(ts)
    #     print('found one: ', folder + file)

    ts_concat = _timeseries.concat(all_ts)
    ts_concat.data.sort_index(inplace=True)
    return ts_concat

def _check_availability(folder, prod_name, time_window=('1990-01-01', '2030-01-01'), verbose = False):
    out = _atm_arm.check_availability(folder, data_product=[prod_name], time_window = time_window,  custom_product_keys=[prod_name], verbose=verbose)
    return out

##########################################################################################
#########################
##### making the products

##############
### some tools
def _splitup_filename(filename):

    fname_split = _os.path.split(filename)[-1].split('.')
    site = fname_split[0][:3]
    date = fname_split[-3]
    out = {'site': site,
           'date': date}
    return out

def _get_other_filenames(filename, others, all_files, verbose = False):
    out_dict = {}
    splitname = _splitup_filename(filename)
    site = splitname['site']
    date = splitname['date']

    all_that_day = all_files[_np.char.find(all_files, date) > -1]

    for prod in others:
        fname = all_that_day[_np.char.find(all_that_day, prod) > -1]
        if fname.shape[0] != 1:
            if verbose:
                print('no corresponding aosacsm found ... continue on %s'%date)
            return False

        else:
            out = {'fname': fname[0]}
            out_dict[prod] = out
    out_dict['date']= date
    out_dict['site'] = site
    return out_dict

################
### The products

##################
## f of RH related
class TdmaapsTdmahygAosacsm2fofrh(object):
    def __init__(self, data_quality='patchy',
                 diameter_cutoff='1um',
                 hygroscopicity_diameter=400,
                 RH_dry=40,
                 RH_wet=85,
                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                 folder='/Users/htelg/data/ARM/SGP/',
                 test = False):

        self.data_quality = data_quality
        self.diameter_cutoff = diameter_cutoff
        self.hygroscopicity_diameter = hygroscopicity_diameter
        self.RH_dry = RH_dry
        self.RH_wet = RH_wet
        self.folder_out = folder_out
        self.folder = folder
        self.test = test
        self.test_file = 'sgptdmaapssizeC1.c1.20120201.002958.cdf'

    def _calculate_one(self, tdmaapssize, tdmahyg, aosacsm, diameter_cutoff,
                       hygroscopicity_diameter, RH_dry, RH_wet):

        tdmaapssize.size_distribution.index_of_refraction = aosacsm.refractive_index
        tdmaapssize.size_distribution.physical_property_density = aosacsm.density

        if diameter_cutoff == '1um':
            dcoff = 1000
        elif diameter_cutoff == '10um':
            dcoff = 10000

        dist = tdmaapssize.size_distribution.zoom_diameter(end=dcoff)

        which_diameters = [hygroscopicity_diameter]

        df = tdmahyg.mean_growth_factor.data.loc[:, :, 'mean']

        #     opt = dist.calculate_optical_properties(550)

        fofrh = _pd.DataFrame(index=dist.data.index)
        for wd in which_diameters:
            #         gf = df.loc[wd, : ]

            #         dist_g = dist.apply_growth(gf)
            kappa = tdmahyg.kappa_values.align_to(dist)
            kappa = kappa.data[wd]

            dist_d = dist.apply_hygro_growth(kappa, RH_dry)
            opt_d = dist_d.calculate_optical_properties(550)

            dist_w = dist.apply_hygro_growth(kappa, RH_wet)
            opt_w = dist_w.calculate_optical_properties(550)

            f_RH_tdma = opt_w.extinction_coeff_sum_along_d.data / opt_d.extinction_coeff_sum_along_d.data
            fofrh[wd] = f_RH_tdma

        fofrh.columns.name = 'diameter regime'
        fofrh = _timeseries.TimeSeries(fofrh)
        fofrh._data_period = dist._data_period
        fofrh._y_label = '$f_{RH = %s/%s}$' % (RH_wet, RH_dry)
        return fofrh

    def calculate_new(self):
        self._calculate_all(False)

    def calculate_all(self):
        self._calculate_all(True)

    def _calculate_all(self, overwrite, time_window=False, verbose=False):
        fofrh_list = []
        all_files = _os.listdir(self.folder)
        all_files = _np.array(all_files)

        all_files_tdmaapssize = all_files[_np.char.find(all_files, 'tdmaapssize') > -1]

        test_done = False

        for e, fname_tdmaapssize in enumerate(all_files_tdmaapssize):
            if self.test:
                if fname_tdmaapssize == self.test_file:
                    test_done = True
                else:
                    if test_done:
                        break
                    else:
                        continue

            if time_window:
                if not _atm_arm._is_in_time_window(fname_tdmaapssize, verbose):
                    continue

            fname_others = _get_other_filenames(fname_tdmaapssize, ['aosacsm', 'tdmahyg'], all_files)
            if not fname_others:
                continue

            name_addon = '%s_hyg%s_rh%sv%s' % (self.diameter_cutoff, self.hygroscopicity_diameter, self.RH_wet, self.RH_dry)
            name_addon  # = '1um_hyg400_rh85v40'

            my_prod_name = self.folder_out + '/' + fname_others['site'] + '_HT_tdmaapshyg_' + name_addon + '.' + fname_others['date'] + '.000000.cdf'

            if not overwrite:
                if _os.path.isfile(my_prod_name):
                    if verbose:
                        print('product %s already exists' % my_prod_name)
                    continue
            verbose = False
            tdmaapssize = _atm_arm.read_cdf(self.folder + fname_tdmaapssize, data_quality=self.data_quality, verbose=verbose)
            aosacsm = _atm_arm.read_cdf(self.folder + fname_others['aosacsm']['fname'], data_quality=self.data_quality, verbose=verbose)
            tdmahyg = _atm_arm.read_cdf(self.folder + fname_others['tdmahyg']['fname'], data_quality=self.data_quality, verbose=verbose)

            fofrh = self._calculate_one(tdmaapssize, tdmahyg, aosacsm, diameter_cutoff=self.diameter_cutoff,
                                        hygroscopicity_diameter=self.hygroscopicity_diameter, RH_dry=self.RH_dry, RH_wet=self.RH_wet)

            fofrh_list.append(fofrh)
            fofrh.save_netCDF(my_prod_name)
            # if self.test:
            #     if len(fofrh_list) == 1:
            #         break

        print(my_prod_name)
        fofrh_cat = _timeseries.concat(fofrh_list)
        self.result = fofrh_cat.close_gaps(verbose=False)
        return  # all_files_tdmaapssize

tdmaaps_tdmahyg_aosacsm2fofrh_1um_hyg400_rh85v40 = TdmaapsTdmahygAosacsm2fofrh(  data_quality='patchy',
                                                                 diameter_cutoff='1um',
                                                                 hygroscopicity_diameter=400,
                                                                 RH_dry=40,
                                                                 RH_wet=85,
                                                                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                                                 folder='/Users/htelg/data/ARM/SGP/',
                                                                 test = False)


class Tdmaaps2scatteringcoeff(object):
    def __init__(self, data_quality='patchy',
                 diameter_cutoff='1um',
                 wavelength=550,  # in nm
                 refractive_index = 'aosacsm',
                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                 folder='/Users/htelg/data/ARM/SGP/',
                 test = False):
        """

        Parameters
        ----------
        data_quality
        diameter_cutoff
        wavelength
        refractive_index: string or float
            Define which data product to get the RI from ['aosacsm']
            or set a fixed value.
        folder_out
        folder
        test
        """
        self.wavelength = wavelength
        self.refractive_index = refractive_index
        self.data_quality = data_quality
        self.diameter_cutoff = diameter_cutoff
        self.folder_out = folder_out
        self.folder = folder
        self.test = test
        self.test_file = 'sgptdmaapssizeC1.c1.20120201.002958.cdf'

    def _calculate_one(self, tdmaapssize, refractive_index):
        if self.diameter_cutoff == '1um':
            dcoff = 1000
        elif self.diameter_cutoff == '10um':
            dcoff = 10000

        dist = tdmaapssize.size_distribution.zoom_diameter(end=dcoff)
        dist.index_of_refraction = refractive_index
        opt = dist.calculate_optical_properties(self.wavelength)
        extcoeff = opt.extinction_coeff
        extcoeff.data.values[extcoeff.data.values == 0] = _np.nan
        extcoeff.data *= 1e6  # conversion from inverse meter to mega meter
        return extcoeff

    def calculate_new(self):
        self._calculate_all(False)

    def calculate_all(self):
        self._calculate_all(True)

    def _calculate_all(self, overwrite, time_window=False, verbose=False):

        extcoeff_list = []
        all_files = _os.listdir(self.folder)
        all_files = _np.array(all_files)

        all_files_tdmaapssize = all_files[_np.char.find(all_files, 'tdmaapssize') > -1]
        test_done = False
        for e, fname_tdmaapssize in enumerate(all_files_tdmaapssize):
            if self.test:
                if fname_tdmaapssize == self.test_file:
                    test_done = True
                else:
                    if test_done:
                        break
                    else:
                        continue
            if time_window:
                if not _atm_arm._is_in_time_window(fname_tdmaapssize, verbose):
                    continue

            splitname = _splitup_filename(fname_tdmaapssize)
            site = splitname['site']
            date = splitname['date']
            name_addon = ('RI%s_%s_%snm' % (self.refractive_index, self.diameter_cutoff, self.wavelength)).replace('.','o')
            my_prod_name = self.folder_out + site + 'tdmaaps2scatteringcoeff_' + name_addon + '.' + date + '.000000.cdf'

            if not overwrite:
                if _os.path.isfile(my_prod_name):
                    if verbose:
                        print('product %s already exists' % my_prod_name)
                    continue

            if self.refractive_index == 'aosacsm':
                fname_others = _get_other_filenames(fname_tdmaapssize, ['aosacsm'], all_files)
                if not fname_others:
                    continue
                else:
                    aosacsm = _atm_arm.read_cdf(self.folder + fname_others['aosacsm']['fname'], data_quality=self.data_quality, verbose=verbose)
                    refractive_index = aosacsm.refractive_index
            elif type(self.refractive_index).__name__ == 'float':
                refractive_index = self.refractive_index

            tdmaapssize = _atm_arm.read_cdf(self.folder + fname_tdmaapssize, data_quality=self.data_quality, verbose=verbose)

            extcoeff = self._calculate_one(tdmaapssize, refractive_index)

            extcoeff_list.append(extcoeff)
            if not self.test:
                extcoeff.save_netCDF(my_prod_name)
        # if len(extcoeff_list) == 2:
        #                 break


        print(my_prod_name)
        if len(extcoeff_list) > 1:
            extcoeff_cat = _timeseries.concat(extcoeff_list)
            self.result = extcoeff_cat.close_gaps(verbose=False)
        else:
            self.result = extcoeff_list[0]

def tdmaaps2scatteringcoeff_RIaosacsm_1um_550nm():
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='1um',
                                  wavelength=550,
                                  refractive_index='aosacsm',
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Users/htelg/data/ARM/SGP/',
                                  test = False)
    return out

def tdmaaps2scatteringcoeff_RI1o5_1um_550nm():
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='1um',
                                  wavelength=550,
                                  refractive_index=1.5,
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Users/htelg/data/ARM/SGP/',
                                  test = False)
    return out