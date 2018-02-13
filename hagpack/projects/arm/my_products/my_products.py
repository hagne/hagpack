import os as _os
from atmPy.data_archives.arm import _read_data as _atm_arm
from atmPy.general import timeseries as _timeseries
from atmPy.aerosols.physics import hygroscopicity as _hygroscopicity
import pandas as _pd
import numpy as _np
import warnings as _warnings
from atmPy import read_file
_warnings.catch_warnings()
_warnings.simplefilter("ignore")

products = {#### hygroscopicity
            'tdmaapstdmahyg2fofrh_1um_hyg400_rh85v0_ior1o5_patchy':
                {'info': ''},
            'tdmaapstdmahyg2fofrh_10um_hyg400_rh85v0_ior1o5_patchy':
                {'info': ''},
            'tdmaapstdmahyg2fofrh_10um_hyg400_rh85v0_ioraosacsm_patchy':
                {'info': ''},
            'noaaaos2hygroscopicity_10um_550nm_patchy':
                {'info': 'calculates hygroscopicity (kapp, gamma, f_rh(85), std of everything) from the dry and wet neph data within the noaaaos product'},
            'noaaaos2hygroscopicity_1um_550nm_patchy':
                {'info': ''},
            #### backscattering ratio
            'tdmaaps2backscatteringratio_RI1o5_1um_550nm':
                {'info': ''},
            'tdmaaps2backscatteringratio_RI1o5_10um_550nm':
                {'info': ''},
            'tdmaaps2backscatteringratio_bc_abs_noaaaos_RI1o5_10um_550nm':
                {'info': ''},
            'tdmaaps2backscatteringratio_RIaosacsm_10um_550nm':
                {'info': ''},
            'tdmaaps2backscatteringratio_bc_abs_noaaaos_RIaosacsm_10um_550nm':
                {'info': ''},
            ####scattering
            'tdmaaps2scatteringcoeff_bc_abs_noaaaos_RIaosacsm_10um_550nm':
                {'info': 'as tdmaaps2scatteringcoeff_bc_abs_noaaaos_RI1o5_10um_550nm but with the index of refraction taken from the acsm'},
            'tdmaaps2scatteringcoeff_bc_abs_noaaaos_RI1o5_10um_550nm':
                {'info': 'as tdmaaps2scatteringcoeff_bc_ratio_0o050_RI1o5_10um_550nm just withthe absorption perscribed by the noaaaos psap measurements'},
            'tdmaaps2scatteringcoeff_bc_ratio_0o050_RI1o5_10um_550nm':
                {'info': 'as tdmaaps2scatteringcoeff_RIaosacsm_1um_550nm with 5% of all particles to be assumed black carbon'},
            'tdmaaps2scatteringcoeff_RIaosacsm_1um_550nm':
                {'info': 'scattering (extinction) coefficient calculated from tdmaaps using refrective indeces from aosacsm. Particle diameters considered up to 1um. '},
            'tdmaaps2scatteringcoeff_RIaosacsm_10um_550nm':
                {'info': ''},
            'tdmaaps2scatteringcoeff_RIaosacsm_1um_550nm_good':
                {'info': 'scattering (extinction) coefficient calculated from tdmaaps using refrective indeces from aosacsm '},
            'tdmaaps2scatteringcoeff_RI1o5_1um_550nm':
                {'info': 'Refractive index is fixed to 1.5. Particle diameters considered up to 1um.'},
            'tdmaaps2scatteringcoeff_RI1o5_10um_550nm':
                {'info': 'scattering (extinction) coefficient calculated from tdmaaps using refrective indeces of 1.5 . Particle diameters considered up to 10um.'},
            #### Kappa
            'aipfitrh2kappa_RH_85_40_tdmaapssize_RI1o5_1um_550nm_patchy':
                {'info': 'kappa from aipfitrh'},
            'tdmahyg2kappa_avg_d200_patchy':
                {'info': 'kappa from tdmahyg'}
            # 'HT_tdmaapshyg_1um_hyg400_rh85v40':     {'info': 'f(RH) calculated from tdmaaps using hygroscopicity from tdmahyg'},
            # 'HT_tdmaapsscattcoeff_1um_550nm':       {'info': 'OLD scattering (extinction) coefficient calculated from tdmaaps using refrective indeces from aosacsm '},
            # 'HT_tdmaapsmass_1um':                   {'info': 'aerosol mass concentration calculated from tdmaaps using densities from aosacsm'},
            # 'HT_tdmaapsbackscatt_1um_550nm':        {'info': 'hemispheric backscattering calculated from tdmaaps using index of refraction from aosacsm'}
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
        if not file.split('.')[0][3:] == prod_name:
            if not file.split('.')[0] == prod_name:
                continue

        fname = folder + file
        #print(fname)
        # ts = _timeseries.load_netCDF(fname)
        ts = read_file.netCDF(fname)
        all_ts.append(ts)
    #     print('found one: ', folder + file)
    if len(all_ts) == 0:
        raise ValueError('no file meets criteria')
    ts_concat = _timeseries.concat(all_ts)
    ts_concat.data.sort_index(inplace=True)
    ts_concat = ts_concat.close_gaps()
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
                print('no corresponding file found ... continue on %s'%date)
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

class TdmaapsTdmahyg2fofrh(object):
    """using an effective kappa from the growthdistributions in tdmahyg"""
    def __init__(self, data_quality='patchy',
                 diameter_cutoff='1um',
                 hygroscopicity_diameter=400,
                 RH_dry=0,
                 RH_wet=85,
                 refractive_index ='aosacsm',
                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                 folder='/Users/htelg/data/ARM/SGP/',
                 keep_data = False,
                 test = False):

        self.data_quality = data_quality
        self.diameter_cutoff = diameter_cutoff
        self.hygroscopicity_diameter = hygroscopicity_diameter
        self.RH_dry = RH_dry
        self.RH_wet = RH_wet
        self.refractive_index = refractive_index
        self.black_carbon = None
        self.folder_out = folder_out
        self.folder = folder
        self.test = test
        self.test_file = 'sgptdmaapssizeC1.c1.20120201.002958.cdf'
        self.keep_data = keep_data
        if self.test:
            self.keep_data = True
            self.intres = {}

    def _calculate_one(self, tdmaapssize, tdmahyg, refractive_index, diameter_cutoff,
                       hygroscopicity_diameter, RH_dry, RH_wet):

        tdmaapssize.size_distribution.index_of_refraction = refractive_index
        # tdmaapssize.size_distribution.physical_property_density = aosacsm.density

        if diameter_cutoff == '1um':
            dcoff = 1000
        elif diameter_cutoff == '10um':
            dcoff = 10000

        dist = tdmaapssize.size_distribution.zoom_diameter(end=dcoff)
        which_diameters = [hygroscopicity_diameter]
        # df = tdmahyg.mean_growth_factor.data.loc[:, :, 'mean']
        #     opt = dist.calculate_optical_properties(550)

        fofrh = _pd.DataFrame(index=dist.data.index)
        for wd in which_diameters:
            #         gf = df.loc[wd, : ]

            #         dist_g = dist.apply_growth(gf)
            kappa = tdmahyg.kappa_values.align_to(dist)
            kappa = kappa._del_all_columns_but(wd)

            # self.intres['dist'] = dist
            # self.intres['kappa'] = kappa
            # self.intres['RH_wet'] = RH_wet
            if RH_dry == 0:
                dist_d = dist.copy()
            else:
                dist_d = dist.apply_hygro_growth(kappa, RH_dry)
            opt_d = dist_d.calculate_optical_properties(550)

            dist_w = dist.apply_hygro_growth(kappa, RH_wet)
            opt_w = dist_w.calculate_optical_properties(550)
            # self.intres['opt_w'] = opt_w
            f_RH_tdma = opt_w.scattering_coeff.data / opt_d.scattering_coeff.data
            fofrh[wd] = f_RH_tdma

        fofrh.columns.name = 'diameter regime'
        fofrh = _timeseries.TimeSeries(fofrh)
        fofrh._data_period = dist._data_period
        fofrh._y_label = '$f_{RH = %s/%s}$' % (RH_wet, RH_dry)
        return fofrh
        # return dist

    def calculate_new(self):
        self._calculate_all(False)

    def calculate_all(self):
        self._calculate_all(True)

    def _calculate_all(self, overwrite, time_window=False, verbose=False):
        fofrh_list = []

        bc_name_signitur = ''
        if self.black_carbon:
            if 'bc_ratio' in self.black_carbon.keys():
                bc_ratio = self.black_carbon['bc_ratio']
                bc_name_signitur = ('_bc_ratio_%.3f'%bc_ratio).replace('.', 'o')
            elif 'perscribed_absorption' in self.black_carbon.keys():
                perscribed_absorption = self.black_carbon['perscribed_absorption']
                if type(perscribed_absorption).__name__ == 'str':
                    pat = perscribed_absorption
                elif type(perscribed_absorption).__name__ in ['int', 'float']:
                    pat = '%.1f'%perscribed_absorption
                bc_name_signitur = ('_bc_abs_%s' %pat).replace('.', 'o')

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
            ior = str(self.refractive_index).replace('.','o')
            name_addon = '{}_hyg{:d}_rh{:d}v{:d}_ior{}_{}'.format(self.diameter_cutoff,
                                                     self.hygroscopicity_diameter,
                                                     self.RH_wet,
                                                     self.RH_dry,
                                                     ior,
                                                     self.data_quality)
            name_addon  # = '1um_hyg400_rh85v40'

            splitname = _splitup_filename(fname_tdmaapssize)
            site = splitname['site']
            date = splitname['date']


            my_prod_name = self.folder_out + '/' + fname_others['site'] + '_HT_tdmaapshyg_' + name_addon + '.' + fname_others['date'] + '.000000.cdf'
            my_prod_name = self.folder_out + site + 'tdmaapstdmahyg2fofrh_' + bc_name_signitur + name_addon + '.' + date +  '.000000.cdf'
            if not overwrite:
                if _os.path.isfile(my_prod_name):
                    if verbose:
                        print('product %s already exists' % my_prod_name)
                    continue
            verbose = False

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

            fname_others = _get_other_filenames(fname_tdmaapssize, ['tdmahyg'], all_files)
            tdmahyg = _atm_arm.read_cdf(self.folder + fname_others['tdmahyg']['fname'], data_quality=self.data_quality, verbose=verbose)

            fofrh = self._calculate_one(tdmaapssize, tdmahyg, refractive_index, diameter_cutoff=self.diameter_cutoff,
                                        hygroscopicity_diameter=self.hygroscopicity_diameter, RH_dry=self.RH_dry, RH_wet=self.RH_wet)

            # fofrh_list.append(fofrh)
            # if not self.test:
            #     fofrh.save_netCDF(my_prod_name)
            # if self.test:
            #     if len(fofrh_list) == 1:
            #         break

        # print(my_prod_name)
        # fofrh_cat = _timeseries.concat(fofrh_list)
        # self.result = fofrh_cat.close_gaps(verbose=False)
        # return  # all_files_tdmaapssize

            if self.keep_data:
                fofrh_list.append(fofrh)
            if not self.test:
                fofrh.save_netCDF(my_prod_name)


# if len(extcoeff_list) == 2:
#                 break


        print(my_prod_name)
        if self.keep_data:
            if len(fofrh_list) > 1:
                fofrh_cat = _timeseries.concat(fofrh_list)
                self.result = fofrh_cat.close_gaps(verbose=False)
            else:
                self.result = fofrh_list[0]

def tdmaapstdmahyg2fofrh_1um_hyg400_rh85v0_ior1o5_patchy(test = False):
    out = TdmaapsTdmahyg2fofrh(data_quality='patchy',
                                               diameter_cutoff='1um',
                                               hygroscopicity_diameter=400,
                                               RH_dry=0,
                                               RH_wet=85,
                                               refractive_index=1.5,
                                               folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                               folder='/Users/htelg/data/ARM/SGP/',
                                               test=test)
    return out

def tdmaapstdmahyg2fofrh_10um_hyg400_rh85v0_ior1o5_patchy(test = False):
    out = TdmaapsTdmahyg2fofrh(data_quality='patchy',
                                               diameter_cutoff='10um',
                                               hygroscopicity_diameter=400,
                                               RH_dry=0,
                                               RH_wet=85,
                                               refractive_index=1.5,
                                               folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                               folder='/Users/htelg/data/ARM/SGP/',
                                               test=test)
    return out

def tdmaapstdmahyg2fofrh_10um_hyg400_rh85v0_ioraosacsm_patchy(test = False):
    out = TdmaapsTdmahyg2fofrh(data_quality='patchy',
                                               diameter_cutoff='10um',
                                               hygroscopicity_diameter=400,
                                               RH_dry=0,
                                               RH_wet=85,
                                               refractive_index='aosacsm',
                                               folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                               folder='/Users/htelg/data/ARM/SGP/',
                                               test=test)
    return out








class Tdmaaps2scatteringcoeff(object):
    def __init__(self, data_quality='patchy',
                 diameter_cutoff='1um',
                 wavelength=550,  # in nm
                 refractive_index = 'aosacsm',
                 apply_growth = None,
                 black_carbon = False,
                 mode_selection = None,
                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                 folder='/Users/htelg/data/ARM/SGP/',
                 keep_data = False,
                 test = False,
                 verbose = False):
        """

        Parameters
        ----------
        data_quality
        diameter_cutoff
        wavelength
        refractive_index: string or float
            Define which data product to get the RI from ['aosacsm']
            or set a fixed value.
        absorption: dictionary
            'perscribed_absorption' in Mm^-1
            'abs_ratio'
            dict: 'noaaaos', will use the absorption measurements from the noaaaos product
            float: fraction to assume to be black carbon, e.g. 0.05 will assume 5% of the accumulation mode to be black carbon
        mode_selection: [None], {'moment':'volume', 'n_accu':'aosacsm', 'n_coarse': complex(1.63, 0.01), 'gf_coarse': 1.1}
            If mode selection is to be performed. Dict is expected with keys: 'moment', 'n_accu', 'n_coarse', where
            moment is the moment at which to perform the mode separation. Does not work for with black_carbon, or
            absorption. Following kwargs will be ignored: refractive_index, black_carbon
        folder_out
        keep_data: bool
            if the data is returned or discarded after saving it. Returning it can result in massive RAM use.
        folder
        test
        """
        self.verbose = verbose
        self.wavelength = wavelength
        self.refractive_index = refractive_index
        self.data_quality = data_quality
        self.diameter_cutoff = diameter_cutoff
        self.folder_out = folder_out
        self.folder = folder
        self.test = test
        if self.test:
            keep_data = True
        self.test_file = 'sgptdmaapssizeC1.c1.20120201.002958.cdf'
        self.keep_data = keep_data

        self.black_carbon = black_carbon
        if self.black_carbon:
            allowed_kwargs_bc = {'perscribed_absorption': [],
                                 'bc_ratio': []}

            for k in self.black_carbon:
                if k not in allowed_kwargs_bc.keys():
                    txt = '%s is an unknown key for kwarg black_carbon. Allowed keys are %s.'%(k, allowed_kwargs_bc.keys())
                    raise KeyError(txt)

        self.mode_sel = mode_selection
        if self.mode_sel:
            if self.black_carbon:
                raise KeyError('If performing mode_sel black_carbon is not allowed.')
            if self.refractive_index:
                raise KeyError('If performing mode_sel refractive_index defined in mode_sel, set refractive_index to None.')
            for k in self.mode_sel.keys():
                if k not in ['moment', 'n_accu', 'n_coarse', 'gf_coarse']:
        # if self.mode_sel.keys() != {'moment':'volume', 'n_accu':'aosacsm', 'n_coarse': 1.63}.keys():
                    raise KeyError("Mode selection dict not correct. Should look like this: {'moment':'volume', 'n_accu':'aosacsm', 'n_coarse': complex(1.63, 0.01)}")

            allowed_moments = ['volume']
            if self.mode_sel['moment'] not in allowed_moments:
                raise ValueError('Moment {} not implemented. Choose from {}.'.format(self.mode_sel['moment'], allowed_moments))

            allowed_n = ['aosacsm']
            if type(self.mode_sel['n_accu']) == str:
                if self.mode_sel['n_accu'] not in allowed_n:
                    raise ValueError('n_accu has to be in {}. It is {}'.format(allowed_n, self.mode_sel['n_accu']))
            if type(self.mode_sel['n_coarse']) == str:
                if self.mode_sel['n_coarse'] not in allowed_n:
                    raise ValueError('n_coarse has to be in {}. It is {}'.format(allowed_n, self.mode_sel['n_coarse']))

        self.apply_growth = apply_growth

        if test:
            self.intres = {}

    def _calculate_one(self, dist, refractive_index, abs_ratio, perscribed_absorption, growthfct):
        """

        Parameters
        ----------
        tdmaapssize: obvious
        refractive_index: obvious
        abs_ratio: ratio of the numberconcentration that is black carbon
        perscribed_absorption: this will result in an iteration in which the abs_ratio will be determined that is necessary so the absorbtion is == abs_goal

        Returns
        -------

        """
        if self.diameter_cutoff == '1um':
            dcoff = 1000
        elif self.diameter_cutoff == '10um':
            dcoff = 10000

        # if mode_sel:


        # dist = tdmaapssize.size_distribution.zoom_diameter(end=dcoff)
        dist = dist.zoom_diameter(end=dcoff)
        dist = dist.convert2numberconcentration()
        if growthfct:

            dist = _hygroscopicity.apply_growth2sizedist(dist, growthfct)
            if self.verbose:
                print('Applied growth factor of {}'.format(growthfct))

        if _np.any(perscribed_absorption):
            dist_bc = dist.copy()
            dist_bc.index_of_refraction = complex(2.26, 1.26)
            opt = dist_bc.calculate_optical_properties(self.wavelength)
            absorption_coeff_bc = opt.absorption_coeff
            absorption_coeff_bc.data.values[absorption_coeff_bc.data.values == 0] = _np.nan
            absorption_coeff_bc.data *= 1e6  # conversion from inverse meter to mega meter

            if type(perscribed_absorption).__name__ == 'TimeSeries':
                abs_ratio = (perscribed_absorption / absorption_coeff_bc).data['perscribed_absorption']
            else:
                abs_ratio = (perscribed_absorption / absorption_coeff_bc.data['abs_coeff_m^1'])


        if _np.any(abs_ratio):
            dist_bc = dist.copy()
            dist_bc.data = dist_bc.data.mul(abs_ratio, axis = 0)

            # dist_bc.index_of_refraction = complex(2.26,1.26)
            # opt = dist_bc.calculate_optical_properties(self.wavelength)
            # scattcoeff_bc = opt.scattering_coeff

            dist_bc.optical_properties.parameters.refractive_index = complex(2.26,1.26)
            dist_bc.optical_properties.parameters.wavelength = self.wavelength
            scattcoeff_bc = dist_bc.optical_properties.scattering_coeff

            scattcoeff_bc.data.values[scattcoeff_bc.data.values == 0] = _np.nan
            scattcoeff_bc.data *= 1e6  # conversion from inverse meter to mega meter
        else:
            abs_ratio = 0

        dist.data = dist.data.mul((1 - abs_ratio), axis = 0)

        # dist.index_of_refraction = refractive_index
        # opt = dist.calculate_optical_properties(self.wavelength)
        dist.optical_properties.parameters.refractive_index = refractive_index
        dist.optical_properties.parameters.wavelength = self.wavelength

        scattcoeff = dist.optical_properties.scattering_coeff
        # scattcoeff = opt.scattering_coeff
        scattcoeff.data.values[scattcoeff.data.values == 0] = _np.nan
        scattcoeff.data *= 1e6  # conversion from inverse meter to mega meter

        if _np.any(abs_ratio):
            scattcoeff = scattcoeff + scattcoeff_bc
        return scattcoeff

    def calculate_new(self):
        self._calculate_all(False)

    def calculate_all(self):
        self._calculate_all(True)

    def _calculate_all(self, overwrite, time_window=False, verbose=False):

        bc_ratio = None
        perscribed_absorption = None
        bc_name_signitur = ''
        if self.black_carbon:
            if 'bc_ratio' in self.black_carbon.keys():
                bc_ratio = self.black_carbon['bc_ratio']
                bc_name_signitur = ('_bc_ratio_%.3f'%bc_ratio).replace('.', 'o')
            elif 'perscribed_absorption' in self.black_carbon.keys():
                perscribed_absorption = self.black_carbon['perscribed_absorption']
                if type(perscribed_absorption).__name__ == 'str':
                    pat = perscribed_absorption
                elif type(perscribed_absorption).__name__ in ['int', 'float']:
                    pat = '%.1f'%perscribed_absorption
                bc_name_signitur = ('_bc_abs_%s' %pat).replace('.', 'o')

        if self.mode_sel:
            if 'aosacsm' in [self.mode_sel['n_accu'], self.mode_sel['n_coarse']]:
                self.refractive_index = 'aosacsm'

        if self.keep_data:
            scattcoeff_list = []

        all_files = _os.listdir(self.folder + 'tdmaaps/')
        all_files = _np.array(all_files)
        all_files_tdmaapssize = all_files[_np.char.find(all_files, 'tdmaapssize') > -1]
        all_files_acsm = _np.array(_os.listdir(self.folder + 'acsm/'))
        test_done = False
        for e, fname_tdmaapssize in enumerate(all_files_tdmaapssize):
            if self.test:
                if fname_tdmaapssize == self.test_file:
                    if self.verbose:
                        print('found test file {}'.format(self.test_file))
                    test_done = True
                else:
                    if test_done:
                        break
                    else:
                        continue
            if time_window:
                if not _atm_arm._is_in_time_window(fname_tdmaapssize, verbose):
                    if self.verbose:
                        print('{} not in time window'.format(fname_tdmaapssize))
                    continue

            splitname = _splitup_filename(fname_tdmaapssize)
            site = splitname['site']
            date = splitname['date']
            if self.mode_sel:
                com2str = lambda x: '{:0.2f}j{:0.3f}'.format(x.real, x.imag) if type(x) != str else x
                refidxname = 'modesel_{}_{}'.format(com2str(self.mode_sel['n_accu']), com2str(self.mode_sel['n_coarse']))
                if 'gf_coarse' in self.mode_sel.keys():
                    refidxname += '_gfc{:0.2f}'.format(self.mode_sel['gf_coarse'])
            else:
                refidxname = self.refractive_index

            name_addon = '_RI%s_%s_%snm' % (refidxname, self.diameter_cutoff, self.wavelength)

            if self.apply_growth:
                name_addon += '_gf{:.2f}'.format(self.apply_growth)

            if self.data_quality != 'patchy':
                name_addon += '{}'.format(self.data_quality)


            name_addon = name_addon.replace('.', 'o')
            my_prod_name = self.folder_out + site + 'tdmaaps2scatteringcoeff' + bc_name_signitur + name_addon + '.' + date + '.000000.cdf'

            if not overwrite:
                if _os.path.isfile(my_prod_name):
                    if self.verbose:
                        print('product %s already exists' % my_prod_name)
                    continue

            if self.refractive_index == 'aosacsm':
                fname_others = _get_other_filenames(fname_tdmaapssize, ['aosacsm'], all_files_acsm)
                if not fname_others:
                    if self.verbose:
                        print('no matching acsm file found')
                    continue
                else:
                    aosacsm = _atm_arm.read_cdf(self.folder + 'acsm/' + fname_others['aosacsm']['fname'], data_quality=self.data_quality, verbose=verbose)
                    refractive_index = aosacsm.refractive_index
            elif type(self.refractive_index).__name__ == 'float':
                refractive_index = self.refractive_index

            if self.black_carbon:
                if 'perscribed_absorption' in self.black_carbon.keys():
                    if self.black_carbon['perscribed_absorption'] == 'noaaaos':
                        fname_others = _get_other_filenames(fname_tdmaapssize, ['noaaaos'], all_files, verbose= True)
                        if not fname_others:
                            continue
                        else:
                            noaaaos = _atm_arm.read_cdf(self.folder + 'noaaaos/' + fname_others['noaaaos']['fname'], data_quality=self.data_quality, verbose=verbose)
                            if self.wavelength == 550 and self.diameter_cutoff == '1um':
                                keep = 'Ba_G_Dry_1um_PSAP1W_1'
                            elif self.wavelength == 550 and self.diameter_cutoff == '10um':
                                keep = 'Ba_G_Dry_10um_PSAP1W_1'
                            else:
                                raise ValueError('Wavelength and diameter_cutoff do not match allowed values.\n\twavelength: %s\n\tdiameter_cutoff: %s'%(self.wavelength, self.diameter_cutoff))
                            perscribed_absorption = noaaaos.abs_coeff._del_all_columns_but(keep)
                            perscribed_absorption.data.rename(columns={keep: 'perscribed_absorption'}, inplace=True)


            tdmaapssize = _atm_arm.read_cdf(self.folder  + 'tdmaaps/' + fname_tdmaapssize, data_quality=self.data_quality, verbose=verbose)

            if self.mode_sel:
                dist = tdmaapssize.size_distribution
                n_accu = refractive_index if self.mode_sel['n_accu'] == 'aosacsm' else self.mode_sel['n_accu']
                n_coarse = refractive_index if self.mode_sel['n_coarse'] == 'aosacsm' else self.mode_sel['n_coarse']
                if self.mode_sel['moment'] ==  'volume':
                    dist = dist.convert2dVdlogDp()
                else:
                    raise KeyError('not possible')
                dist_accu = dist.mode_analysis.size_dist_accu
                dist_coarse = dist.mode_analysis.size_dist_coarse
                scattaccu = self._calculate_one(dist_accu, n_accu, bc_ratio, perscribed_absorption, None)
                if 'gf_coarse' in self.mode_sel.keys():
                    gf = self.mode_sel['gf_coarse']
                else:
                    gf = None
                scattcoarse = self._calculate_one(dist_coarse, n_coarse, bc_ratio, perscribed_absorption, gf)
                scattcoeff = scattaccu + scattcoarse
            else:
                scattcoeff = self._calculate_one(tdmaapssize.size_distribution, refractive_index, bc_ratio, perscribed_absorption, self.apply_growth)

            if self.keep_data:
                scattcoeff_list.append(scattcoeff)
            if not self.test:
                scattcoeff.save_netCDF(my_prod_name)
        # if len(extcoeff_list) == 2:
        #                 break


        print(my_prod_name)
        if self.keep_data:
            if len(scattcoeff_list) > 1:
                extcoeff_cat = _timeseries.concat(scattcoeff_list)
                self.result = extcoeff_cat.close_gaps(verbose=False)
            else:
                self.result = scattcoeff_list[0]

#########################
###### apply growth
#########################
#### to coarse mode only

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o60j0o010_gfc1o40_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index= None,
                                  mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                  'n_coarse': complex(1.6, 0.01), 'gf_coarse': 1.4},
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                  test=test,
                                  keep_data=keep_data,
                                  verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o60j0o010_gfc1o30_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index= None,
                                  mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                  'n_coarse': complex(1.6, 0.01), 'gf_coarse': 1.3},
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                  test=test,
                                  keep_data=keep_data,
                                  verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o60j0o010_gfc1o20_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index= None,
                                  mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                  'n_coarse': complex(1.6, 0.01), 'gf_coarse': 1.2},
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                  test=test,
                                  keep_data=keep_data,
                                  verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o60j0o010_gfc1o10_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index= None,
                                  mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                  'n_coarse': complex(1.6, 0.01), 'gf_coarse': 1.1},
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                  test=test,
                                  keep_data=keep_data,
                                  verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o60j0o010_gfc0o70_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index= None,
                                  mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                  'n_coarse': complex(1.6, 0.01), 'gf_coarse': 0.7},
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                  test=test,
                                  keep_data=keep_data,
                                  verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o60j0o010_gfc0o90_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index= None,
                                  mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                  'n_coarse': complex(1.6, 0.01), 'gf_coarse': 0.9},
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                  test=test,
                                  keep_data=keep_data,
                                  verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o60j0o010_gfc0o80_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index= None,
                                  mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                  'n_coarse': complex(1.6, 0.01), 'gf_coarse': 0.8},
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                  test=test,
                                  keep_data=keep_data,
                                  verbose=verbose)
    return out

#### overall


def tdmaaps2scatteringcoeff_RI1o5_10um_550nm_gf0o70(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index= 1.5,
                                  apply_growth = 0.7,
#                                               mode_selection= {'moment':'volume', 'n_accu':'aosacsm', 'n_coarse': complex(1.6, 0.01)},
#                                               mode_selection= {'moment':'volume', 'n_accu':'acsm', 'n_coarse': complex(1.63, 0.01)},
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                  test=test,
                                  keep_data=keep_data,
                                  verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RI1o5_10um_550nm_gf0o80(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index= 1.5,
                                  apply_growth = 0.8,
#                                               mode_selection= {'moment':'volume', 'n_accu':'aosacsm', 'n_coarse': complex(1.6, 0.01)},
#                                               mode_selection= {'moment':'volume', 'n_accu':'acsm', 'n_coarse': complex(1.63, 0.01)},
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                  test=test,
                                  keep_data=keep_data,
                                  verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RI1o5_10um_550nm_gf0o90(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index= 1.5,
                                  apply_growth = 0.9,
#                                               mode_selection= {'moment':'volume', 'n_accu':'aosacsm', 'n_coarse': complex(1.6, 0.01)},
#                                               mode_selection= {'moment':'volume', 'n_accu':'acsm', 'n_coarse': complex(1.63, 0.01)},
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                  test=test,
                                  keep_data=keep_data,
                                  verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RI1o5_10um_550nm_gf1o10(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index= 1.5,
                                  apply_growth = 1.1,
#                                               mode_selection= {'moment':'volume', 'n_accu':'aosacsm', 'n_coarse': complex(1.6, 0.01)},
#                                               mode_selection= {'moment':'volume', 'n_accu':'acsm', 'n_coarse': complex(1.63, 0.01)},
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                  test=test,
                                  keep_data=keep_data,
                                  verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RI1o5_10um_550nm_gf1o20(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index= 1.5,
                                  apply_growth = 1.2,
#                                               mode_selection= {'moment':'volume', 'n_accu':'aosacsm', 'n_coarse': complex(1.6, 0.01)},
#                                               mode_selection= {'moment':'volume', 'n_accu':'acsm', 'n_coarse': complex(1.63, 0.01)},
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                  test=test,
                                  keep_data=keep_data,
                                  verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RI1o5_10um_550nm_gf1o30(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index= 1.5,
                                  apply_growth = 1.3,
#                                               mode_selection= {'moment':'volume', 'n_accu':'aosacsm', 'n_coarse': complex(1.6, 0.01)},
#                                               mode_selection= {'moment':'volume', 'n_accu':'acsm', 'n_coarse': complex(1.63, 0.01)},
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                  test=test,
                                  keep_data=keep_data,
                                  verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RI1o5_10um_550nm_gf1o40(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index= 1.5,
                                  apply_growth = 1.4,
#                                               mode_selection= {'moment':'volume', 'n_accu':'aosacsm', 'n_coarse': complex(1.6, 0.01)},
#                                               mode_selection= {'moment':'volume', 'n_accu':'acsm', 'n_coarse': complex(1.63, 0.01)},
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                  test=test,
                                  keep_data=keep_data,
                                  verbose=verbose)
    return out
###################
### mode sel all n fixed
###################

def tdmaaps2scatteringcoeff_RImodesel_1o50j0o00_1o63j0o01_1um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                      diameter_cutoff='1um',
                                      wavelength=550,
                                      refractive_index=None,
                                      mode_selection={'moment': 'volume', 'n_accu': 1.5,
                                                      'n_coarse': complex(1.63, 0.01)},
                                      folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                      folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                      test=test,
                                      keep_data=keep_data,
                                      verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_1o50j0o00_1o63j0o01_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                      diameter_cutoff='10um',
                                      wavelength=550,
                                      refractive_index=None,
                                      mode_selection={'moment': 'volume', 'n_accu': 1.5,
                                                      'n_coarse': complex(1.63, 0.01)},
                                      folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                      folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                      test=test,
                                      keep_data=keep_data,
                                      verbose=verbose)
    return out

########
### coarse mode n dependence
#######

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o66j0o01_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                      diameter_cutoff='10um',
                                      wavelength=550,
                                      refractive_index=None,
                                      mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                      'n_coarse': complex(1.66, 0.01)},
                                      folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                      folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                      test=test,
                                      keep_data=keep_data,
                                      verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o57j0o01_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                      diameter_cutoff='10um',
                                      wavelength=550,
                                      refractive_index=None,
                                      mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                      'n_coarse': complex(1.57, 0.01)},
                                      folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                      folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                      test=test,
                                      keep_data=keep_data,
                                      verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o54j0o01_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                      diameter_cutoff='10um',
                                      wavelength=550,
                                      refractive_index=None,
                                      mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                      'n_coarse': complex(1.54, 0.01)},
                                      folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                      folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                      test=test,
                                      keep_data=keep_data,
                                      verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o60j0o01_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                      diameter_cutoff='10um',
                                      wavelength=550,
                                      refractive_index=None,
                                      mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                      'n_coarse': complex(1.6, 0.01)},
                                      folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                      folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                      test=test,
                                      keep_data=keep_data,
                                      verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o63j0o01_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                      diameter_cutoff='10um',
                                      wavelength=550,
                                      refractive_index=None,
                                      mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                      'n_coarse': complex(1.63, 0.01)},
                                      folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                      folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                      test=test,
                                      keep_data=keep_data,
                                      verbose=verbose)
    return out

############
### Coarse mode imaginary part of n dependence
##########

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o63j0o00_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                      diameter_cutoff='10um',
                                      wavelength=550,
                                      refractive_index=None,
                                      mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                      'n_coarse': 1.63},
                                      folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                      folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                      test=test,
                                      keep_data=keep_data,
                                      verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o63j0o005_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                      diameter_cutoff='10um',
                                      wavelength=550,
                                      refractive_index=None,
                                      mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                      'n_coarse': complex(1.63, 0.005)},
                                      folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                      folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                      test=test,
                                      keep_data=keep_data,
                                      verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o63j0o015_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                      diameter_cutoff='10um',
                                      wavelength=550,
                                      refractive_index=None,
                                      mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                      'n_coarse': complex(1.63, 0.015)},
                                      folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                      folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                      test=test,
                                      keep_data=keep_data,
                                      verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o63j0o025_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                      diameter_cutoff='10um',
                                      wavelength=550,
                                      refractive_index=None,
                                      mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                      'n_coarse': complex(1.63, 0.025)},
                                      folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                      folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                      test=test,
                                      keep_data=keep_data,
                                      verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o63j0o02_10um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                      diameter_cutoff='10um',
                                      wavelength=550,
                                      refractive_index=None,
                                      mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                      'n_coarse': complex(1.63, 0.02)},
                                      folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                      folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                      test=test,
                                      keep_data=keep_data,
                                      verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RImodesel_aosacsm_1o63j0o01_1um_550nm(test = False):
    keep_data = verbose = True if test else False
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                      diameter_cutoff='1um',
                                      wavelength=550,
                                      refractive_index=None,
                                      mode_selection={'moment': 'volume', 'n_accu': 'aosacsm',
                                                      'n_coarse': complex(1.63, 0.01)},
                                      folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                      folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                      test=test,
                                      keep_data=keep_data,
                                      verbose=verbose)
    return out

def tdmaaps2scatteringcoeff_RIaosacsm_1um_550nm(test = False):
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='1um',
                                  wavelength=550,
                                  refractive_index='aosacsm',
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Users/htelg/data/ARM/SGP/',
                                  test = test)
    return out

def tdmaaps2scatteringcoeff_RIaosacsm_10um_550nm(test = False):
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index='aosacsm',
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Users/htelg/data/ARM/SGP/',
                                  test = test)
    return out

def tdmaaps2scatteringcoeff_RIaosacsm_1um_550nm_good():
    out = Tdmaaps2scatteringcoeff(data_quality='good',
                                  diameter_cutoff='1um',
                                  wavelength=550,
                                  refractive_index='aosacsm',
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Users/htelg/data/ARM/SGP/',
                                  test = False)
    return out

def tdmaaps2scatteringcoeff_RIaosacsm_1um_550nm_bad():
    out = Tdmaaps2scatteringcoeff(data_quality='bad',
                                  diameter_cutoff='1um',
                                  wavelength=550,
                                  refractive_index='aosacsm',
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Users/htelg/data/ARM/SGP/',
                                  test = False)
    return out

def tdmaaps2scatteringcoeff_RI1o5_1um_550nm(test = False):
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='1um',
                                  wavelength=550,
                                  refractive_index=1.5,
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Users/htelg/data/ARM/SGP/',
                                  test = test)
    return out


def tdmaaps2scatteringcoeff_RI1o5_10um_550nm(test = False):
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index=1.5,
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Users/htelg/data/ARM/SGP/',
                                  test = test)
    return out


def tdmaaps2scatteringcoeff_bc_ratio_0o050_RI1o5_10um_550nm(test = False):
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                 diameter_cutoff='10um',
                                 wavelength=550,
                                 refractive_index=1.5,
                                 black_carbon={'bc_ratio': 0.05},
                                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                 folder='/Volumes/HTelg_4TB_Backup/arm_data/SGP/',
                                 test=test)
    return out

def tdmaaps2scatteringcoeff_bc_abs_noaaaos_RI1o5_10um_550nm(test = False):
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                 diameter_cutoff='10um',
                                 wavelength=550,
                                 refractive_index=1.5,
                                 black_carbon={'perscribed_absorption': 'noaaaos'},
                                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                 folder='/Users/htelg/data/ARM/SGP/',
                                 test=test)
    return out

def tdmaaps2scatteringcoeff_bc_abs_noaaaos_RIaosacsm_10um_550nm(test = False):
    out = Tdmaaps2scatteringcoeff(data_quality='patchy',
                                 diameter_cutoff='10um',
                                 wavelength=550,
                                  refractive_index='aosacsm',
                                  black_carbon={'perscribed_absorption': 'noaaaos'},
                                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                 folder='/Users/htelg/data/ARM/SGP/',
                                 test=test)
    return out


####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

class Tdmaaps2backscatteringratio(object):
    def __init__(self, data_quality='patchy',
                 diameter_cutoff='1um',
                 wavelength=550,  # in nm
                 refractive_index = 'aosacsm',
                 black_carbon = False,
                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                 folder='/Users/htelg/data/ARM/SGP/',
                 keep_data = False,
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
        absorption: dictionary
            'perscribed_absorption' in Mm^-1
            'abs_ratio'
            dict: 'noaaaos', will use the absorption measurements from the noaaaos product
            float: fraction to assume to be black carbon, e.g. 0.05 will assume 5% of the accumulation mode to be black carbon
        folder_out
        keep_data: bool
            if the data is returned or discarded after saving it. Returning it can result in massive RAM use.
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
        if self.test:
            keep_data = True
        self.test_file = 'sgptdmaapssizeC1.c1.20120201.002958.cdf'
        self.keep_data = keep_data

        self.black_carbon = black_carbon
        if self.black_carbon:
            allowed_kwargs_bc = {'perscribed_absorption': [],
                                 'bc_ratio': []}

            for k in self.black_carbon:
                if k not in allowed_kwargs_bc.keys():
                    txt = '%s is an unknown key for kwarg black_carbon. Allowed keys are %s.'%(k, allowed_kwargs_bc.keys())
                    raise KeyError(txt)

        if test:
            self.intres = {}

    def _calculate_one(self, tdmaapssize, refractive_index, abs_ratio, perscribed_absorption):
        """

        Parameters
        ----------
        tdmaapssize: obvious
        refractive_index: obvious
        abs_ratio: ratio of the numberconcentration that is black carbon
        perscribed_absorption: this will result in an iteration in which the abs_ratio will be determined that is necessary so the absorbtion is == abs_goal

        Returns
        -------

        """
        if self.diameter_cutoff == '1um':
            dcoff = 1000
        elif self.diameter_cutoff == '10um':
            dcoff = 10000

        dist = tdmaapssize.size_distribution.zoom_diameter(end=dcoff)
        dist = dist.convert2numberconcentration()

        if _np.any(perscribed_absorption):
            dist_bc = dist.copy()
            dist_bc.index_of_refraction = complex(2.26, 1.26)
            opt = dist_bc.calculate_optical_properties(self.wavelength)
            absorption_coeff_bc = opt.absorption_coeff
            absorption_coeff_bc.data.values[absorption_coeff_bc.data.values == 0] = _np.nan
            absorption_coeff_bc.data *= 1e6  # conversion from inverse meter to mega meter

            if type(perscribed_absorption).__name__ == 'TimeSeries':
                abs_ratio = (perscribed_absorption / absorption_coeff_bc).data['perscribed_absorption']
            else:
                abs_ratio = (perscribed_absorption / absorption_coeff_bc.data['abs_coeff_m^1'])


        if _np.any(abs_ratio):
            dist_bc = dist.copy()
            dist_bc.data = dist_bc.data.mul(abs_ratio, axis = 0)
            dist_bc.index_of_refraction = complex(2.26,1.26)
            opt = dist_bc.calculate_optical_properties(self.wavelength)
            backscatt_bc = opt.hemispheric_backscattering_ratio
            backscatt_bc.data.values[backscatt_bc.data.values == 0] = _np.nan
        else:
            abs_ratio = 0

        dist.data = dist.data.mul((1 - abs_ratio), axis = 0)

        dist.index_of_refraction = refractive_index
        opt = dist.calculate_optical_properties(self.wavelength)
        backscatt = opt.hemispheric_backscattering_ratio
        backscatt.data.values[backscatt.data.values == 0] = _np.nan

        if _np.any(abs_ratio):
            backscatt_bc.data = backscatt_bc.data.multiply(abs_ratio, axis=0)
            backscatt.data = backscatt.data.multiply((1 - abs_ratio), axis=0)
            backscatt = backscatt + backscatt_bc
        return backscatt

    def calculate_new(self):
        self._calculate_all(False)

    def calculate_all(self):
        self._calculate_all(True)

    def _calculate_all(self, overwrite, time_window=False, verbose=False):

        bc_ratio = None
        perscribed_absorption = None
        bc_name_signitur = ''
        if self.black_carbon:
            if 'bc_ratio' in self.black_carbon.keys():
                bc_ratio = self.black_carbon['bc_ratio']
                bc_name_signitur = ('_bc_ratio_%.3f'%bc_ratio).replace('.', 'o')
            elif 'perscribed_absorption' in self.black_carbon.keys():
                perscribed_absorption = self.black_carbon['perscribed_absorption']
                if type(perscribed_absorption).__name__ == 'str':
                    pat = perscribed_absorption
                elif type(perscribed_absorption).__name__ in ['int', 'float']:
                    pat = '%.1f'%perscribed_absorption
                bc_name_signitur = ('_bc_abs_%s' %pat).replace('.', 'o')

        if self.keep_data:
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
            if self.data_quality == 'patchy':
                name_addon = ('_RI%s_%s_%snm' % (self.refractive_index, self.diameter_cutoff, self.wavelength)).replace('.','o')
            else:
                name_addon = ('_RI%s_%s_%snm_%s' % (self.refractive_index, self.diameter_cutoff, self.wavelength, self.data_quality)).replace('.', 'o')
            my_prod_name = self.folder_out + site + 'tdmaaps2backscatteringratio' + bc_name_signitur + name_addon + '.' + date + '.000000.cdf'

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

            if self.black_carbon:
                if 'perscribed_absorption' in self.black_carbon.keys():
                    if self.black_carbon['perscribed_absorption'] == 'noaaaos':
                        fname_others = _get_other_filenames(fname_tdmaapssize, ['noaaaos'], all_files, verbose= True)
                        if not fname_others:
                            continue
                        else:
                            noaaaos = _atm_arm.read_cdf(self.folder + fname_others['noaaaos']['fname'], data_quality=self.data_quality, verbose=verbose)
                            if self.wavelength == 550 and self.diameter_cutoff == '1um':
                                keep = 'Ba_G_Dry_1um_PSAP1W_1'
                            elif self.wavelength == 550 and self.diameter_cutoff == '10um':
                                keep = 'Ba_G_Dry_10um_PSAP1W_1'
                            else:
                                raise ValueError('Wavelength and diameter_cutoff do not match allowed values.\n\twavelength: %s\n\tdiameter_cutoff: %s'%(self.wavelength, self.diameter_cutoff))
                            perscribed_absorption = noaaaos.abs_coeff._del_all_columns_but(keep)
                            perscribed_absorption.data.rename(columns={keep: 'perscribed_absorption'}, inplace=True)


            tdmaapssize = _atm_arm.read_cdf(self.folder + fname_tdmaapssize, data_quality=self.data_quality, verbose=verbose)

            extcoeff = self._calculate_one(tdmaapssize, refractive_index, bc_ratio, perscribed_absorption)

            if self.keep_data:
                extcoeff_list.append(extcoeff)
            if not self.test:
                extcoeff.save_netCDF(my_prod_name)
        # if len(extcoeff_list) == 2:
        #                 break


        print(my_prod_name)
        if self.keep_data:
            if len(extcoeff_list) > 1:
                extcoeff_cat = _timeseries.concat(extcoeff_list)
                self.result = extcoeff_cat.close_gaps(verbose=False)
            else:
                self.result = extcoeff_list[0]

def tdmaaps2backscatteringratio_RI1o5_1um_550nm(test=False):
    out = Tdmaaps2backscatteringratio(data_quality='patchy',
                                  diameter_cutoff='1um',
                                  wavelength=550,
                                  refractive_index=1.5,
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Users/htelg/data/ARM/SGP/',
                                  test=test)
    return out

def tdmaaps2backscatteringratio_RI1o5_10um_550nm(test = False):
    out = Tdmaaps2backscatteringratio(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index=1.5,
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Users/htelg/data/ARM/SGP/',
                                  test = test)
    return out

def tdmaaps2backscatteringratio_RIaosacsm_10um_550nm(test = False):
    out = Tdmaaps2backscatteringratio(data_quality='patchy',
                                  diameter_cutoff='10um',
                                  wavelength=550,
                                  refractive_index='aosacsm',
                                  folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                  folder='/Users/htelg/data/ARM/SGP/',
                                  test = test)
    return out

def tdmaaps2backscatteringratio_bc_abs_noaaaos_RI1o5_10um_550nm(test = False):
    out = Tdmaaps2backscatteringratio(data_quality='patchy',
                                 diameter_cutoff='10um',
                                 wavelength=550,
                                 refractive_index=1.5,
                                 black_carbon={'perscribed_absorption': 'noaaaos'},
                                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                 folder='/Users/htelg/data/ARM/SGP/',
                                 test=test)
    return out

def tdmaaps2backscatteringratio_bc_abs_noaaaos_RIaosacsm_10um_550nm(test = False):
    out = Tdmaaps2backscatteringratio(data_quality='patchy',
                                 diameter_cutoff='10um',
                                 wavelength=550,
                                  refractive_index='aosacsm',
                                  black_carbon={'perscribed_absorption': 'noaaaos'},
                                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                 folder='/Users/htelg/data/ARM/SGP/',
                                 test=test)
    return out
##########################################################################################
##########################################################################################
##########################################################################################


class Aipfitrh2Kappa(object):
    def __init__(self, data_quality='patchy',
                 diameter_cutoff='1um',
                 wavelength=550,  # in nm
                 sizedistribution = 'tdmaapssize',
                 refractive_index = 1.5, #'aosacsm',
                 f_of_rh_product = 'f_RH_scatt_2p_85_40',
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
        self.sizedistribution = sizedistribution
        self.f_of_rh_product = f_of_rh_product
        self.data_quality = data_quality
        self.diameter_cutoff = diameter_cutoff
        self.folder_out = folder_out
        self.folder = folder
        self.test = test
        self.test_file = 'sgpaipfitrh1ogrenC1.c1.20120302.000000.cdf'

    def _calculate_one(self, aipfitrh, sizedist, refractive_index):
        if self.diameter_cutoff == '1um':
            dcoff = 1000
        elif self.diameter_cutoff == '10um':
            dcoff = 10000

        sizedist = sizedist.zoom_diameter(end = dcoff)
        sizedist.index_of_refraction = refractive_index

        aipfitrh.sup_kappa_sizedist = sizedist
        aipfitrh.sup_kappa_wavelength = self.wavelength

        if self.f_of_rh_product == 'f_RH_scatt_2p_85_40':
            aipfitrh.f_RH_scatt_2p_85_40._del_all_columns_but('ratio_85by40_Bs_G_1um_2p', inplace=True)
            return aipfitrh.kappa_85_40

        else:
            raise ValueError("Don't know what to do with self.f_of_rh_product = %s"%(self.f_of_rh_product))

    def calculate_new(self):
        self._calculate_all(False)

    def calculate_all(self):
        self._calculate_all(True)

    def _calculate_all(self, overwrite, time_window=False, verbose=False):

        if self.f_of_rh_product == 'f_RH_scatt_2p_85_40':
            RH_for_name = 'RH_85_40'

        kappa_list = []
        all_files = _os.listdir(self.folder)
        all_files = _np.array(all_files)

        all_files_aipfitrh = all_files[_np.char.find(all_files, 'aipfitrh1ogrenC1') > -1]
        test_done = False

        for e, fname_aipfitrh in enumerate(all_files_aipfitrh):
            if self.test:
                if fname_aipfitrh == self.test_file:
                    test_done = True
                else:
                    if test_done:
                        break
                    else:
                        continue
            if time_window:
                if not _atm_arm._is_in_time_window(fname_aipfitrh, verbose):
                    continue

            splitname = _splitup_filename(fname_aipfitrh)
            site = splitname['site']
            date = splitname['date']

            name_addon = ('%s_%s_RI%s_%s_%snm_%s' % (RH_for_name, self.sizedistribution, self.refractive_index, self.diameter_cutoff, self.wavelength, self.data_quality)).replace('.', 'o')
            my_prod_name = self.folder_out + site + 'aipfitrh2kappa_' + name_addon + '.' + date + '.000000.cdf'

            if not overwrite:
                if _os.path.isfile(my_prod_name):
                    if verbose:
                        print('product %s already exists' % my_prod_name)
                    continue

            if self.sizedistribution == 'tdmaapssize':
                fname_others = _get_other_filenames(fname_aipfitrh, ['tdmaapssize'], all_files)
                if not fname_others:
                    continue
                else:
                    tdmaaps = _atm_arm.read_cdf(self.folder + fname_others['tdmaapssize']['fname'], data_quality=self.data_quality, verbose=verbose)
                    sizedist = tdmaaps.size_distribution

            else:
                txt = "Unknown sizedistribution type (%s). Try one of these: 'tdmaapssize'"%(self.sizedistribution)
                raise ValueError(txt)

            if self.refractive_index == 'aosacsm':
                fname_others = _get_other_filenames(fname_aipfitrh, ['aosacsm'], all_files)
                if not fname_others:
                    continue
                else:
                    aosacsm = _atm_arm.read_cdf(self.folder + fname_others['aosacsm']['fname'], data_quality=self.data_quality, verbose=verbose)
                    refractive_index = aosacsm.refractive_index
            elif type(self.refractive_index).__name__ == 'float':
                refractive_index = self.refractive_index

            aipfitrh = _atm_arm.read_cdf(self.folder + fname_aipfitrh, data_quality=self.data_quality, verbose=verbose)

            kappa = self._calculate_one(aipfitrh, sizedist, refractive_index)

            kappa_list.append(kappa)
            if not self.test:
                kappa.save_netCDF(my_prod_name)
        # if len(extcoeff_list) == 2:
        #                 break


        print(my_prod_name)
        if len(kappa_list) > 1:
            extcoeff_cat = _timeseries.concat(kappa_list)
            self.result = extcoeff_cat.close_gaps(verbose=False)
        else:
            self.result = kappa_list[0]


def sgpaipfitrh2kappa_RH_85_40_tdmaapssize_RI1o5_1um_550nm_patchy(test = False):
    product = Aipfitrh2Kappa(data_quality='patchy',
                             diameter_cutoff='1um',
                             wavelength=550,
                             sizedistribution='tdmaapssize',
                             refractive_index=1.5,
                             f_of_rh_product='f_RH_scatt_2p_85_40',
                             folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                             folder='/Users/htelg/data/ARM/SGP/',
                             test=test)
    return product




class Tdmahyg2Kappa(object):
    def __init__(self, data_quality='patchy',
                 method = 'avg',
                 diameter = 200,
                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                 folder='/Users/htelg/data/ARM/SGP/',
                 test = False):
        """

        Parameters
        ----------
        data_quality
        method = 'avg' #methode of how kapppa was calculated, 'avg' simply take the average of gf
        diameter = 200 #center diameter for which the kappa is calculated for

        folder_out
        folder
        test
        """
        self.data_quality = data_quality
        self.method = method
        self.diameter = diameter
        self.folder_out = folder_out
        self.folder = folder
        self.test = test
        self.test_file = 'sgptdmahygC1.b1.20120301.003735.cdf'

    def _calculate_one(self, tdmahyg):
        kappas = tdmahyg.kappa_values._del_all_columns_but(self.diameter)
        kappas1D = _timeseries.TimeSeries(kappas.data)
        kappas1D._data_period = kappas._data_period
        return kappas1D

    def calculate_new(self):
        self._calculate_all(False)

    def calculate_all(self):
        self._calculate_all(True)

    def _calculate_all(self, overwrite, time_window=False, verbose=False):


        kappa_list = []
        all_files = _os.listdir(self.folder)
        all_files = _np.array(all_files)

        all_files_tdmahyg = all_files[_np.char.find(all_files, 'tdmahyg') > -1]
        test_done = False

        for e, fname_tdmahyg in enumerate(all_files_tdmahyg):
            if self.test:
                if fname_tdmahyg == self.test_file:
                    test_done = True
                else:
                    if test_done:
                        break
                    else:
                        continue
            if time_window:
                if not _atm_arm._is_in_time_window(fname_tdmahyg, verbose):
                    continue

            splitname = _splitup_filename(fname_tdmahyg)
            site = splitname['site']
            date = splitname['date']

            name_addon = ('%s_d%s_%s' % (self.method, self.diameter, self.data_quality)).replace('.', 'o')
            my_prod_name = self.folder_out + site + 'tdmahyg2kappa_' + name_addon + '.' + date + '.000000.cdf'

            if not overwrite:
                if _os.path.isfile(my_prod_name):
                    if verbose:
                        print('product %s already exists' % my_prod_name)
                    continue

            tdmahyg = _atm_arm.read_cdf(self.folder + fname_tdmahyg, data_quality=self.data_quality, verbose=verbose)

            kappa = self._calculate_one(tdmahyg)

            kappa_list.append(kappa)
            if not self.test:
                kappa.save_netCDF(my_prod_name)
        # if len(extcoeff_list) == 2:
        #                 break


        print(my_prod_name)
        if len(kappa_list) > 1:
            extcoeff_cat = _timeseries.concat(kappa_list)
            self.result = extcoeff_cat.close_gaps(verbose=False)
        else:
            self.result = kappa_list[0]

# def sgptdmahyg2kappa_avg_d200_patchy(test = False):
#     prod = Tdmahyg2Kappa(data_quality='patchy',
#                                      method='avg',
#                                      diameter=200,
#                                      folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
#                                      folder='/Users/htelg/data/ARM/SGP/',
#                                      test=test)
#     return prod


class Noaaaos2Hygroscopicity(object):
    def __init__(self,
                 data_quality='patchy',
                 diameter_cutoff='10um',
                 wavelength=550,  # in nm
                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                 folder='/Users/htelg/data/ARM/SGP/',
                 keep_data = False,
                 test = False,
                 error_bad_file = False,
                 verbose = False):
        """

        Parameters
        ----------
        data_quality
        diameter_cutoff
        wavelength
        refractive_index: string or float
            Define which data product to get the RI from ['aosacsm']
            or set a fixed value.
        absorption: dictionary
            'perscribed_absorption' in Mm^-1
            'abs_ratio'
            dict: 'noaaaos', will use the absorption measurements from the noaaaos product
            float: fraction to assume to be black carbon, e.g. 0.05 will assume 5% of the accumulation mode to be black carbon
        folder_out
        keep_data: bool
            if the data is returned or discarded after saving it. Returning it can result in massive RAM use.
        folder
        test
        """
        self.wavelength = wavelength
        self.data_quality = data_quality
        self.diameter_cutoff = diameter_cutoff
        self.folder_out = folder_out
        self.folder = folder
        self.test = test
        if self.test:
            keep_data = True

        if type(test).__name__ == 'str':
            self.test_file = test
        else:
            self.test_file = 'sgpnoaaaosC1.b1.20120201.000000.cdf'
        self.keep_data = keep_data

        if test:
            self.intres = {}

        self.verbose = verbose
        self._error_bad_file = error_bad_file

    def _calculate_one(self, noaaaos):
        """

        Parameters
        ----------
        tdmaapssize: obvious
        refractive_index: obvious
        abs_ratio: ratio of the numberconcentration that is black carbon
        perscribed_absorption: this will result in an iteration in which the abs_ratio will be determined that is necessary so the absorbtion is == abs_goal

        Returns
        -------

        """
        if self.diameter_cutoff == '1um':
            out = noaaaos.hygroscopicity_1um
        elif self.diameter_cutoff == '10um':
            out = noaaaos.hygroscopicity_10um
        else:
            raise ValueError('not implemented yet ... programming required')


        return out

    def calculate_new(self):
        self._calculate_all(False)

    def calculate_all(self):
        self._calculate_all(True)

    def _calculate_all(self, overwrite, time_window=False, verbose=False):

        if self.keep_data:
            hygro_list = []

        all_files = _os.listdir(self.folder)
        all_files = _np.array(all_files)
        all_files_noaaaos = all_files[_np.char.find(all_files, 'noaaaos') > -1]
        test_done = False
        for e, fname_noaaaos in enumerate(all_files_noaaaos):
            if self.verbose:
                print('===========')
                print('file: ', fname_noaaaos)
            if self.test:
                if fname_noaaaos == self.test_file:
                    test_done = True
                else:
                    if test_done:
                        break
                    else:
                        continue
            if time_window:
                if not _atm_arm._is_in_time_window(fname_noaaaos, verbose):
                    continue

            splitname = _splitup_filename(fname_noaaaos)
            site = splitname['site']
            date = splitname['date']
            name_addon = ('_%s_%snm_%s' % (self.diameter_cutoff, self.wavelength, self.data_quality)).replace('.', 'o')
            my_prod_name = self.folder_out + site + 'noaaaos2hygroscopicity' + name_addon + '.' + date + '.000000.cdf'

            if not overwrite:
                if _os.path.isfile(my_prod_name):
                    if verbose:
                        print('product %s already exists' % my_prod_name)
                    continue


            noaaaos = _atm_arm.read_cdf(self.folder + fname_noaaaos, data_quality=self.data_quality, verbose=verbose, error_bad_file = self._error_bad_file)
            if noaaaos._parsing_error:
                if verbose:
                    print('Parsing error, file skipped!')
                continue

            hygroscopicity = self._calculate_one(noaaaos)

            if self.keep_data:
                hygro_list.append(hygroscopicity)
            if not self.test:
                hygroscopicity.save_netCDF(my_prod_name)
        # if len(extcoeff_list) == 2:
        #                 break


        print(my_prod_name)
        if self.keep_data:
            if len(hygro_list) > 1:
                hygro_cat = _timeseries.concat(hygro_list)
                self.result = hygro_cat.close_gaps(verbose=False)
            else:
                self.result = hygro_list[0]

def noaaaos2hygroscopicity_1um_550nm_patchy(test = False, verbose = False):
    out = Noaaaos2Hygroscopicity(data_quality='patchy',
                                   diameter_cutoff='1um',
                                   wavelength=550,
                                   folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                   folder='/Users/htelg/data/ARM/SGP/',
                                   keep_data=False,
                                   test=test,
                                   error_bad_file=False,
                                   verbose=verbose)
    return out


def noaaaos2hygroscopicity_10um_550nm_patchy(test = False, verbose = False):
    out = Noaaaos2Hygroscopicity(data_quality='patchy',
                                 diameter_cutoff='10um',
                                 wavelength=550,
                                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                 folder='/Users/htelg/data/ARM/SGP/',
                                 keep_data=False,
                                 test=test,
                                 verbose = verbose)
    return out


# todo: below is not working yet!!!
class Acsm2RefrectiveIndex(object):
    def __init__(self,
                 data_quality='patchy',
                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                 folder='/Users/htelg/data/ARM/SGP/',
                 keep_data = False,
                 test = False,
                 error_bad_file = False,
                 verbose = False):
        """

        Parameters
        ----------
        data_quality
        diameter_cutoff
        wavelength
        refractive_index: string or float
            Define which data product to get the RI from ['aosacsm']
            or set a fixed value.
        absorption: dictionary
            'perscribed_absorption' in Mm^-1
            'abs_ratio'
            dict: 'noaaaos', will use the absorption measurements from the noaaaos product
            float: fraction to assume to be black carbon, e.g. 0.05 will assume 5% of the accumulation mode to be black carbon
        folder_out
        keep_data: bool
            if the data is returned or discarded after saving it. Returning it can result in massive RAM use.
        folder
        test
        """
        self.data_quality = data_quality
        self.folder_out = folder_out
        self.folder = folder
        self.test = test
        if self.test:
            keep_data = True

        if type(test).__name__ == 'str':
            self.test_file = test
        else:
            self.test_file = 'sgpnoaaaosC1.b1.20120201.000000.cdf'
        self.keep_data = keep_data

        if test:
            self.intres = {}

        self.verbose = verbose
        self._error_bad_file = error_bad_file

    def _calculate_one(self, noaaaos):
        """

        Parameters
        ----------
        tdmaapssize: obvious
        refractive_index: obvious
        abs_ratio: ratio of the numberconcentration that is black carbon
        perscribed_absorption: this will result in an iteration in which the abs_ratio will be determined that is necessary so the absorbtion is == abs_goal

        Returns
        -------

        """
        # if self.diameter_cutoff == '1um':
        #     raise ValueError (value err)
        if self.diameter_cutoff == '10um':
            dcoff = 10000
        else:
            raise ValueError('not implemented yet ... programming required')


        return noaaaos.hygroscopicity

    def calculate_new(self):
        self._calculate_all(False)

    def calculate_all(self):
        self._calculate_all(True)

    def _calculate_all(self, overwrite, time_window=False, verbose=False):

        if self.keep_data:
            hygro_list = []

        all_files = _os.listdir(self.folder)
        all_files = _np.array(all_files)
        all_files_noaaaos = all_files[_np.char.find(all_files, 'noaaaos') > -1]
        test_done = False
        for e, fname_noaaaos in enumerate(all_files_noaaaos):
            if self.verbose:
                print('===========')
                print('file: ', fname_noaaaos)
            if self.test:
                if fname_noaaaos == self.test_file:
                    test_done = True
                else:
                    if test_done:
                        break
                    else:
                        continue
            if time_window:
                if not _atm_arm._is_in_time_window(fname_noaaaos, verbose):
                    continue

            splitname = _splitup_filename(fname_noaaaos)
            site = splitname['site']
            date = splitname['date']
            name_addon = ('_%s_%snm_%s' % (self.diameter_cutoff, self.wavelength, self.data_quality)).replace('.', 'o')
            my_prod_name = self.folder_out + site + 'noaaaos2hygroscopicity' + name_addon + '.' + date + '.000000.cdf'

            if not overwrite:
                if _os.path.isfile(my_prod_name):
                    if verbose:
                        print('product %s already exists' % my_prod_name)
                    continue


            noaaaos = _atm_arm.read_cdf(self.folder + fname_noaaaos, data_quality=self.data_quality, verbose=verbose, error_bad_file = self._error_bad_file)
            if noaaaos._parsing_error:
                if verbose:
                    print('Parsing error, file skipped!')
                continue

            hygroscopicity = self._calculate_one(noaaaos)

            if self.keep_data:
                hygro_list.append(hygroscopicity)
            if not self.test:
                hygroscopicity.save_netCDF(my_prod_name)
        # if len(extcoeff_list) == 2:
        #                 break


        print(my_prod_name)
        if self.keep_data:
            if len(hygro_list) > 1:
                hygro_cat = _timeseries.concat(hygro_list)
                self.result = hygro_cat.close_gaps(verbose=False)
            else:
                self.result = hygro_list[0]


def acsm2refrectiveindex(test = False, verbose = False):
    out = Noaaaos2Hygroscopicity(data_quality='patchy',
                                 diameter_cutoff='10um',
                                 wavelength=550,
                                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                                 folder='/Users/htelg/data/ARM/SGP/',
                                 keep_data=False,
                                 test=test,
                                 verbose = verbose)
    return out
