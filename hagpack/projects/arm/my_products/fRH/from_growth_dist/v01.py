import numpy as _np
import os as _os
from atmPy.data_archives.arm import _read_data as _atm_arm
from atmPy.general import timeseries as _timeseries

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


class Product(object):
    """Growth distributionsn (tdmahyg) are applied to size distributions (tdmaaps) where each growth mode is applied
    separately according to its fraction of the entire size distribution. Optical prperties are calculated for each
    separately grown sizedistribution fraction and then add up.

    ToDo
    ----
    separate accumulation from coarse mode and treat separately ... maybe take difference between 200 and 600 growthmode
    as the growthmode for the coarse mode?!?

    Change log
    ----------
    2017-02-22: initially completed"""

    def __init__(self, data_quality='patchy',
                 diameter_cutoff='1um',
                 hygroscopicity_diameter=200,
                 wavelength = 550,
                 RH_dry=0,
                 RH_wet=85,
                 refractive_index ='aosacsm',
                 folder_out='/Users/htelg/data/ARM/myproducts/SGP/',
                 folder_tdmaaps='/Users/htelg/data/ARM/SGP/tdmaaps/',
                 folder_acsm='/Users/htelg/data/ARM/SGP/acsm/',
                 folder_tdmahyg='/Users/htelg/data/ARM/SGP/tdmahyg/',
                 keep_data = False,
                 test = False):

        if hygroscopicity_diameter != 200:
            raise ValueError('hygroscopicity diameters other then 200 have not been implemented yet')
        self.data_quality = data_quality
        self.diameter_cutoff = diameter_cutoff
        self.hygroscopicity_diameter = hygroscopicity_diameter
        self.RH_dry = RH_dry
        self.RH_wet = RH_wet
        self.wavelength = wavelength
        self.refractive_index = refractive_index
        self.black_carbon = None
        self.folder_out = folder_out
        self.folder = folder_tdmaaps
        self.folder_acsm = folder_acsm
        self.folder_tdmahyg = folder_tdmahyg
        self.test = test
        self.test_file = 'sgptdmaapssizeC1.c1.20120201.002958.cdf'
        self.keep_data = keep_data
        # self.verbose = verbose
        if self.test:
            self.keep_data = True
            self.intres = {}

    def _calculate_one(self, tdmaapssize, tdmahyg, refractive_index, diameter_cutoff,
                       hygroscopicity_diameter, wavelength, RH_dry, RH_wet):

        if diameter_cutoff == '1um':
            dcoff = 1000
        elif diameter_cutoff == '10um':
            dcoff = 10000

        dist = tdmaapssize.size_distribution.zoom_diameter(end=dcoff)
        dist.parameters4reductions.refractive_index = refractive_index
        dist.parameters4reductions.wavelength = wavelength
        dist.parameters4reductions.growth_distribution = tdmahyg.hyg_distributions_d200nm
        fRH = dist.hygroscopicity._get_f_RH(RH_dry, RH_wet)
        return fRH

    def calculate_new(self):
        self._calculate_all(False)

    def calculate_all(self):
        self._calculate_all(True)

    def _calculate_all(self, overwrite, time_window=False, verbose=False):
        if overwrite:
            fname = self.folder_out + 'README.txt'
            readme = open(fname, 'w')
            readme.write(self.__doc__)
            readme.close()

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

        all_files = _os.listdir(self.folder_acsm)
        all_files_acsm = _np.array(all_files)
        all_files = _os.listdir(self.folder_tdmahyg)
        all_files_tdmahyg = _np.array(all_files)

        test_done = False
        for e, fname_tdmaapssize in enumerate(all_files_tdmaapssize):
            if verbose:
                print('processing:',  fname_tdmaapssize)
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

            fname_others_acsm = _get_other_filenames(fname_tdmaapssize, ['aosacsm'], all_files_acsm)
            if not fname_others_acsm:
                if self.refractive_index == 'aosacsm':
                    if verbose:
                        print('acsm file missing .... continue')
                    continue

            fname_others_tdmahyg = _get_other_filenames(fname_tdmaapssize, ['tdmahyg'], all_files_tdmahyg)
            if not fname_others_tdmahyg:
                if verbose:
                    print('tdmahyg file missing .... continue')
                continue

            ior = str(self.refractive_index).replace('.','o')
            name_addon = '{}_hyg{:d}_rh{:d}v{:d}_ior{}_{}'.format(self.diameter_cutoff,
                                                     self.hygroscopicity_diameter,
                                                     self.RH_wet,
                                                     self.RH_dry,
                                                     ior,
                                                     self.data_quality)

            splitname = _splitup_filename(fname_tdmaapssize)
            site = splitname['site']
            date = splitname['date']

            my_prod_name = self.folder_out + site + 'tdmaapstdmahyg2fofrh_' + bc_name_signitur + name_addon + '.' + date +  '.000000.cdf'
            if not overwrite:
                if _os.path.isfile(my_prod_name):
                    if verbose:
                        print('product %s already exists' % my_prod_name)
                    continue
            # verbose = False

            if self.refractive_index == 'aosacsm':
                fname_others = _get_other_filenames(fname_tdmaapssize, ['aosacsm'], all_files_acsm)
                if not fname_others:
                    continue
                else:
                    aosacsm = _atm_arm.read_cdf(self.folder_acsm + fname_others['aosacsm']['fname'], data_quality=self.data_quality, verbose=verbose)
                    refractive_index = aosacsm.refractive_index

            elif type(self.refractive_index).__name__ == 'float':
                refractive_index = self.refractive_index

            tdmaapssize = _atm_arm.read_cdf(self.folder + fname_tdmaapssize, data_quality=self.data_quality, verbose=verbose)

            fname_others = _get_other_filenames(fname_tdmaapssize, ['tdmahyg'], all_files_tdmahyg)
            tdmahyg = _atm_arm.read_cdf(self.folder_tdmahyg + fname_others['tdmahyg']['fname'], data_quality=self.data_quality, verbose=verbose)

            fofrh = self._calculate_one(tdmaapssize, tdmahyg, refractive_index, diameter_cutoff=self.diameter_cutoff, wavelength= self.wavelength,
                                        hygroscopicity_diameter=self.hygroscopicity_diameter, RH_dry=self.RH_dry, RH_wet=self.RH_wet)

            if self.keep_data:
                fofrh_list.append(fofrh)
            if not self.test:
                fofrh.save_netCDF(my_prod_name)

        print(my_prod_name)
        if self.keep_data:
            if len(fofrh_list) > 1:
                fofrh_cat = _timeseries.concat(fofrh_list)
                self.result = fofrh_cat.close_gaps(verbose=False)
            else:
                self.result = fofrh_list[0]


def sgptdmaapstdmahyg2fofrh_1um_hyg200_rh85v0_ioraosacsm_patchy(test = False):
    out = Product(data_quality='patchy',
                       diameter_cutoff='1um',
                       hygroscopicity_diameter=200,
                       wavelength=550,
                       RH_dry=0,
                       RH_wet=85,
                       refractive_index='aosacsm',
                       folder_out='/Volumes/HTelg_4TB_Backup/arm_data/my_products/fRH/tdmaaps_tdmahyg_acsm/1um_hyg200_rh85v0_ioraosacsm_patchy/',
                       folder_tdmaaps='/Users/htelg/data/ARM/SGP/tdmaaps/',
                       folder_acsm='/Users/htelg/data/ARM/SGP/acsm/',
                       folder_tdmahyg='/Users/htelg/data/ARM/SGP/tdmahyg/',
                       keep_data=False,
                       test=test)
    return out

def sgptdmaapstdmahyg2fofrh_10um_hyg200_rh85v0_ioraosacsm_patchy(test = False):
    out = Product(data_quality='patchy',
                    diameter_cutoff='10um',
                    hygroscopicity_diameter=200,
                    wavelength=550,
                    RH_dry=0,
                    RH_wet=85,
                    refractive_index='aosacsm',
                    folder_out='/Volumes/HTelg_4TB_Backup/arm_data/my_products/fRH/tdmaaps_tdmahyg_acsm/10um_hyg200_rh85v0_ioraosacsm_patchy/',
                    folder_tdmaaps='/Users/htelg/data/ARM/SGP/tdmaaps/',
                    folder_acsm='/Users/htelg/data/ARM/SGP/acsm/',
                    folder_tdmahyg='/Users/htelg/data/ARM/SGP/tdmahyg/',
                    keep_data=False,
                    test=test)
    return out

def sgptdmaapstdmahyg2fofrh_10um_hyg200_rh85v0_ior1o5_patchy(test = False):
    out =  Product(data_quality='patchy',
                    diameter_cutoff='10um',
                    hygroscopicity_diameter=200,
                    wavelength=550,
                    RH_dry=0,
                    RH_wet=85,
                    refractive_index=1.5,
                    folder_out='/Volumes/HTelg_4TB_Backup/arm_data/my_products/fRH/tdmaaps_tdmahyg_acsm/10um_hyg200_rh85v0_ior1o5_patchy/',
                    folder_tdmaaps='/Users/htelg/data/ARM/SGP/tdmaaps/',
                    folder_acsm='/Users/htelg/data/ARM/SGP/acsm/',
                    folder_tdmahyg='/Users/htelg/data/ARM/SGP/tdmahyg/',
                    keep_data=False,
                    test=test)
    return out

def sgptdmaapstdmahyg2fofrh_1um_hyg200_rh85v0_ior1o5_patchy(test = False):
    out =  Product(data_quality='patchy',
                        diameter_cutoff='1um',
                        hygroscopicity_diameter=200,
                        wavelength=550,
                        RH_dry=0,
                        RH_wet=85,
                        refractive_index=1.5,
                        folder_out='/Volumes/HTelg_4TB_Backup/arm_data/my_products/fRH/tdmaaps_tdmahyg_acsm/1um_hyg200_rh85v0_ior1o5_patchy/',
                        folder_tdmaaps='/Users/htelg/data/ARM/SGP/tdmaaps/',
                        folder_acsm='/Users/htelg/data/ARM/SGP/acsm/',
                        folder_tdmahyg='/Users/htelg/data/ARM/SGP/tdmahyg/',
                        keep_data=False,
                        test=test)
    return out