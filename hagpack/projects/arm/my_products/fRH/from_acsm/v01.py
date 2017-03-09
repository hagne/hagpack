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
    """ Kappa is estimated from the chemical composition. This kappa is applied to sizedistributions ... here from
    TDMAAPS and growth thus fRH is calculated from those.

    ToDo
    ----
    use growthmode difference between 200 and 600 growthmode as the growthmode for the coarse mode?!?

    Change log
    ----------
    2017-02-28: initially completed
    2017-03-02: mode analysis implemented"""

    def __init__(self, data_quality='patchy',
                 diameter_cutoff='1um',
                 wavelength = 550,
                 RH_dry=0,
                 RH_wet=85,
                 refractive_index='aosacsm',
                 mode_analysis = False,
                 kappa_coarse = None,
                 folder_out= '/Volumes/HTelg_4TB_Backup/arm_data/my_products/fRH/',
                 folder_tdmaaps= '/Users/htelg/data/ARM/SGP/tdmaaps/',
                 folder_acsm= '/Users/htelg/data/ARM/SGP/acsm/',
                 keep_data = False,
                 test = False):

        self.data_quality = data_quality
        self.diameter_cutoff = diameter_cutoff
        self.RH_dry = RH_dry
        self.RH_wet = RH_wet
        self.wavelength = wavelength
        self.mode_analysis = mode_analysis
        self.kappa_coarse = kappa_coarse
        if mode_analysis:
            if type(refractive_index) != tuple:
                raise TypeError('if mode_analysis is true refractive_index must be a tuple, e.g. ("aosacsm", 1.5)')
            if not kappa_coarse:
                raise ValueError('if mode analysis is true kappa_coarse has to be set')
            self.refractive_index = refractive_index[0]
            self.refractive_index_coarse =  refractive_index[1]
        else:
            self.refractive_index = refractive_index
        self.folder_out = folder_out
        self.folder = folder_tdmaaps
        self.folder_acsm = folder_acsm
        self.test = test
        self.test_file = 'sgptdmaapssizeC1.c1.20120201.002958.cdf'
        self.keep_data = keep_data
        # self.verbose = verbose
        if self.test:
            self.keep_data = True
            self.intres = {}

    def _calculate_one(self, tdmaapssize, kappa, refractive_index, diameter_cutoff,
                       wavelength, RH_dry, RH_wet, mode_analysis = None):

        if diameter_cutoff == '1um':
            dcoff = 1000
        elif diameter_cutoff == '10um':
            dcoff = 10000

        dist = tdmaapssize.size_distribution.copy()

        # mode_analysis = True
        if mode_analysis:
            refractive_index_aiken_accu, refractive_index_coarse = refractive_index
            kappa_aiken_accu, kappa_coarse = kappa
            # mode selection and parameter setting
            dist = dist.convert2dVdlogDp()
            ## aiken accu
            dist_aiken_accu = dist.mode_analysis.size_dist_accu + dist.mode_analysis.size_dist_aiken
            dist_aiken_accu = dist_aiken_accu.zoom_diameter(end=dcoff)
            dist_aiken_accu.parameters4reductions.refractive_index = refractive_index_aiken_accu
            dist_aiken_accu.parameters4reductions.wavelength = wavelength
            dist_aiken_accu.parameters4reductions.kappa = kappa_aiken_accu

            ## coarse
            dist_coarse = dist.mode_analysis.size_dist_coarse
            dist_coarse = dist_coarse.zoom_diameter(end=dcoff)
            dist_coarse.parameters4reductions.refractive_index = refractive_index_coarse
            dist_coarse.parameters4reductions.wavelength = wavelength
            dist_coarse.parameters4reductions.kappa = kappa_coarse

            #scattcoeff
            ## aiken_accu
            ### dry
            dist_aiken_accu.parameters4reductions.RH = RH_dry
            scattcoeff_dry_aiken_accu = dist_aiken_accu.hygroscopicity.grown_size_distribution.optical_properties.scattering_coeff
            ### wet
            dist_aiken_accu.parameters4reductions.RH = RH_wet
            scattcoeff_wet_aiken_accu = dist_aiken_accu.hygroscopicity.grown_size_distribution.optical_properties.scattering_coeff

            ## coarse
            ### dry
            dist_coarse.parameters4reductions.RH = RH_dry
            scattcoeff_dry_coarse = dist_coarse.hygroscopicity.grown_size_distribution.optical_properties.scattering_coeff
            ### wet
            dist_coarse.parameters4reductions.RH = RH_wet
            scattcoeff_wet_coarse = dist_coarse.hygroscopicity.grown_size_distribution.optical_properties.scattering_coeff

            fRH = (scattcoeff_wet_coarse + scattcoeff_wet_aiken_accu) / (scattcoeff_dry_coarse + scattcoeff_dry_aiken_accu)


        else:
            dist = dist.zoom_diameter(end=dcoff)
            dist.parameters4reductions.refractive_index = refractive_index
            dist.parameters4reductions.wavelength = wavelength
            dist.parameters4reductions.kappa = kappa
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


        all_files = _os.listdir(self.folder)
        all_files = _np.array(all_files)
        all_files_tdmaapssize = all_files[_np.char.find(all_files, 'tdmaapssize') > -1]

        all_files = _os.listdir(self.folder_acsm)
        all_files_acsm = _np.array(all_files)


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
                # if self.refractive_index == 'aosacsm':
                if verbose:
                    print('acsm file missing .... continue')
                continue

            if self.mode_analysis:
                ior = str(self.refractive_index).replace('.','o')
                name_addon = '{}_rh{:d}v{:d}_ior{}_{}_kc{}_{}'.format(self.diameter_cutoff,
                                                         self.RH_wet,
                                                         self.RH_dry,
                                                         ior,
                                                         '{}j{}'.format(self.refractive_index_coarse.real, self.refractive_index_coarse.imag),
                                                         self.kappa_coarse,
                                                         self.data_quality)
                name_addon = name_addon.replace('.','o')
            else:
                ior = str(self.refractive_index).replace('.', 'o')
                name_addon = '{}_rh{:d}v{:d}_ior{}_{}'.format(self.diameter_cutoff,
                                                              self.RH_wet,
                                                              self.RH_dry,
                                                              ior,
                                                              self.data_quality)


            splitname = _splitup_filename(fname_tdmaapssize)
            site = splitname['site']
            date = splitname['date']

            my_prod_name = self.folder_out + site + 'acsm2fofrh_'  + name_addon + '.' + date +  '.000000.cdf'
            if not overwrite:
                if _os.path.isfile(my_prod_name):
                    if verbose:
                        print('product %s already exists' % my_prod_name)
                    continue
            # verbose = False

            fname_others = _get_other_filenames(fname_tdmaapssize, ['aosacsm'], all_files_acsm)
            if not fname_others:
                continue
            aosacsm = _atm_arm.read_cdf(self.folder_acsm + fname_others['aosacsm']['fname'], data_quality=self.data_quality, verbose=verbose)
            kappa = aosacsm.kappa
            if self.refractive_index == 'aosacsm':
                refractive_index = aosacsm.refractive_index

            elif type(self.refractive_index).__name__ == 'float':
                refractive_index = self.refractive_index

            if self.mode_analysis:
                refractive_index = (refractive_index, self.refractive_index_coarse)
                kappa = (kappa, self.kappa_coarse)

            tdmaapssize = _atm_arm.read_cdf(self.folder + fname_tdmaapssize, data_quality=self.data_quality, verbose=verbose)

            fofrh = self._calculate_one(tdmaapssize, kappa, refractive_index, diameter_cutoff=self.diameter_cutoff,
                                        wavelength= self.wavelength, RH_dry=self.RH_dry, RH_wet=self.RH_wet, mode_analysis=self.mode_analysis)
            # (self, tdmaapssize, kappa, refractive_index, diameter_cutoff, wavelength, RH_dry, RH_wet)
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

### Functions
###############################


def sgpacsm2fofrh_1um_rh85v0_ior1o5_patchy(test = False):
    out = Product(data_quality='patchy',
                      diameter_cutoff='1um',
                      wavelength=550,
                      RH_dry=0,
                      RH_wet=85,
                      refractive_index= 1.5,
                      folder_out='/Volumes/HTelg_4TB_Backup/arm_data/my_products/fRH/sgpacsm2fofrh/1um_rh85v0_ior1o5_patchy/',
                      folder_tdmaaps='/Users/htelg/data/ARM/SGP/tdmaaps/',
                      folder_acsm='/Users/htelg/data/ARM/SGP/acsm/',
                      keep_data=False,
                      test=test)
    return out

def sgpacsm2fofrh_10um_rh85v0_ior1o5_patchy(test = False):
    out = Product(data_quality='patchy',
                          diameter_cutoff='10um',
                          wavelength=550,
                          RH_dry=0,
                          RH_wet=85,
                          refractive_index=1.5,  # ('aosacsm', complex(1.63,0.01)),#"aosacsm"
                          mode_analysis=False,
                          kappa_coarse=0.0001,
                          folder_out='/Volumes/HTelg_4TB_Backup/arm_data/my_products/fRH/sgpacsm2fofrh/10um_rh85v0_ior1o5_patchy/',
                          folder_tdmaaps='/Users/htelg/data/ARM/SGP/tdmaaps/',
                          folder_acsm='/Users/htelg/data/ARM/SGP/acsm/',
                          keep_data=False,
                          test=test)
    return out

def sgpacsm2fofrh_1um_rh85v0_ioraosacsm_patchy(test = False):
    out = Product(data_quality='patchy',
                          diameter_cutoff='1um',
                          wavelength=550,
                          RH_dry=0,
                          RH_wet=85,
                          refractive_index= 'aosacsm',
                          folder_out='/Volumes/HTelg_4TB_Backup/arm_data/my_products/fRH/sgpacsm2fofrh/1um_rh85v0_ioraosacsm_patchy/',
                          folder_tdmaaps='/Users/htelg/data/ARM/SGP/tdmaaps/',
                          folder_acsm='/Users/htelg/data/ARM/SGP/acsm/',
                          keep_data=False,
                          test=test)
    return out

def sgpacsm2fofrh_10um_rh85v0_ioraosacsm_patchy(test = False):
    out = Product(data_quality='patchy',
                      diameter_cutoff='10um',
                      wavelength=550,
                      RH_dry=0,
                      RH_wet=85,
                      refractive_index= 'aosacsm',
                      folder_out='/Volumes/HTelg_4TB_Backup/arm_data/my_products/fRH/sgpacsm2fofrh/10um_rh85v0_ioraosacsm_patchy/',
                      folder_tdmaaps='/Users/htelg/data/ARM/SGP/tdmaaps/',
                      folder_acsm='/Users/htelg/data/ARM/SGP/acsm/',
                      keep_data=False,
                      test=test)
    return out



def sgpacsm2fofrh_1um_rh85v0_ior1o5_1o63j0o01_kc0o0001_patchy(test = False):
    out = Product(data_quality='patchy',
                          diameter_cutoff='1um',
                          wavelength=550,
                          RH_dry=0,
                          RH_wet=85,
                          refractive_index=(1.5, complex(1.63, 0.01)),  # "aosacsm"
                          mode_analysis=True,
                          kappa_coarse=0.0001,
                          folder_out='/Volumes/HTelg_4TB_Backup/arm_data/my_products/fRH/sgpacsm2fofrh/1um_rh85v0_ior1o5_1o63j0o01_kc0o0001_patchy/',
                          folder_tdmaaps='/Users/htelg/data/ARM/SGP/tdmaaps/',
                          folder_acsm='/Users/htelg/data/ARM/SGP/acsm/',
                          keep_data=False,
                          test=test)
    return out


def sgpacsm2fofrh_10um_rh85v0_ior1o5_1o63j0o01_kc0o0001_patchy(test = False):
    out = Product(data_quality='patchy',
                          diameter_cutoff='10um',
                          wavelength=550,
                          RH_dry=0,
                          RH_wet=85,
                          refractive_index=(1.5, complex(1.63, 0.01)),  # "aosacsm"
                          mode_analysis=True,
                          kappa_coarse=0.0001,
                          folder_out='/Volumes/HTelg_4TB_Backup/arm_data/my_products/fRH/sgpacsm2fofrh/10um_rh85v0_ior1o5_1o63j0o01_kc0o0001_patchy/',
                          folder_tdmaaps='/Users/htelg/data/ARM/SGP/tdmaaps/',
                          folder_acsm='/Users/htelg/data/ARM/SGP/acsm/',
                          keep_data=False,
                          test=test)
    return out


def sgpacsm2fofrh_1um_rh85v0_ioraosacsm_1o63j0o01_kc0o0001_patchy(test = False):
    out = Product(data_quality='patchy',
                          diameter_cutoff='1um',
                          wavelength=550,
                          RH_dry=0,
                          RH_wet=85,
                          refractive_index=('aosacsm', complex(1.63, 0.01)),  # "aosacsm"
                          mode_analysis=True,
                          kappa_coarse=0.0001,
                          folder_out='/Volumes/HTelg_4TB_Backup/arm_data/my_products/fRH/sgpacsm2fofrh/1um_rh85v0_ioraosacsm_1o63j0o01_kc0o0001_patchy/',
                          folder_tdmaaps='/Users/htelg/data/ARM/SGP/tdmaaps/',
                          folder_acsm='/Users/htelg/data/ARM/SGP/acsm/',
                          keep_data=False,
                          test=test)
    return out

def sgpacsm2fofrh_10um_rh85v0_ioraosacsm_1o63j0o01_kc0o0001_patchy(test= False):
    out = Product(data_quality='patchy',
                          diameter_cutoff='10um',
                          wavelength=550,
                          RH_dry=0,
                          RH_wet=85,
                          refractive_index=('aosacsm', complex(1.63, 0.01)),  # "aosacsm"
                          mode_analysis=True,
                          kappa_coarse=0.0001,
                          folder_out='/Volumes/HTelg_4TB_Backup/arm_data/my_products/fRH/sgpacsm2fofrh/10um_rh85v0_ioraosacsm_1o63j0o01_kc0o0001_patchy/',
                          folder_tdmaaps='/Users/htelg/data/ARM/SGP/tdmaaps/',
                          folder_acsm='/Users/htelg/data/ARM/SGP/acsm/',
                          keep_data=False,
                          test=test)
    return out