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