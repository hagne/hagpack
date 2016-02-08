# -*- coding: utf-8 -*-
"""
Created on Tue May 26 09:07:22 2015

@author: htelg
"""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LightSource
import os
import pandas as pd
from osgeo import gdal
import pylab as plt
from mpl_toolkits.basemap import Basemap 
from scipy import ndimage 
from atmPy.tools import plt_tools
from atmPy import vertical_profile
from atmPy.aerosols import sampling_efficiency
from atmPy import sizedistribution

import colorsys


def read_gruvebadet_SMPS(fname):
    df = pd.read_csv(fname, sep='\t', skiprows=10, skip_blank_lines=True, error_bad_lines = False, encoding='iso8859_5')

    data = df.loc[:,[' 10.4', ' 11.1',
           ' 12.0', ' 12.9', ' 13.8', ' 14.9', ' 16.0', ' 17.2', ' 18.4', ' 19.8',
           ' 21.3', ' 22.9', ' 24.6', ' 26.4', ' 28.4', ' 30.5', ' 32.8', ' 35.2',
           ' 37.9', ' 40.7', ' 43.7', ' 47.0', ' 50.5', ' 54.2', ' 58.3', ' 62.6',
           ' 67.3', ' 72.3', ' 77.7', ' 83.5', ' 89.8', ' 96.5', '103.7', '111.4',
           '119.7', '128.6', '138.2', '148.6', '159.6', '171.5', '184.3', '198.1',
           '212.9', '228.8', '245.8', '264.2', '283.9', '305.1', '327.8', '352.3',
           '378.6', '406.8', '437.1', '469.8']]

    hk = df.drop(data.columns,axis=1)

    hk.index = pd.to_datetime(hk.Date + ' ' + hk['Start Time'], format = '%m/%d/%y %H:%M:%S')
    hk.index.name = 'Datetime_UTC'
    hk = hk.drop(['Date', 'Start Time'], axis = 1)

    data.index = hk.index

    cols = data.columns.copy()

    cols = cols.astype(float)

    avg_dist = np.array((np.log10(cols[1:]) - np.log10(cols[:-1]))).mean()

    bins = [cols[0] - avg_dist/2.]

    for e,c in enumerate(cols):
    #     delta = np.log10(cols[e]) - np.log10(bins[e])
        newbin = 10**(np.log10(cols[e]) + avg_dist/2.)
        bins.append(newbin)
    bins = np.array(bins)

    sd = sizedistribution.SizeDist_TS(data, bins, 'numberConcentration')
    return sd, hk

def read_gruvebadet_APS(fname):
    df = pd.read_csv(fname, sep='\t', skiprows=12, skip_blank_lines=True, error_bad_lines = False, encoding='iso8859_5')

    data = df.loc[:,['<0.523', '0.542', '0.583', '0.626', '0.673', '0.723', '0.777', '0.835', '0.898',
           '0.965', '1.037', '1.114', '1.197', '1.286', '1.382', '1.486', '1.596',
           '1.715', '1.843', '1.981', '2.129', '2.288', '2.458', '2.642', '2.839',
           '3.051', '3.278', '3.523', '3.786', '4.068', '4.371', '4.698', '5.048',
           '5.425', '5.829', '6.264', '6.732', '7.234', '7.774', '8.354', '8.977',
           '9.647', '10.37', '11.14', '11.97', '12.86', '13.82', '14.86', '15.96',
           '17.15', '18.43', '19.81']]
    hk = df.drop(data.columns,axis=1)

    hk.index = pd.to_datetime(hk.Date + ' ' + hk['Start Time'], format = '%m/%d/%y %H:%M:%S')
    hk.index.name = 'Datetime_UTC'
    hk = hk.drop(['Date', 'Start Time'], axis = 1)

    data.index = hk.index



    cols = data.columns.copy()

    bins = [float(cols[0].strip('<'))]
    cols = cols[1:]
    cols = cols.astype(float)

    avg_dist = np.array((np.log10(cols[1:]) - np.log10(cols[:-1]))).mean()

    for e,c in enumerate(cols):
    #     delta = np.log10(cols[e]) - np.log10(bins[e])
        newbin = 10**(np.log10(cols[e]) + avg_dist/2.)
        bins.append(newbin)
    bins = np.array(bins)
    bins *= 1e3
    data = data.drop('<0.523', axis = 1)

    sd = sizedistribution.SizeDist_TS(data, bins, 'numberConcentration')
    return sd, hk

def manta_sample_efficiency(particle_diameters = np.logspace(np.log10(0.14), np.log10(2.5),100),
                            manta_speed = 30, # m/s
                            pressure = 67., #kPa
                            main_inlet_diameter = 4.65 * 1e-3,
                            pick_off_diameter = 2.15 * 1e-3,
                            pick_off_flow_rate = 3,
                            lfe_diameter = 0.7 * 1e-3,
                            verbose = False):

    """Returns the manta sample efficiency for the POPS instrument (up most inlet)"""



    main_inlet_bent = sampling_efficiency.loss_in_a_bent_section_of_circular_tubing(
                                                  pressure = pressure,         # kPa
                                                  particle_diameter = particle_diameters,  # µm
                                                  tube_air_velocity = manta_speed, # m/s
                                                  tube_diameter = main_inlet_diameter,   # m
                                                  angle_of_bend = 90,       # degrees
                                                  flow_type = 'auto',
                                                  verbose = False)

    t_pick_of = sampling_efficiency.loss_in_a_T_junction(particle_diameter=particle_diameters,
                                       particle_velocity=30, 
                                       pick_of_tube_diameter=pick_off_diameter, 
                                       verbose=False)

    laminar_flow_element = sampling_efficiency.loss_at_an_abrupt_contraction_in_circular_tubing(pressure=pressure,  # kPa
                                                         particle_diameter=particle_diameters,  # µm
                                                         tube_air_velocity=False,  # m/s
                                                         flow_rate_in_inlet=pick_off_flow_rate,  # cc/s
                                                         tube_diameter=pick_off_diameter,  # m
                                                         contraction_diameter=lfe_diameter,  # m
                                                         contraction_angle=90,  # degrees
                                                         verbose=False,
                                                         )

    bent_before_pops = sampling_efficiency.loss_in_a_bent_section_of_circular_tubing(
                                                  pressure = pressure,         # kPa
                                                  particle_diameter = particle_diameters,  # µm
                                                  tube_air_velocity = False, # m/s
                                                  tube_air_flow_rate = pick_off_flow_rate,
                                                  tube_diameter = pick_off_diameter,   # m
                                                  angle_of_bend = 90,       # degrees
                                                  flow_type = 'auto',
                                                  verbose = False)

    gravitational_loss = sampling_efficiency.gravitational_loss_in_circular_tube(pressure=101.3,  # kPa
                                            particle_diameter=particle_diameters,  # µm
                                            tube_diameter=pick_off_diameter,  # m
                                            tube_length=0.25,  # m
                                            incline_angle=0,  # degrees from horizontal (0-90)
                                            flow_rate=3,  # cc/s
                                            mean_flow_velocity=False,  # 0.1061    # m/s)
                                            flow_type='auto',
                                            verbose=False)


    loss_list = [main_inlet_bent, t_pick_of, laminar_flow_element, bent_before_pops, gravitational_loss]
    names  = ['all_losses', 'main_inlet_bent', 't_pick_of', 'laminar_flow_element', 'bent_before_pops', 'gravitational_loss']

    all_losses = 1
    for l in loss_list:
        all_losses *= l

    loss_list.insert(0,all_losses)



    df = pd.DataFrame(np.array(loss_list).transpose(), columns = names, index = particle_diameters*1e3)
    df.index.name = 'diameters_nm'
    # if 0:
    #     fname = '/Users/htelg/data/20150414_Svalbard/particle_loss_in_inlet/sampling_efficiency.csv'
    #     df.to_csv(fname)
    return df



def PMEL_extCoeff_TS_get_plot_save(dist, fname = False):
#     fname = False
    wavelength = [450., 525., 624.]
    n = 1.455
    # AOD_list_PMEL = []
    AOD_dict_PMEL = {}
    for w in wavelength:
        AOD = dist.calculate_optical_properties(wavelength=w, n = n)
    #     calculate_AOD(wavelength=w, n = n)
    #     opt= sizedistribution.OpticalProperties(AOD, dist.bins)
    #     AOD_list.append({'wavelength':w, 'opt_inst': AOD})
        AOD_dict_PMEL['%.1f'%w] = AOD



    folder = '/Users/htelg/data/20150414_Svalbard/forPMEL/'

    keys = list(AOD_dict_PMEL.keys())
    keys.sort()
    df = pd.DataFrame()
    f,a = plt.subplots()

    for e, key in enumerate(keys):
        opt_tmp = AOD_dict_PMEL[key]
        ext_coeff = opt_tmp['extCoeff_perrow']
        df[key] = ext_coeff.data['ext_coeff']

        data = ext_coeff.data
        a.plot(data.index,data.values*1e5, color = plt_tools.wavelength_to_rgb(key))
    #     ext_coeff.data.plot(ax = a, color = plt_tools.wavelength_to_rgb(key))
        g = a.get_lines()[-1]
        g.set_label(key)
    a.set_ylabel('Extinction coefficient (10$^{-5}$ m$^{-1}$)')
    a.set_xlabel('Time (UTC)')
    # for k in keys:
    f.autofmt_xdate()
    a.legend()
    df = df.sort_index()
    df.index.name = 'TimeUTC'
    if fname:
        df.to_csv(folder+fname + '_ts' + '.csv')  
        f.savefig(folder+fname + '_ts' + '.png', dpi = 300)
    return a,df

def PMEL_extCoeff_get_plot_save(dist_LS, fname = False, as_time_series = False):

    wavelength = [450., 525., 624.]
    n = 1.455
    # AOD_list_PMEL = []
    AOD_dict_PMEL = {}
    for w in wavelength:
        AOD = dist_LS.calculate_optical_properties(wavelength=w, n = n)
    #     calculate_AOD(wavelength=w, n = n)
    #     opt= sizedistribution.OpticalProperties(AOD, dist_LS.bins)
    #     AOD_list.append({'wavelength':w, 'opt_inst': AOD})
        AOD_dict_PMEL['%.1f'%w] = AOD

    folder = '/Users/htelg/data/20150414_Svalbard/forPMEL/'

    keys = list(AOD_dict_PMEL.keys())
    keys.sort()
    df = pd.DataFrame()
    f,a = plt.subplots()

    for e, key in enumerate(keys):
        opt_tmp = AOD_dict_PMEL[key]
        ext_coeff = opt_tmp.get_extinction_coeff_verticle_profile()
        df[key] = ext_coeff.data['ext. coeff.']
        ext_coeff.plot(ax = a, color = plt_tools.wavelength_to_rgb(key))
        g = a.get_lines()[-1]
        g.set_label(key)
    # for k in keys:
    a.legend()
    if fname:
        f.savefig(folder+fname + '_vp' + '.png', dpi = 300)
        df.to_csv(folder+fname + '_vp' + '.csv')
    
    vp = vertical_profile.VerticalProfile(df)
    
    if type(as_time_series) == bool:
        if not as_time_series:
            ts = None
    else:    
        ts = vp.convert2timeseries(as_time_series)        
        aa = ts.data.plot()
        aa.set_ylabel('Extinction coefficient (m$^{-1}$)')
        ff = aa.get_figure()    
        if fname:
            ff.savefig(folder+fname + '_ts' + '.png', dpi = 300)
            ts.data.to_csv(folder+fname + '_ts' + '.csv')        
            
    return a,aa,vp, ts


def plot_POPS_v_mSASP_OD_corr(sun_int_su, 
                              aods_corr, 
                              offset=[3.295,3.44,3.995,4.16], 
                              additional_axes = False,
                              rayleigh = False,
                              error = False,
                              color = False):
    """plots the OD from miniSASP versus that which we simulate from the POPS results.
    Data is corrected for airmass factor.
    Arguments
    ---------
    sun_int_su: sun_int_su instance
    aod_corr: output of miniSASP.simulate_from_size_dist_LS
    rayleigh: bool
        if rayleigh is included or not.
    color: bool, or matplotlib color
        if bool, the color of the msasp channel is used"""
    
    airmassfct = False
    if not rayleigh:
        rayleigh_corr = aods_corr
    else:
        rayleigh_corr = rayleigh
    a = sun_int_su.plot(offset=offset, airmassfct = airmassfct, move_max = False, additional_axes = additional_axes, rayleigh = rayleigh_corr)

    f = a[0].get_figure()
    f.set_figwidth(15)
    # f.set_figheight(10)
    if rayleigh:
        a[0].set_xlabel('Optical depth')
    else:
        a[0].set_xlabel('Arosol optical depth')
    colors_fill = []
    colors_l = []
    for e,aa in enumerate(a):
    #     aa.set_xlim((-0.06,0.39))
        if e > 3:
            break
        aa.set_xlim((-0.01,0.14))
        aa.set_ylim((0,3200))
        aa.grid()
        g = aa.get_lines()[-1]
        txt = g.get_label()
        aa.set_title(txt)
        rgb_l = plt_tools.wavelength_to_rgb(float(txt.strip(' nm'))) * 0.8
        hls = list(colorsys.rgb_to_hls(*rgb_l))
        hls_sum = np.array(hls).sum()
        # hls[1] *=2. # lightness ... makes it brighter
        hls[1] = 0.8 # lightness ... makes it brighter
        # hls[2] *=0.4 #saturation ... makes it more gray
        hls[2] = 0.5 #saturation ... makes it more gray
        if hls_sum == 0:
            hls[1] = 0.6
            hls[2] = 0.0
        rgb = colorsys.hls_to_rgb(*hls)
        colors_fill.append(rgb)
        colors_l.append(rgb_l)
        # print('rgb', rgb_l)
        g.set_label('miniSASP')
        aa.locator_params(axis = 'x', nbins = 6)
    # print(colors_l)
    # print(colors_fill)

    #     aa.vlines(0,0,3200)
    keys = list(aods_corr.keys())
    keys.sort()
#    ms1,ms2,ms3,ms4 = (a[0].get_lines()[-1],a[1].get_lines()[-1],a[2].get_lines()[-1],a[3].get_lines()[-1])
    for e,k in enumerate(keys):
        out = aods_corr[k]

        if not np.any(color):
            col_l = colors_l[e]

        if rayleigh:
            g, = a[e].plot(out['sum'].values, out['sum'].index.values)
            g.set_linewidth(2)
            g.set_color(col_l)
            g.set_label('POPS')
        x = out['aerosol'].values
        y = out['aerosol'].index.values
        if np.any(error):
            col_fill = colors_fill[e]
            a[e].fill_betweenx(y, x * (1 - (error[0] / 100)), x * (1 + (error[1]/100)), color = col_fill)

        g2, = a[e].plot(x,y )
        g2.set_linewidth(2)
        g2.set_color(col_l)
        if rayleigh:
            g2.set_linestyle('--')
            g2.set_label('AOD only')
        else:

            g2.set_label('POPS')



        leg = a[e].legend(numpoints = 1,handlelength=1,  prop={'size':17})
        if e == 0:
            leg.draw_frame(True)
            txt = leg.get_texts()[0]
        #     txt.set_fontsize(10)
        else:
            leg.set_visible(False)
    
    return a



def get_colorMap_tundra():
    """elevation map according to a tundra climate """
    colors = []
#     color.append(np.array([0.,0.,0.])/255.) #white for ice
    colors.append(np.array([ 95., 93.,94.])/300.)
    colors.append(np.array([88.,120.,97.])/255.)
    colors.append(np.array([39., 62., 44.])/255.)
    colors.append(np.array([77.,102.,70.])/255.)
    colors.append(np.array([126., 129., 110.])/255.)
    colors.append(np.array([ 95., 93.,94.])/255.)
    colors.append(np.array([1.,1.,1.])) #white for ice
    
    steps = np.linspace(0,1,len(colors))
#     print(len(colors))
#    print(steps)
    red = []
    green = []
    blue = []
    
    for e,c in enumerate(colors):
        red.append((steps[e],c[0],c[0])) 
        green.append((steps[e],c[1],c[1])) 
        blue.append((steps[e],c[2],c[2])) 
        
    cdict = {'red':  red,
             'green': green,
             'blue':  blue
            }
    
    hag_cmap  = LinearSegmentedColormap('svalbard',cdict)
    hag_cmap.set_bad(np.array([10.,15.,100.])/255)
    return hag_cmap


def plot_map(town_label_loc = (50,50), town_label_alpha = 0.3):
    dem_in = DEM_map('/Users/htelg/projecte/svalbard/data/USGS_Nyalgen/Bulk_Order_462679/test/', 
              allInFolder=True
             )
             
    latMin, latMax = (78.8,79.1)
    lonMin, lonMax = (11.4,13.0)
    dem_in.zoom(latMin = latMin, latMax = latMax, lonMin = lonMin, lonMax = lonMax)
    
    #### some preparations
    dem_in.DataFrame.values[np.isnan(dem_in.DataFrame.values)] = 0
    newFrame = dem_in.DataFrame.copy()
    newFrame = newFrame.sort()
    
    #### create water
    water = newFrame.copy()
    water[water > 1.1] = np.nan
    sigma = 0.2
    water = ndimage.gaussian_filter(water, (sigma, 4*sigma) ,mode='nearest')
    water_masked = np.ma.array(water, mask=np.isnan(water))
    
    #### create relief
    sigma = 4
    altitude = ndimage.gaussian_filter(newFrame.values,  (sigma, 4*sigma) ,mode='nearest')
    ls = LightSource(azdeg=360-70,altdeg=65)
    relief = ls.shade(altitude,plt.cm.gray)
    
        
    #### Plot
    lat_center = 78.943736 
    lon_center = 12.2
    width = 34000
    height = 22000
    f,a = plt.subplots()
    bmap = Basemap(projection='laea',
                   lat_0=lat_center,
                   lon_0=lon_center,
                   width=width,
                   height=height,
    #                resolution='h'
                   ax = a
                  )
    # bmap
    # Fill the globe with a blue color
    
    
    lon, lat = np.meshgrid(newFrame.columns.values, newFrame.index.values)
    x,y = bmap(lon, lat)
    
    bmap.pcolormesh(x,y,relief[:,:,0], cmap = plt.cm.gray, vmin = -.4)
    # img.set_alpha(0.1)
    # img.set_edgecolor('face')
    
    bmap.pcolormesh(x,y,water_masked, 
                    cmap = get_colorMap_water()
                   )
    # img.set_alpha(0.1)
    # img.set_edgecolor('face')
    
    parallels = [78.9,78.95,79.0]
    g = bmap.drawparallels(parallels, labels=[True,True,True])
    for l in g.keys():
        g[l][0][0].set_linewidth(1)
    
    
    meridians = [11.6,12.0,12.4, 12.8]
    g = bmap.drawmeridians(meridians, labels=[True,True, True, True])
    for l in g.keys():
        g[l][0][0].set_linewidth(1)
    
    lon, lat = 11.922222, 78.925 # Location Ny-Alesund
    # convert to map projection coords.
    # Note that lon,lat can be scalars, lists or numpy arrays.
    xpt,ypt = bmap(lon,lat)
    p = bmap.plot(xpt,ypt,'bo')
    color = plt_tools.color_cycle[2]
    
    p[0].set_markerfacecolor('None')
    p[0].set_markeredgecolor(color)
    p[0].set_markeredgewidth(2)
    p[0].set_markersize(15)
    
    bmap.drawmapscale(lon_center - 0.35 , lat_center-0.075 , lon_center, lat_center, 10, barstyle='fancy', fontsize = 20)
    # a = img.get_axes()
    
    alpha = 0.3
    a.annotate('Ny-Ålesund', xy=(xpt, ypt),  
    #                 xycoords='data',
                    xytext=town_label_loc, 
                    size = 22,
                    ha="left",
                    va = 'top',
                    textcoords='offset points', 
                    bbox=dict(boxstyle="round", fc=[1,1,1,alpha], ec=color),
                    arrowprops=dict(#arrowstyle = "fancy", 
                                    arrowstyle="wedge,tail_width=0.7",
                                    fc=color, 
                                    ec="none",
    #                                 patchB=el,
                                    alpha = town_label_alpha,
                                    relpos = (-0.05,0.5),
                                    connectionstyle="angle3,angleA=0,angleB=-90"),
                    )
    # a.text(xpt,ypt,'bla', color = 'b')
    # f = img.get_figure()
    f.set_size_inches((12,12))
    out = {}
    out['f'] = f
    out['a'] = a
    out['bmap'] = bmap
    out['data_xy'] = (x,y)
    out['data_land_relief'] = relief[:,:,0]
    out['data_land_altitude'] = altitude
    out['data_water'] = water_masked
    return out
    
    
def add2map(bmap,hk):
    x, y = bmap(hk.data.Lon.values, hk.data.Lat.values)
    bmap.plot(x, y,
          color= plt_tools.color_cycle[1])
    

def get_colorMap_water():
    """elevation map according to a tundra climate """
    colors = []
#     color.append(np.array([0.,0.,0.])/255.) #white for ice
#     blue = np.array([ 0., 0., 50])/255.
    
    blue = np.array([161., 190., 255.]) / 255.
    colors.append(blue)
    colors.append(blue)
#     colors.append(np.array([39., 62., 44.])/255.)
#     colors.append(np.array([77.,102.,70.])/255.)
#     colors.append(np.array([126., 129., 110.])/255.)
#     colors.append(np.array([ 95., 93.,94.])/255.)
#     colors.append(np.array([1.,1.,1.])) #white for ice
    
    steps = np.linspace(0,1,len(colors))
#     print(len(colors))
#    print(steps)
    red = []
    green = []
    blue = []
    
    for e,c in enumerate(colors):
        red.append((steps[e],c[0],c[0])) 
        green.append((steps[e],c[1],c[1])) 
        blue.append((steps[e],c[2],c[2])) 
        
    cdict = {'red':  red,
             'green': green,
             'blue':  blue
            }
    
    hag_cmap  = LinearSegmentedColormap('svalbard',cdict)
    hag_cmap.set_bad(np.array([ 0., 0.,0.,0]))
    return hag_cmap

class DEM_map():
    def __init__(self,filename, allInFolder = False):
        if allInFolder:
            self.openAllInFolder(filename)
        else:
            self.open_NED_img_file(filename)
            self.create_transformationData()
            self.create_DataFrame()
            self.set_TransformationExtent()
#         self.create_DataTabel()

    def openAllInFolder(self, foldername, which = 'dem'):
        
        os.chdir(foldername)
        inDir = os.listdir('.')
#        tiles = []
        output = []
        for f in inDir:
            dirNf = os.path.split(f)
            fclean = os.path.splitext(dirNf[1])
            if fclean[1] == '.tif':
                if fclean[0].split('_')[2] == which:
                    self.open_NED_img_file(f)
                    self.create_transformationData()
                    self.create_DataFrame()
                    out = {'DataFrame': self.DataFrame.copy()}
                    out['transformations'] = self.get_transformationData()
                    output.append(out)
                    
        # create the concat datFrame
        # create the correct extent

        
        self.transformations['extent'] = get_commenExtent(output)

        newIndex = np.array([])
        newCols = np.array([])
        for o in output:
            newIndex = np.concatenate((newIndex, o['DataFrame'].index))
            newCols = np.concatenate((newCols, o['DataFrame'].columns))
        
        newIndex = np.unique(newIndex)[::-1]
        newCols = np.unique(newCols)
        newDat = pd.DataFrame(index=newIndex,columns=newCols, dtype=float)
        newDat.values[:] = 0.
        
        for o in output:
            newDat.loc[o['DataFrame'].index,o['DataFrame'].columns] = o['DataFrame'].values
            
        self.DataFrame = newDat
        self.DataFrameOriginal = self.DataFrame.copy()
        self.gdalInstance = False #because it would only represent the file opened last
        return
    
    def create_DataFrame(self):
        elevationArray = self.gdalInstance.ReadAsArray().astype(float)
#         elevationArray[np.where(elevationArray<0)] = np.nan
        elevationArray[elevationArray == -9999.0] = np.nan
        trans = self.get_transformationData()
        columns = np.linspace(trans['LonZero'],trans['LonMax'], trans['LonNoOfPts'])
        index = np.linspace(trans['LatZero'],trans['LatMax'], trans['LatNoOfPts'])
        self.DataFrame = pd.DataFrame(elevationArray, index = index, columns=columns)
#         self.DataFrame = self.DataFrame.astype(float)
        self.DataFrameOriginal = self.DataFrame.copy()
        return
    
    def get_DataFrame(self):
        return self.DataFrame
    
    def get_LatitudeExtend(self):
        minLat = self.DataFrame.index.values.min()
        maxLat = self.DataFrame.index.values.max()
        return (minLat,maxLat)
    
    def get_LongitudeExtend(self):
        minLon = self.DataFrame.columns.values.min()
        maxLon = self.DataFrame.columns.values.max()
        return (minLon, maxLon)
        
    def create_transformationData(self):
        trans = self.gdalInstance.GetGeoTransform()
        XZero = trans[0]
        XStepWidth = trans[1]
        XNoOfPts = self.gdalInstance.RasterXSize
        XMax = XZero+((XNoOfPts-1) * XStepWidth)
        YZero = trans[3]
        YStepWidth = trans[5]
        YNoOfPts = self.gdalInstance.RasterYSize
        YMax = YZero + ((YNoOfPts-1) * YStepWidth)
#         lat = self.get_LatitudeExtend()
#         lon = self.get_LongitudeExtend()
        self.transformations = {'LonZero' :     XZero, 
                                'LonStepWidth': XStepWidth,
                                'LonNoOfPts':   XNoOfPts,
                                'LonMax':       XMax,
                                'LatZero':      YZero, 
                                'LatStepWidth': YStepWidth,
                                'LatNoOfPts':   YNoOfPts,
                                'LatMax':       YMax,
                                'extent_orig': (XZero,XMax,YMax,YZero),}
#                                 'extent':     (XZero,XMax,YMax,YZero)}
        return self.transformations
    
    def get_transformationData(self):
        return self.transformations
        
    def open_NED_img_file(self,filename):
        """takes the .img file created from NED so it can be plotted by imshow
        Returns
            elevation array
            extend tupel
        """
        self.gdalInstance = gdal.Open(filename)
#        drv = self.gdalInstance.GetDriver()
        return self.gdalInstance
        
#     def get_elevationPdTable(self):
    def get_gdalInstance(self):
        return self.gdalInstance
        
#     def get_elevationArrayAndExtent(self):
        
#         elevationArray[np.where(elevationArray<0)] = np.nan
# #         trans = gdalInstance.GetGeoTransform()
#         extent = (trans[0],trans[0]+(gdalInstance.RasterXSize * trans[1]), trans[3], trans[3] + (gdalInstance.RasterYSize * trans[5]))
#         return elevationArray,extent
    
    def zoom(self, latMin = None, latMax = None, lonMin = None, lonMax = None):
        self.DataFrame = self.DataFrameOriginal.copy()
        if latMin:
            self.DataFrame = self.DataFrame.truncate(after = latMin)
        if latMax:
            self.DataFrame = self.DataFrame.truncate(before= latMax)
        if lonMax:
            self.DataFrame = self.DataFrame.truncate(after= lonMax, axis = 1)
        if lonMin:
            self.DataFrame = self.DataFrame.truncate(before = lonMin, axis = 1)
        self.set_TransformationExtent()

        return
    
    def set_TransformationExtent(self):
        lat = self.get_LatitudeExtend()
        lon = self.get_LongitudeExtend()
        self.transformations['extent'] = (lon[0],lon[1],lat[0],lat[1])
    
    def plot_Elevation(self, waterdepth = False, axes = False, cmap = get_colorMap_tundra()):
        """ axes: pass axes to print on """
        if type(axes).__name__ == 'bool':
            if not axes:
                f,a = plt.subplots()
        else:
            a = axes
            f = None
            
        if waterdepth:
            water = self.DataFrame.values.copy()
            water[water > 0] = np.nan

            surface = self.DataFrame.values.copy()
            surface[surface < 0] = np.nan

    #         a.imshow(water, extent = self.transformations['extent'], cmap = plt.cm.terrain)# (df.columns.max(),df.columns.min(),df.index.min(),df.index.max()))
            a.imshow(surface, extent = self.transformations['extent'], cmap = plt.cm.terrain)# (df.columns.max(),df.columns.min(),df.index.min(),df.index.max()))
        else:
            a.imshow(self.DataFrame.values, extent = self.transformations['extent'], cmap = cmap)# (df.columns.max(),df.columns.min(),df.index.min(),df.index.max()))
        
        return f,a
    
    def plot_Slope(self):
        f,a = plt.subplots()
        ims = a.imshow(self.slopeMap.values, extent = self.transformations['extent'], cmap = plt.cm.terrain)# (df.columns.max(),df.columns.min(),df.index.min(),df.index.max()))
#         f.colorbar(ims)
        return f,a,ims
    
    def plot_contourf(self):
        pass
    
    def get_dLatInMeter(self):
        dLat = self.transformations['LatStepWidth']
        dLat_m = 4.0075e7/360.*dLat
        return dLat_m

    def get_dLonInMeter(self):
        dLon = self.transformations['LonStepWidth']
        Lat = self.DataFrame.index.values.mean()
        dLong_m =np.cos(np.deg2rad(Lat)) * 4.0075e7/360. * dLon
        return dLong_m
    
    def create_slopeMap(self):
        dLon_m = self.get_dLonInMeter()
        dLat_m = self.get_dLatInMeter()
        dDataFrame_dLat, dDataFrame_dLong = np.gradient(self.DataFrame.values, *[dLat_m,dLon_m])
        self.slopeMap = pd.DataFrame(np.rad2deg(np.arctan(np.sqrt(dDataFrame_dLat**2 + dDataFrame_dLong**2))), index = self.DataFrame.index.values, columns=self.DataFrame.columns.values)
        return 
    
    def get_slopeMap(self):
        return self.slopeMap
        
def get_commenExtent(gdalList):
    extentAll = np.zeros((len(gdalList),4))
    for e,i in enumerate(gdalList):
        t = i['transformations']
        extentAll[e] = np.array(t['extent_orig'])

    extentFull = np.array([extentAll[:,0].min(), extentAll[:,1].max(), extentAll[:,2].min(), extentAll[:,3].max()])
    return extentFull