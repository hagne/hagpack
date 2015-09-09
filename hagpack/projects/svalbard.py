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
                              save_fig = '/Users/htelg/tmp/POPS_versus_miniSASP_OD.png'):
    """plots the OD from miniSASP versus that which we simulate from the POPS results.
    Data is corrected for airmass factor.
    Arguments
    ---------
    sun_int_su: sun_int_su instance
    aod_corr: output of miniSASP.simulate_from_size_dist_LS"""
    
    airmassfct = False
    a = sun_int_su.plot(offset=offset, airmassfct = airmassfct, move_max = False)
    f = a[0].get_figure()
    f.set_figwidth(15)
    # f.set_figheight(10)
    a[0].set_xlabel('Optical depth')
    for aa in a:
    #     aa.set_xlim((-0.06,0.39))
        aa.set_xlim((-0.01,0.14))
        aa.set_ylim((0,3200))
        aa.grid()
        g = aa.get_lines()[-1]
        txt = g.get_label()
        aa.set_title(txt)
        g.set_label('miniSASP')
        aa.locator_params(axis = 'x', nbins = 6)

    #     aa.vlines(0,0,3200)
    keys = list(aods_corr.keys())
    keys.sort()
#    ms1,ms2,ms3,ms4 = (a[0].get_lines()[-1],a[1].get_lines()[-1],a[2].get_lines()[-1],a[3].get_lines()[-1])
    for e,k in enumerate(keys):
        out = aods_corr[k]
        g1 = a[e].get_lines()[-1]
        g, = a[e].plot(out['sum'].values, out['sum'].index.values)
        g.set_linewidth(2)
        g.set_color('r')
        g.set_label('POPS')

        g2, = a[e].plot(out['aerosol'].values, out['aerosol'].index.values)
        g2.set_linewidth(2)
        g2.set_color('r')
        g2.set_linestyle('--')
        g2.set_label('AOD only')


        leg = a[e].legend(numpoints = 1,handlelength=1,  prop={'size':17})
        if e == 0:
            leg.draw_frame(True)
            txt = leg.get_texts()[0]
        #     txt.set_fontsize(10)
        else:
            leg.set_visible(False)

    if save_fig:
        f.savefig(save_fig, dpi = 300)
    
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
    relief = ndimage.gaussian_filter(newFrame.values,  (sigma, 4*sigma) ,mode='nearest')
    ls = LightSource(azdeg=360-70,altdeg=65)
    relief = ls.shade(relief,plt.cm.gray)
    
        
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
    bmap
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
    
#    alpha = 0.3
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
    return f,a,bmap
    
    
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