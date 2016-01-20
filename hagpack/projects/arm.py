import pandas as pd
import numpy as np
from netCDF4 import Dataset
import os
from atmPy import sizedistribution
from atmPy.instruments.tools import diameter_binning
from atmPy import timeseries
from atmPy.aerosols import hygroscopic_growth as hg
from atmPy.tools import math_functions
from scipy.optimize import curve_fit



class ArmDict(dict):
    def __init__(self, plottable = [], plot_kwargs = {}, *args):
        super(ArmDict,self).__init__(self,*args)
        self.plottable = plottable
        self.plot_kwargs = plot_kwargs

    def plot(self, which = 'all', fig_size = None):
        if which == 'all':
            for item in self.plottable:
#                 f,a,b,c = self[item].plot(xaxis=0, yaxis = 2, sub_set=5)#*self.plot_kwargs)
#                 print(self.plot_kwargs)
                f,a,b,c = self[item].plot(**self.plot_kwargs)
                if fig_size:
                    f.set_size_inches((fig_size))
                return f,a,b,c


##################


def _read_tdmasize(file_obj):

    bt = file_obj.variables['base_time']
    toff = file_obj.variables['time_offset']
    index = pd.to_datetime(0) + pd.to_timedelta(bt[:].flatten()[0], unit = 's') + pd.to_timedelta(toff[:], unit = 's')

    sd = file_obj.variables['number_concentration']
    df = pd.DataFrame(sd[:])
    df.index = index

    d = file_obj.variables['diameter']
    bins, colnames = diameter_binning.bincenters2binsANDnames(d[:]*1000)

    dist = sizedistribution.SizeDist_TS(df,bins,'dNdlogDp')
    dist = dist.convert2dVdlogDp()

    out = ArmDict(plottable = ['size_distribution'])
    out['size_distribution'] = dist

    return out

def _concat_tdmasize(files):
    dist = files[0]['size_distribution']
    dist.data = pd.concat([i['size_distribution'].data for i in files])
    files[0]['size_distribution'] = dist
    return files[0]

####################

def _read_tdmaapssize(file_obj):

    index = _get_time(file_obj)

    sd = file_obj.variables['number_concentration_DMA_APS']
    df = pd.DataFrame(sd[:])
    df.index = index

    d = file_obj.variables['diameter']
    bins, colnames = diameter_binning.bincenters2binsANDnames(d[:]*1000)

    dist = sizedistribution.SizeDist_TS(df,bins,'dNdlogDp')
    dist = dist.convert2dVdlogDp()

    out = ArmDict(plottable = ['size_distribution'])
    out['size_distribution'] = dist
    return out


def _concat_tdmaapssize(files):
    dist = files[0]['size_distribution']
    dist.data = pd.concat([i['size_distribution'].data for i in files])
    files[0]['size_distribution'] = dist
    return files[0]

#####################

class Tdmahyg(ArmDict):
    def __init__(self,*args,**kwargs):
        print('loading  sdf tdmahy')
        super(Tdmahyg,self).__init__(*args,**kwargs)

    def calc_mean_growth_factor(self):
        """Calculates the mean growthfactor of the particular size bin."""
        def mean_linewise(gf_dist):
            growthfactors = self['hyg_distributions'].data.minor_axis.values
            meanl = ((gf_dist[~ gf_dist.mask] * np.log10(growthfactors[~ gf_dist.mask])).sum()/gf_dist[~gf_dist.mask].sum())
            stdl = np.sqrt((gf_dist[~ gf_dist.mask] * (np.log10(growthfactors[~ gf_dist.mask]) - meanl)**2).sum()/gf_dist[~gf_dist.mask].sum())
            return np.array([10**meanl,stdl])
        data = self['hyg_distributions'].data
        allmeans = timeseries.TimeSeries_3D(pd.Panel(items=data.items, major_axis=data.major_axis, minor_axis= ['mean', 'std_log']))
        for i,time in enumerate(data.values):
            for e,size in enumerate(time):
                allmeans.data.iloc[i,e] = mean_linewise(size)
        self['allmeans'] = allmeans
        return allmeans

    def calc_kappa_values(self):
        if not 'allmeans' in self.keys():
#             print('calc allmeans')
            allmeans = self.calc_mean_growth_factor()
        else:
#             print('don"t calc allmeans')
            allmeans = self['allmeans']

        RH = self['RH_interDMA']
        kappa_values = hg.kappa_simple(allmeans.data.values[:,:,0],RH, inverse = True)

        kappa_values = pd.DataFrame(kappa_values,columns=allmeans.data.major_axis, index = allmeans.data.items)
        out = timeseries.TimeSeries_2D(kappa_values)
        self['kappa_values'] = out
        self.plottable.append('kappa_values')
        return out

    def _fit_growth_factor(self, data):
        """Not recommended, probably not working"""
        def fit_linewise(gf_dist):
            fkt = math_functions.gauss
            amp = gf_dist.max()
            growthfactors = self['hyg_distributions'].data.minor_axis.values
            pos = growthfactors[gf_dist.argmax()]
            sigma = 0.4
            cof,varm = curve_fit(fkt,growthfactors[~gf_dist.mask], gf_dist[~gf_dist.mask], p0=[amp,pos,sigma])
            return np.array(cof)

        shape = list(data.shape)
        shape[-1] = 3
        fitres = np.zeros(shape)
        for i,time in enumerate(data):
            for e,size in enumerate(time):
                fitres[i,e] = fit_linewise(size)
        return fitres

def _read_tdmahyg(file_obj):
    "returns a dictionary, with panels in it"
    index = _get_time(file_obj)
    data = file_obj.variables['hyg_distributions'][:]
    growthfactors = file_obj.variables['growthfactors'][:]
    size_bins = file_obj.variables['size_bins'][:]* 1000
    RH_interDMA = pd.DataFrame(file_obj.variables['RH_interDMA'][:], index = index, columns=size_bins)
    RH_interDMA.columns.name = 'size_bin_center_nm'

    data = pd.Panel(data, items= index, major_axis = size_bins, minor_axis = growthfactors)
    data.items.name = 'Time'
    data.major_axis.name = 'size_bin_center_nm'
    data.minor_axis.name = 'growthfactors'

    out = Tdmahyg(plottable = ['hyg_distributions'], plot_kwargs =  dict(xaxis=0, yaxis = 2, sub_set=5, kwargs = dict(vmin = 0)))

#     data = timeseries.TimeSeries_3D(data)
    print('shape vor concat ', data.shape)
    data = timeseries.TimeSeries_3D(data)
#     data.RH_interDMA = RH_interDMA
    out['hyg_distributions'] = data
    out['RH_interDMA'] = RH_interDMA
#     out['growthfactors'] = growthfactors
#     out['size_bins'] = size_bins
    return out

def _concat_tdmahyg(files):
    out = Tdmahyg(plottable = ['hyg_distributions'], plot_kwargs =  dict(xaxis=0, yaxis = 2, sub_set=5, kwargs = dict(vmin = 0)))
#     data = timeseries.TimeSeries_3D(pd.concat([i['hyg_distributions'].data for i in files]))
    data = pd.concat([i['hyg_distributions'].data for i in files])
    data.iloc[:] = np.ma.masked_array(data.values, mask = data.values == -9999.0, fill_value = -9999.0)
    ts = timeseries.TimeSeries_3D(data)
#     data.RH_interDMA = RH_interDMA
    print('shape after concat ', ts.data.shape)
    out['hyg_distributions'] = ts
    out['RH_interDMA'] = pd.concat([i['RH_interDMA'] for i in files])
    return out


#################
arm_products = {'tdmasize':   {'read': _read_tdmasize,    'concat': _concat_tdmasize},
                'tdmaapssize':{'read': _read_tdmaapssize, 'concat': _concat_tdmaapssize},
                'tdmahyg':    {'read': _read_tdmahyg,     'concat': _concat_tdmahyg}
              }

def _get_time(file_obj):
    bt = file_obj.variables['base_time']
    toff = file_obj.variables['time_offset']
    time = pd.to_datetime(0) + pd.to_timedelta(bt[:].flatten()[0], unit = 's') + pd.to_timedelta(toff[:], unit = 's')
    return time

def read_arm_cdf(fname, concat = True, ignore_unknown = False, verbose = True):

    # list or single file
    if type(fname) == str:
        fname = [fname]
    products = {}

    #loop throuh files
    for f in fname:
        if verbose:
            print('\n', f)

        # error handling: test for netCDF file format
        if os.path.splitext(f)[-1] != '.cdf':
            txt = '\t %s is not a netCDF file ... skipping'%f
            if verbose:
                print(txt)
            continue

        # open file
        file_obj = Dataset(f)

        # unfortunatly Arm data is not very uniform, therefore the following mess
        if 'platform_id' in file_obj.ncattrs():
            product_id = file_obj.platform_id
        else:
            fnt = os.path.split(f)[-1].split('.')[0]
            foundit = False
            for prod in arm_products.keys():
                if prod in fnt:
                    product_id = prod
                    foundit = True
            if not foundit:
                txt = '\t has no ncattr named platform_id. Guess from file name failed ... skip'
                if verbose:
                    print(txt)
                continue

        # Error handling: if product_id not in products
        if product_id not in arm_products.keys():
            txt = 'Platform id %s is unknown.'%product_id
            if ignore_unknown:
                if verbose:
                    print(txt + '... skipping')
                continue
            else:
                raise KeyError(txt)

        if product_id not in products.keys():
            products[product_id] = []

        out = arm_products[product_id]['read'](file_obj)
        products[product_id].append(out)



    if len(fname) == 1:
        return out

    else:
        if concat:
            for pf in products.keys():
                products[pf] = arm_products[pf]['concat'](products[pf])
        return products




################