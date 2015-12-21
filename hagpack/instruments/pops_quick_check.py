from atmPy.instruments.POPS import calibration,peaks,housekeeping
from atmPy.atmosphere_standards import standard_atmosphere
from atmPy import sizedistribution
import os
from IPython.display import Markdown, display

def printmd(string):
    display(Markdown(string))

def load_plot_hk(fname):
    hk = housekeeping.read_csv(fname)

    f,a = plt.subplots()
    hk.data.Barometric_pressure.plot(ax = a)
    a.legend(loc = 'best')

    f,a = plt.subplots()
    hk.data.POPS_Flow.plot(ax = a)
    a.set_ylim((0,0.15))
    a.legend(loc = 'best')

    f,a = plt.subplots()
    hk.data.Particle_rate_nops.plot(ax = a)
    a.legend(loc = 'best')

    f,a = plt.subplots()
    hk.data.L_LP_MON.plot(ax = a)
    a.legend(loc = 'best')

    f,a = plt.subplots()
    hk.data['12V_Mon'].plot(ax = a)
    a.legend(loc = 'best')
    return hk

def load_dist(fname, cal = False,
              flow = 107./60.,
              skip = 2.,
              from_scratch = False,
              dir_out = '/Users/htelg/data/20151203_China/output/',
              bins = np.logspace(np.log10(140), np.log10(2500), 30),
              ):
    p, f = os.path.split(fname)
    flo = dir_out + p.split('/')[-1]
    fn = flo + '/' + f
    fo = fn+'_dist_TS.csv'
#     print('')
#     print(fname)
#     print(flo)
#     print(fn)
#     print(fo)
    if from_scratch:
        if not cal:
            raise ValueError('cal has to be a real calibration instance')
        dist_TS = peaks.read_cal_process_peakFile(fname,cal, bins, normalize = (1 / (skip+1)) * flow)
        try:
            dist_TS.save_csv(fo)
        except FileNotFoundError:
            os.mkdir(flo)
            dist_TS.save_csv(fo)
    else:
        dist_TS = sizedistribution.read_csv(fo)
    out = dist_TS.fit_normal()
    avg = dist_TS.average_overAllTime()
    return dist_TS, avg

def plot_dist(dist_TS,avg):
    f,a1,pc,cb = dist_TS.plot(norm='log', fit_pos = False)
    a1.set_title('Size distribution time series')

    a = dist_TS.plot_particle_concentration()
    # a.set_ylim((0.1,1200))
#     a.set_yscale('log')
    a.set_title('Particle concentration')

    f,a2 = avg.plot()
    a2.set_xlim((140,2500))
    a2.set_yscale('log')
    a2.set_title('Average over all time')
    return a,a2

def do_all(fname,
           bins,
           cal = False,
           from_scratch = True,
           dir_out = '/Users/htelg/data/20151216_china_ping/',
           verbose = False,
           flow = 1, # in cc
           skip = 0,
#            prefix = ''
          ):
    p,f = os.path.split(fname)
    bn = f[::-1].split('_',1)[1][::-1]
    bnt = p.split('/')[-1] + '/' + bn[::-1].split('_',1)[1][::-1]
    printmd('# %s'%bnt)

#     if 1:
    def sub_do():
        hkn = p + '/' + bn + '_HK.csv'
        dn = p + '/' + bn + '_Peak.bin'
    #     printmd('## Housekeeping')
        load_plot_hk(hkn)
        dist_TS, avg  = load_dist(dn, bins = bins, cal = cal, from_scratch= from_scratch, dir_out = dir_out,
                                  flow = flow,
                                  skip = skip )
    #     printmd('## Sizedistributions')
        plot_dist(dist_TS,avg)
    if verbose:
        sub_do()
    else:
        try:
            sub_do()
        except:
            printmd('## Failed')