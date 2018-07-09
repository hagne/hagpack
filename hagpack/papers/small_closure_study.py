import numpy as _np
from scipy.optimize import curve_fit as _curve_fit
from scipy import interpolate as _interpolate
import pandas as _pd
import matplotlib.pylab as _plt
from matplotlib.ticker import MaxNLocator as _MaxNLocator
import plt_tools
from atmPy.general import timeseries as _timeseries
import atmPy as _atmPy
from matplotlib.colors import ListedColormap as _ListedColormap



_colors = _plt.rcParams['axes.prop_cycle'].by_key()['color']

def exp_decay(x, k, c):
    y = _np.exp(- x / (k / _np.log(2))) - c
    return y


def dofit(data, xmax=200):
    df = data.copy()
    df = df.iloc[1:, :]
    df[df.index > xmax] = _np.nan
    df = df.dropna()
    y = df.values.transpose()[0]
    x = df.index.values
    oa, ob = _curve_fit(exp_decay, x, y, p0=[500, 0.01])

    x = _np.append(x[::-1], 0)[::-1]
    y = exp_decay(x, oa[0], oa[1])
    df = _pd.DataFrame(y, index=x)
    out = {}
    out['df'] = df
    out['thalf'] = oa[0]
    out['offset'] = oa[1]
    return out


def plot_auto_corr_fct(auto_corr_fct_neph, auto_corr_fct_tdma,
                       auto_corr_fct_acsm=None,
                       auto_corr_fct_other=None,
                       plot_short_term = True,
                       plot_long_term = True,
                       fitres=True,
                       splabelpos = (0.06, 0.9),
                       cycle_start = 0,
                       fitrestxtpos=((0.35, 0.1), (0.75, 0.1)),
                       fitres_fontsize = 'medium',
                       txtpos=(0.4, 0.2),
                       txtoffset=(30, 50),
                       insert=True,
                       ax = None):
    colors = _plt.rcParams['axes.prop_cycle'].by_key()['color']
    if insert:
        f, a_acft = _plt.subplots()
        left, bottom, width, height = [0.6, 0.6, 0.2, 0.2]
        a_acfb = f.add_axes([left, bottom, width, height])
        a = [a_acfb, a_acft]
    elif ax:
        a = ax
        f = a.get_figure()
        a_acfb = a
        a_acft = None
        a = (a_acfb, a_acft)
    else:
        if plot_short_term and plot_long_term:
            f, a = _plt.subplots(2, )
            a_acfb = a[0]
            a_acft = a[1]
            f.set_figheight(f.get_figheight() * 1.5)
        elif plot_short_term:
            f, a = _plt.subplots(1, )
            a_acfb = a
            a_acft = None
            a = (a_acfb, a_acft)
        else:
            raise ValueError('plot_long_term only not implemented yet')
    # f.set_figwidth(width_tc)

    # a_empty = a[2]
    # a = a[:-1]
    ##################################################################################################
    auto_corr_fct_neph.plot(ax=a_acfb, marker='.')
    g = a_acfb.get_lines()[-1]
    g.set_linestyle('')
    g.set_marker('o')
    g.set_markersize(4)
    #     g.set_markeredgewidth(None)
    g.set_label('TNEPH')

    showfit = 250
    if showfit:
        fres = dofit(auto_corr_fct_neph, showfit)
        df = fres['df']
        df.plot(ax=a_acfb)
        g = a_acfb.get_lines()[-1]
        g.set_label(None)
        g.set_color(colors[0])
        g.set_linestyle('--')

    if fitres:
        txtl = []
        txtl.append('$t_{1/2} = %i$ min.' % (fres['thalf']))
        const = fres['offset']
        digs = int(_np.abs(_np.floor(_np.log10(const))))
        txtl.append('$c = %s$' % (round(const, digs)))
        onelag = auto_corr_fct_neph.iloc[1, 0]
        digs = int(_np.abs(_np.floor(_np.log10(1 - onelag))))
        txtl.append('$r_{1lag} = %s$' % (round(onelag, digs)))
        txt = '\n'.join(txtl)
        bbox = {}
        bbox['fc'] = [1, 1, 1, 0.7]
        bbox['ec'] = colors[0]
        bbox['boxstyle'] = 'round'
        a_acfb.text(fitrestxtpos[0][0], fitrestxtpos[1][0], txt, transform=a_acfb.transAxes, bbox=bbox, fontsize=fitres_fontsize)



    ##
    auto_corr_fct_tdma.plot(ax=a_acfb, marker='>', linestyle='')
    g = a_acfb.get_lines()[-1]
    g.set_label('tdma')
    g.set_markersize(4)
    g.set_color(colors[1])

    if showfit:
        fres = dofit(auto_corr_fct_tdma, showfit)
        df = fres['df']
        df.plot(ax=a_acfb)
        g = a_acfb.get_lines()[-1]
        g.set_label(None)
        g.set_color(colors[1])
        g.set_linestyle('--')

    if fitres:
        txtl = []
        txtl.append('$t_{1/2} = %i$ min.' % (fres['thalf']))
        const = fres['offset']
        digs = int(_np.abs(_np.floor(_np.log10(const))))
        txtl.append('$c = %s$' % (round(const, digs)))
        onelag = auto_corr_fct_tdma.iloc[1, 0]
        digs = int(_np.abs(_np.floor(_np.log10(1 - onelag))))
        txtl.append('$r_{1lag} = %s$' % (round(onelag, digs)))
        txt = '\n'.join(txtl)
        bbox = {}
        bbox['fc'] = [1, 1, 1, 0.7]
        bbox['ec'] = colors[1]
        bbox['boxstyle'] = 'round'
        a_acfb.text(fitrestxtpos[0][1], fitrestxtpos[1][1], txt, transform=a_acfb.transAxes, bbox=bbox, fontsize=fitres_fontsize)

    ##
    if not type(auto_corr_fct_acsm) == type(None):
        auto_corr_fct_acsm.plot(ax=a_acfb)
        g = a_acfb.get_lines()[-1]
        g.set_linestyle('')
        g.set_marker('s')
        g.set_markersize(4)
        g.set_label('ACSM')
        g.set_color(colors[2])

        if showfit:
            fres = dofit(auto_corr_fct_acsm, showfit)
            df = fres['df']
            df.plot(ax=a_acfb)
            g = a_acfb.get_lines()[-1]
            g.set_label(None)
            g.set_color(colors[2])
            g.set_linestyle('--')
        if fitres:
            txtl = []
            txtl.append('$t_{1/2} = %i$ min.' % (fres['thalf']))
            const = fres['offset']
            digs = int(_np.abs(_np.floor(_np.log10(const))))
            txtl.append('$c = %s$' % (round(const, digs)))
            onelag = auto_corr_fct_tdma.iloc[1, 0]
            digs = int(_np.abs(_np.floor(_np.log10(1 - onelag))))
            txtl.append('$r_{1lag} = %s$' % (round(onelag, digs)))
            txt = '\n'.join(txtl)
            bbox = {}
            bbox['fc'] = [1, 1, 1, 0.7]
            bbox['ec'] = colors[2]
            bbox['boxstyle'] = 'round'
            a_acfb.text(fitrestxtpos[0][2], fitrestxtpos[1][2], txt, transform=a_acfb.transAxes, bbox=bbox)
    ##
    if 1:
        if not type(auto_corr_fct_other) == type(None):
            auto_corr_fct_other.plot(ax=a_acfb)
            g = a_acfb.get_lines()[-1]
            g.set_label('Scatt. coeff.')
            g.set_color('0.4')
            g.set_linestyle(':')
            #             g.set_marker('.')

    # a_acfb.set_xlim(right = 210)
    # a_acfb.set_ylim(bottom=0.83)
    a_acfb.set_xlabel('Time lag (min.)')
    a_acfb.set_ylabel('Correlation coefficient r')
    a_acfb.legend(loc=1, fontsize='small')

    ##########################################################################################
    ##########################################################################################
    if plot_long_term:
        auto_corr_fct_neph = auto_corr_fct_neph.iloc[1:, :]
        auto_corr_fct_neph.plot(ax=a_acft)
        g = a_acft.get_lines()[-1]
        # g.set_linewidth(2)

        auto_corr_fct_tdma = auto_corr_fct_tdma.iloc[1:, :]
        auto_corr_fct_tdma.plot(ax=a_acft)
        #     g = a.get_lines()[-1]
        #     g.set_label('TDMAAPS')


        if not type(auto_corr_fct_acsm) == type(None):
            auto_corr_fct_acsm = auto_corr_fct_acsm.iloc[1:, :]
            auto_corr_fct_acsm.plot(ax=a_acft)
        # g = a.get_lines()[-1]
        #     g.set_label('ACSM')


        if not type(auto_corr_fct_other) == type(None):
            auto_corr_fct_other = auto_corr_fct_other.iloc[1:, :]
            auto_corr_fct_other.plot(ax=a_acft)
            g = a_acft.get_lines()[-1]
            g.set_label('Scatt. coeff.')
            g.set_color('0.4')
            g.set_linestyle('--')

        leg = a_acft.legend()
        leg.remove()
        a_acft.set_ylabel('Correlation coefficient r')
        a_acft.set_xlabel('Time lag (min.)')
        # a_acft.set_ylim(bottom=-0.03)

        a_acft.annotate('24 hours', xy=(24 * 60, txtpos[0]), xycoords='data',
                        xytext=txtoffset, textcoords='offset points',
                        arrowprops=dict(arrowstyle="->",
                                        connectionstyle="angle,angleA=0,angleB=90,rad=10"),
                        )
        a_acft.annotate('48 hours', xy=(48 * 60, txtpos[1]), xycoords='data',
                        xytext=txtoffset, textcoords='offset points',
                        arrowprops=dict(arrowstyle="->",
                                        connectionstyle="angle,angleA=0,angleB=90,rad=10"),
                        )

        # a_acft.xaxis.set_major_locator(MultipleLocator(1000))
        # a_acft.xaxis.set_minor_locator(MultipleLocator(500))
        ######
        abc = ['(a)', '(b)', '(c)']
        for e, at in enumerate(a):
            bbox = {}
            bbox['fc'] = [1, 1, 1, 0.7]
            bbox['ec'] = 'black'
            bbox['boxstyle'] = 'round'
            txtat = at.text(0.06, 0.9, abc[e], fontsize='large', transform=at.transAxes, bbox=bbox)
            bbp = txtat.get_bbox_patch()
            bbp.set_linewidth(_plt.rcParams['axes.linewidth'])

    if splabelpos:
        plt_tools.axes.text.set_subplots_labels(a, pos=splabelpos, cycle_start=cycle_start)
    return f, a


def plot_auto_corr_1lag_variability(roll1lag_neph, roll1lag_tdmaaps, contour_neph, contour_tdmaaps,
                                    roll1lag_acsm=None, roll1lag_other=None, contour_acsm=None, contour_other=None):
    f, aa = _plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0})
    a_1l, a_av = aa

    ################
    # 1lag
    roll1lag_neph.plot(ax=a_1l)
    g = a_1l.get_lines()[-1]
    g.set_label('TNeph')

    roll1lag_tdmaaps.plot(ax=a_1l)
    g = a_1l.get_lines()[-1]
    g.set_label('TDMAAPS')
    if not type(roll1lag_acsm) == type(None):
        roll1lag_acsm.plot(ax=a_1l)
        g = a_1l.get_lines()[-1]
        g.set_label('ACSM')
    if not type(roll1lag_other) == type(None):
        roll1lag_other.plot(ax=a_1l)
        g = a_1l.get_lines()[-1]
        g.set_color('0.4')
        g.set_linestyle('--')
        g.set_label('scatt. coeff.')
    # a_1l.set_ylim(top=1.1)
    a_1l.legend(  # loc = 4,
        #          frameon = True,
        #          fontsize = 'small',
    )
    a_1l.set_ylabel('$r_{1-lag}$')

    ######################
    # aerosol varibility
    contour_neph.plot(ax=a_av, label='neph')
    g = a_av.get_lines()[-1]
    g.set_label('neph')

    contour_tdmaaps.plot(ax=a_av, label='tdma')
    g = a_av.get_lines()[-1]
    g.set_label('tdma')

    if type(contour_acsm) != type(None):
        contour_acsm.plot(ax=a_av, label='acsm')
        g = a_av.get_lines()[-1]
        g.set_label('acsm')

    if type(contour_other) != type(None):
        contour_other.plot(ax=a_av, label='scatt. coeff.')
        g = a_av.get_lines()[-1]
        g.set_color("0.4")
        g.set_linestyle('--')
        g.set_label('scatt. coeff')

    # a_av.legend(loc = 1,
    #              frameon = True,
    #              fontsize = 'small',
    #             )
    a_av.set_ylabel('$\Delta t(r=0.8)$ (min.)')
    a_av.set_xlabel('')
    leg = a_av.legend()
    leg.remove()
    #######
    # axis decorations and limits
    #     a_1l.set_ylim(top = 1.02, bottom=0.76)
    for at in aa:
        at.set_xlim('2012', '2013')
        at.yaxis.set_major_locator(_MaxNLocator(6, prune='both'))

    # txts, bxs = set_subplots_labels(aa, pos = (0.9, 0.84), bbox_lw=0.5)
    plt_tools.axes.text.set_subplots_labels(aa, pos=(0.9, 0.84))
    return f, aa


def plot_data(a, data_neph, data_tdma,
              data_acsm=None,
              average=None,
              set_ylabel = None,
              labelpad=0.025,
              legend=None,
              width_sc = None):
    if average:
        #     avg = 10
        data_neph_avg = data_neph.average_time((average, 'm'))
        _timeseries.plot_wrapped(data_neph_avg, ax=a, lw=0.7)  # , label = r'Neph $t_{int} = 10\,$min')
        g = a[0].get_lines()[-1]
        col = g.get_color()
        _timeseries.plot_wrapped(data_neph, ax=a, lw=0.7, alpha=0.3, color=col, label="_nolegend_")
    else:
        _timeseries.plot_wrapped(data_neph,
                                 frequency='M',
                                 max_wraps=14,
                                 ax=a,
                                 zorder=1)

    a = _timeseries.plot_wrapped(data_tdma,
                                 frequency='M',
                                 max_wraps=14,
                                 ax=a,
                                 zorder=1)

    if type(data_acsm) != type(None):
        a = _timeseries.plot_wrapped(data_acsm,
                                     frequency='M',
                                     max_wraps=14,
                                     ax=a,
                                     zorder=0)

    f = a[0].get_figure()
    if width_sc:
        f.set_figheight(6 * 0.3 * width_sc)
    #     data_neph.plot_wrapped(ax = a)



    if set_ylabel:
        _atmPy.tools.plt_tools.set_shared_label(a, set_ylabel, axis='y', labelpad=labelpad)

    if legend:
        at = a[0]
        lines = at.get_lines()
        for e,l in enumerate(lines):
            l.set_label(legend[e])
        # gtdmaaps.set_label(legend[1])
        leg = at.legend(loc=1)
        leg.set_zorder(20)
    for at in a:
        #         at.set_ylim(0.08,0.25)
        at.yaxis.set_major_locator(_MaxNLocator(5, prune='both'))
        pass

    return f, a


def plot_the_ratio(theratio_neph_tdma,
                   ax=None,
                   labels = ['$\mathcal{Q}$'],
                   # instruments=('TNeph', 'TDMAAPS'),
                   txt_pos=[0.93, 0.93],
                   leg_orientation='vertical',
                   splabelpos = (0.9, 0.92),
                   cycle_start = 0
                   ):

    colors = _plt.rcParams['axes.prop_cycle'].by_key()['color']

    if ax:
        a = ax
        f = a.get_figure()
    else:
        f, a = _plt.subplots()

    ################
    # nephVtdma

    theratio_neph_tdma.plot(ax=a, label=labels[0])#r'$g_{{%s}}/g_{{%s}}$' % (instruments[0], instruments[1]))
    roll = theratio_neph_tdma.rolling((30, 'D'), min_good_ratio=0.3)
    mean = roll.median()
    std = roll.std()
    color = colors[1]

    mean.plot(ax=a, label=r'$\bar x_{30\,days}$')
    g = a.get_lines()[-1]
    g.set_color(color)

    a.fill_between(std.index, mean + std, mean - std, lw=0, linestyle=':', color=color, alpha=0.4, zorder=10,
                   label=r'$s_{30\,days}$')
    # a.set_ylim((0,1.5))

    props = dict(boxstyle='round',
                 facecolor=[1, 1, 1, 0.7],
                 linewidth=0.8
                 )
    # if
    if splabelpos:
        plt_tools.axes.text.set_subplots_labels(a, pos=splabelpos, cycle_start=cycle_start)
    txt = [r'$\bar x = {:0.2f}$'.format(theratio_neph_tdma.data.dropna().values.mean())]
    txt.append('$s = {:0.2f}$'.format(theratio_neph_tdma.data.dropna().values.std()))
    txt = '\n'.join(txt)

    a.text(txt_pos[0], txt_pos[1], txt, transform=a.transAxes, horizontalalignment='right', verticalalignment='top',
           bbox=props)

    ###############
    ## Legend
    loc = 2
    if leg_orientation == 'horizontal':
        leg = a.legend(loc=loc, ncol=3, columnspacing=1)
    elif leg_orientation == 'vertical':
        leg = a.legend(loc=loc)

    # leg.set_alpha(1)
    fr = leg.get_frame()
    fr.set_facecolor([0.9, 0.9, 0.9])
    fr.set_alpha(0.8)

    a.set_ylabel('Ratio')
    a.set_xlabel('')
    a.minorticks_off()

    return f, a, leg




def plot_correlation_main(corr,
                          corr_alt = None,
                          lims=(1, 4),
                          gridsize=40,
                          fr_pos=(0.75, 0.3),
                          fr_pos_alt = (0.75, 0.6),
                          label_addon = ',',
                          uncertainty = None,
                          uncertainty_text = '$1\sigma$',
                          sigma2 = False,
                          uncertainty_pos = [175, 150],
                          ax=None,
                          cycle_start=0,
                          splabels = True,
                          splabelpos=(0.9, 0.92),
                          show_residual = False,
                          residual_gridsize=40):
    if ax:
        at = ax
        f = at.get_figure()
        a = at
    else:
        if show_residual:
            f, a = _plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [3, 1],
                                                             'hspace': 0})
            at, ab = a
            f.set_figheight(f.get_figwidth() * 1.2)
        else:
            f, at = _plt.subplots()
            a = at

    #######
    ## corr
    cm = _plt.cm.gnuplot
    cm.set_under([0, 0, 0, 0])
    # if corr_alt:
    col = _colors[1]
    col_alt = _colors[5]
    # else:
    #     col = [0,0,0]
    if corr_alt:
        colt = col
    else:
        colt = _colors[2]#  [0,0,0,1]
    fit_res_kwargs = {'pos': fr_pos, 'show_params': ['r', 'm', 'c', 's'], 'bb_fc': [1, 1, 1, 0.8], 'bb_ec': colt}
    at, hb = corr.plot_regression(ax=at, reg_type='odr', gridsize=gridsize, xlim=lims, ylim=lims,
                              fit_res_kwargs=fit_res_kwargs, hexbin_kwargs={},
                              aspect='auto',
                              show_slope_1=False)
    hb.set_cmap(cm)
    hb.set_clim(vmin=0.001)
    hb.set_linewidths(0.2)
    hb.set_zorder(2)
    g = at.get_lines()[-1]
    # col = _plt.rcParams['axes.prop_cycle'].by_key()['color'][2]
    g.set_color(colt)
    g.set_zorder(11)
    #     g.set_linewidth(lw)
    g.set_label('fit$_{odr%s}$' % (label_addon))

    #######
    ## corr alt

    if corr_alt:
        alpha = 1
        cm_alt = _plt.cm.YlGn_r
        cm_alt = cm_alt(_np.arange(cm_alt.N))
        cm_alt[:, -1] = alpha
        cm_alt = _ListedColormap(cm_alt)
        cm_alt.set_under([1, 1, 1, 0])
        # col = _colors[0]
        fit_res_kwargs = {'pos': fr_pos_alt, 'show_params': ['r', 'm', 'c', 's'], 'bb_fc': [1, 1, 1, 0.8], 'bb_ec': col_alt}
        at, hb = corr_alt.plot_regression(ax=at, reg_type='odr', gridsize=gridsize, xlim=lims, ylim=lims,
                                          aspect='auto',
                                          fit_res_kwargs=fit_res_kwargs, hexbin_kwargs={}, show_slope_1=False)
        hb.set_cmap(cm_alt)
        hb.set_clim(vmin=0.001)
        hb.set_linewidths(0.2)
        hb.set_zorder(1)
        g = at.get_lines()[-1]
        #         col = plt.rcParams['axes.prop_cycle'].by_key()['color'][0]
        g.set_color(col_alt)
        g.set_zorder(10)
        #     g.set_linewidth(lw)
        #         g.set_label('fit$_{odr%s}$'%(label_addon))
        g.set_label('fit$_{odr}$')

    #######
    ## other stuff
    g, = at.plot([0, 300], [0, 300])
    #     g.set_linewidth(lw)
    col11 = 'r'#_colors[1]
    g.set_color(col11)
    g.set_linestyle('--')
    g.set_dashes((4, 4))
    g.set_label('1:1')
    at.set_xlim(lims)
    at.set_ylim(lims)
    leg = at.legend(loc=2, handlelength=2, fontsize='medium')

    #     splabelpos = (0.9,0.92)
    if splabels:
        plt_tools.axes.text.set_subplots_labels(at, pos=splabelpos, cycle_start=cycle_start)

    if uncertainty:
    # uncertainties
        x = _np.linspace(corr._data.min(), corr._data.max() * 2, 10)
        lin = lambda x, param, dm=1: param[0] + (param[1] * dm * x)
        corrout = corr.orthogonla_distance_regression['output']
        uncm, uncp = lin(x, corrout.beta, 1 - uncertainty), lin(x, corrout.beta, 1 + uncertainty)

    ## 1 sigma
        col = [0, 0, 0]
        alpha = 0.5
        ls = ''
        gm, = at.plot(x, uncm, color=col, alpha=alpha, ls=ls)
        gp, = at.plot(x, uncp, color=col, alpha=alpha, ls=ls)
        at.fill_between(x, uncm, uncp, zorder=0, color=col, alpha=0.15)

        tax = 'y'
        pos = uncertainty_pos
        # uncertainty_text = '$1\sigma$'
        if uncertainty_text:
            plt_tools.text.add_text_along_graph(gm, uncertainty_text, pos[1], axes=tax)
            plt_tools.text.add_text_along_graph(gp, uncertainty_text, pos[0], axes=tax)

    ## 2 sigma
        if sigma2:
            uncertain2s = uncertainty * 2
            uncm, uncp = lin(x, corrout.beta, 1 - uncertain2s), lin(x, corrout.beta, 1 + uncertain2s)
            alpha = 0.5
            ls = ''
            gm, = at.plot(x, uncm, color=col, alpha=alpha, ls=ls)
            gp, = at.plot(x, uncp, color=col, alpha=alpha, ls=ls)
            at.fill_between(x, uncm, uncp, zorder=0, color=col, alpha=0.1)

            tax = 'y'
            pos = uncertainty_pos[1]
            txt = '$2\sigma$'
            plt_tools.text.add_text_along_graph(gm, txt, pos, axes=tax)
            plt_tools.text.add_text_along_graph(gp, txt, pos, axes=tax)

    ## Residual
    if show_residual:
        ab, hb = corr.plot_residual(ax=ab, gridsize=residual_gridsize[0])
        ab.yaxis.set_major_locator(_MaxNLocator(4, prune='both'))
        hb.set_cmap(cm)
        hb.set_clim(vmin=0.001)
        hb.set_linewidths(0.2)
        hb.set_zorder(2)
        ab.set_ylabel('Residual')
        if corr_alt:
            ab, hb = corr_alt.plot_residual(ax=ab, gridsize=residual_gridsize[1])
            hb.set_cmap(cm_alt)
            hb.set_clim(vmin=0.001)
            hb.set_linewidths(0.2)
            hb.set_zorder(1)
    f = at.get_figure()
    return f, a


def plot_correlation_rolling(corroll_tdmaVneph, tdma_1lag_rollcorr, corroll_acsmVneph=None, acsm_1lag_rollcorr=None,
                             corroll_acsmVtdma=None, cycle_start=1, good_window = None):
    colors = _plt.rcParams['axes.prop_cycle'].by_key()['color']
    f, a = _plt.subplots()
    f.set_figheight(f.get_figheight() * 0.7)

    ###################
    ## tdma versus neph
    color = colors[0]
    corroll_tdmaVneph.plot(ax=a, zorder=1)
    g = a.get_lines()[-1]
    g.set_label('$r_{TNEPH, TDMAAPS}$')

    x = tdma_1lag_rollcorr.data.index.values
    y = tdma_1lag_rollcorr.data.iloc[:, 0].values

    col = plt_tools.colors.Color(color, model='hex')
    col.brightness = 0.9
    col.saturation = 0.2

    cole = plt_tools.colors.Color(color, model='hex')
    cole.brightness = 0.8
    cole.saturation = 0.3

    a.fill_between(x, y, 1,
                   #                    color = col.rgb,
                   #                    edgecolor = cole.rgb,
                   color=color,
                   alpha=0.3,
                   lw=0.5, label='$r_{1-lag,TDMAAPS}$')

    ###################
    ## acsm versus neph
    if type(corroll_acsmVneph) != type(None):
        color = colors[1]
        corroll_acsmVneph.plot(ax=a, zorder=1)
        g = a.get_lines()[-1]
        g.set_label('$r_{TNEPH, ACSM}$')
        g.set_color(color)

    if type(acsm_1lag_rollcorr) != type(None):
        color = colors[1]
        x = acsm_1lag_rollcorr.data.index.values
        y = acsm_1lag_rollcorr.data.iloc[:, 0].values

        col = plt_tools.colors.Color(color, model='hex')
        col.brightness = 0.9
        col.saturation = 0.2

        cole = plt_tools.colors.Color(color, model='hex')
        cole.brightness = 0.8
        cole.saturation = 0.3

        a.fill_between(x, y, 1,
                       #                        color = col.rgb,
                       #                        edgecolor = cole.rgb,
                       color=color,
                       alpha=0.3,
                       zorder=0,
                       lw=0.5, label='$r_{1-lag,ACSM}$')
    ###################
    ## acsm versus neph
    if type(corroll_acsmVtdma) != type(None):
        color = colors[2]
        corroll_acsmVtdma.plot(ax=a, zorder=1)
        g = a.get_lines()[-1]
        g.set_label('$r_{TDMAAPS, ACSM}$')
        g.set_color(color)

    ###################
    a.set_xlim(good_window)
    a.set_ylabel('Correlation coefficient $r$')
    a.set_xlabel('')
    a.legend()

    if cycle_start:
        splabelpos = (0.9, 0.88)
        txt, [box] = plt_tools.axes.text.set_subplots_labels(a, pos=splabelpos, cycle_start=cycle_start)
        box.set_facecolor([1, 1, 1, 0.5])

    return f, a


def plot_variability(sensibility, xlabel = 'Growth factor', text_prefix = 'gf = ', maxticks_r = 6, maxticks_m = 8, ax = None):
    lw = 1

    if type(ax) != type(None):
        a = ax
        f = a[0].get_figure()
    else:
        f, a = _plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0.05})

    a_r, a_m = a
    ## r
    x, y = sensibility.index, sensibility.pearson_r
    a_r.plot(x, y)

    fct = _interpolate.interp1d(x, y, kind='quadratic')
    xnew = _np.linspace(x.min(), x.max(), 500)
    ynew = fct(xnew)
    idxatmax = ynew.argmax()
    xatmax = xnew[idxatmax]

    a_r.plot(xnew, ynew)
    a_r.axvline(xatmax, color=_colors[1], lw=lw)

    bbox = {}
    txt = '{}{:0.2f}'.format(text_prefix, xatmax)
    t = a_r.text(xatmax, y.min() + (ynew.max() - y.min()) / 2, txt, bbox=bbox)
    t.set_verticalalignment('center')
    t.set_horizontalalignment('center')
    bb = t.get_bbox_patch()
    bb.set_boxstyle('round')
    bb.set_facecolor([1, 1, 1, 1])
    bb.set_edgecolor(_colors[1])
    bb.set_linewidth(lw)

    a_r.set_ylabel('r')
    a_r.yaxis.set_major_locator(_MaxNLocator(maxticks_r))

    ###############
    ## m
    x, y = sensibility.index, sensibility.odr_m
    a_m.plot(x, y)

    fct = _interpolate.interp1d(y, x, kind='quadratic')
    ynew = _np.linspace(y.max(), y.min(), 500)
    xnew = fct(ynew)
    xone = fct(1)

    a_m.plot(xnew, ynew)
    a_m.axvline(xone, color=_colors[1], lw=lw)
    a_m.axhline(1, color='0.6', ls=':')

    bbox = {}
    txt = '{}{:0.2f}'.format(text_prefix, xone)
    t = a_m.text(xone, 1 + (y.max() - 1) / 2, txt, bbox=bbox)
    t.set_verticalalignment('center')
    t.set_horizontalalignment('center')
    bb = t.get_bbox_patch()
    bb.set_boxstyle('round')
    bb.set_facecolor([1, 1, 1, 1])
    bb.set_edgecolor(_colors[1])
    bb.set_linewidth(lw)

    a_m.set_ylabel('m')
    a_m.yaxis.set_major_locator(_MaxNLocator(maxticks_m, prune='upper'))
    #######
    for at in a:
        g = at.get_lines()[0]
        g.set_linestyle('')
        g.set_marker('o')
        g.set_markersize(3)

        g = at.get_lines()[1]
        g.set_linestyle('--')
        g.set_color(_colors[0])

    at = a[-1]
    at.set_xlabel(xlabel)

    return f, a