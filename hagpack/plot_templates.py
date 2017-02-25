from atmPy.tools import plt_tools as _plt_tools

def plot_rolling_time_laps_corr(rtlc_2D, rtlc_dtatmax =None, roll_corr = None, clim = (0.8, 1), title = None,
                                rtlc_dtatmax_label = 'dt at max corr.',
                                roll_corr_label = 'corr. coeff at dt = 0',
                                twin_y_label = 'corr. coeff at dt = 0',
                                cb_pad = 0.1,
                                save = False):
    # tl_max, tl_full = time_laps
    f ,a ,pc ,cb = rtlc_2D.plot(cb_kwargs={'pad' : cb_pad})
    f.set_figwidth(15)
    a.set_ylabel('dt (min)')

    cb.set_label('Correlation coeff. norm2max')
    cm = _plt_tools.get_colorMap_intensity_r()
    cm.set_bad('w')
    cm.set_under(color = 'w')
    pc.set_cmap(cm)
    pc.set_clim(clim  )  # 0.8,1)

    if rtlc_dtatmax:
        rtlc_dtatmax.plot(ax = a)
        g = a.get_lines()[-1]
        g.set_color('magenta')
        g.set_marker('.')
        g.set_label(rtlc_dtatmax_label)
        g.set_linewidth(2)
        leg = a.legend(fancybox = True, framealpha = 0.6)
        x = 0.05
        y = 0.05
        leg._loc = (x ,y)

    if roll_corr:
        a2 = a.twinx()
        roll_corr.plot(ax = a2)
        g = a2.get_lines()[-1]
        g.set_color('red')
        g.set_label(roll_corr_label)
        g.set_linewidth(2)
        a2.set_ylabel(twin_y_label)
        a2.spines['right'].set_color('red')
        a2.yaxis.label.set_color('red')
        a2.tick_params(axis='y', color = 'red', labelcolor = 'red')
        a.spines['right'].set_alpha(0)
        # leg2 = a2.legend(fancybox = True, framealpha = 0.6)
        # leg2._loc = (x+0.51,y)
        # lfram2 = leg2.get_frame()
        # lfram2.set_edgecolor('red')
    else:
        a2 = None

    #     a.set_xlim(right=735025)
    a.grid()

    if title:
        a.set_title(title)
    if save:
        f.savefig(save, dpi = 300, transparent = True)
    # f.savefig('time_laps_corr_all.png', dpi = 300, transparent = True)
    out = {'f':f,
           'a': a,
           'a2': a2,
           'pc': pc,
           'cb': cb,}
    return a, out

def plot_rolling_time_laps_corr_wrapped(rtlc_2D,
                                        rtlc_dtatmax=None,
                                        roll_corr=None,
                                        clim=(0.8, 1),
                                        title=None,
                                        rtlc_dtatmax_label='dt at max corr.',
                                        roll_corr_label='corr. coeff at dt = 0',
                                        twin_y_label='corr. coeff at dt = 0',
                                        save=False):


    # tl_max, tl_full = time_laps

    # rtlc_2D, rtlc_dtatmax =None,
    # roll_corr = None,
    # clim = (0.8, 1), title = None,
    #     rtlc_dtatmax_label = 'dt at max corr.'
    #     roll_corr_label = 'corr. coeff at dt = 0'
    #     twin_y_label = 'corr. coeff at dt = 0'
    # save = False):
    leg = 2

    frequency = 'Y'

    cm = _plt_tools.get_colorMap_intensity_r()
    cm.set_bad('w')
    cm.set_under(color='w')

    a = rtlc_2D.plot_wrapped(frequency=frequency, ylim=clim,
                             #                          cb_kwargs={'pad' :0.1},
                             pc_kwargs={'cmap': cm},
                             )
    f = a[0].get_figure()
    # f.set_figwidth(15)
    a[-1].set_ylabel('dt (min)')
    # cb = a.cb
    # cb.set_label('Correlation coeff. norm2max')

    if rtlc_dtatmax:
        rtlc_dtatmax.plot_wrapped(ax=a)
        for at in a:
            g = at.get_lines()[-1]
            g.set_color('magenta')
            g.set_marker('.')
            g.set_label(rtlc_dtatmax_label)
            g.set_linewidth(2)
    # leg = at.legend(fancybox = True, framealpha = 0.6)
    #         x = 0.05
    #         y = 0.05
    #         leg._loc = (x ,y)

    if roll_corr:
        # if 0:
        #     a2 = a.twinx()
        a = roll_corr.plot_wrapped(ax=a, twin_x=True)
        for e, at in enumerate(a.twins_x):
            at.set_ylim((0, 1))
            g = at.get_lines()[-1]
            g.set_color('red')
            #         g.set_label(roll_corr_label)
            g.set_linewidth(2)
            if e == 1:
                at.set_ylabel(twin_y_label)
            at.spines['right'].set_color('red')
            at.yaxis.label.set_color('red')
            at.tick_params(axis='y', color='red', labelcolor='red')
    # a.spines['right'].set_alpha(0)
    #     # leg2 = a2.legend(fancybox = True, framealpha = 0.6)
    #     # leg2._loc = (x+0.51,y)
    #     # lfram2 = leg2.get_frame()
    #     # lfram2.set_edgecolor('red')
    else:
        a2 = None

    # #     a.set_xlim(right=735025)
    for at in a:
        at.grid()

    # if title:
    #     a.set_title(title)
    # if save:
    #     f.savefig(save, dpi = 300, transparent = True)
    # # f.savefig('time_laps_corr_all.png', dpi = 300, transparent = True)
    # out = {'f':f,
    #        'a': a,
    #        'a2': a2,
    #        'pc': pc,
    #        'cb': cb,}
    #     return a, out
    # a[0].twinx()
    return f, a