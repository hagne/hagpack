# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 09:11:55 2015

@author: htelg
"""



def plot_hk_up_down_verticle(hk_up, hk_down, what, xlim = False, ylim = False, 
                             ylabel = False, xlabel = False, title = False, txt = False):
    """plot up and down of hk
    Arguments
    ---------
    hk_up, hk_down: obvious
    what: key of hk
    """
    ax = hk_up.plot_versus_altitude([what])
    g = ax[0].get_lines()[-1]
    g.set_label('up')
    ax =hk_down.plot_versus_altitude([what],ax = ax)
    g = ax[0].get_lines()[-1]
    g.set_label('down')
    f = ax[0].get_figure()
    ax[0].legend()
    f.set_size_inches((12,8))
    a = ax[0]
    if xlim:
        a.set_xlim(xlim)
    if ylim:
        a.set_ylim(ylim)
    if ylabel:
        a.set_ylabel(ylabel)
    if xlabel:
        a.set_xlabel(xlabel)
    if title:
        a.set_title(title)
    if txt:
        print(txt)
    return a