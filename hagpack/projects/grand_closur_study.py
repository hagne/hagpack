from atmPy.general import timeseries

def save_figur(fname, fold_name, ax, save = False):
    a = ax
    f = a.get_figure()
    if save:
        f.tight_layout()
        f.savefig(fold_name+ '/' + fname, dpi = 300, transparent= True)

def do_correlate(data, correlant, data_name, correlant_name, whats_correlated, data_column=False,
                         correlant_column=False, xmax=False, ymax=False, zero_intersect=False, text_pos=(0.1, 0.9)):
    # out = data.correlate_to(correlant)


out = timeseries.correlate(data, correlant, data_column=data_column, correlant_column=correlant_column)
out._x_label_correlation = data_name
out._y_label_correlation = correlant_name
out._y_label_orig_correlant = correlant_name
out._y_label_orig_data = data_name
a, a1, a2 = out.plot_pearsonANDoriginal_data(zero_intersect=zero_intersect, xlim=xmax,
                                             ylim=ymax, corr_kwargs={'text_pos': text_pos})

title = 'Correlation %s  %s/%s' % (whats_correlated, data_name, correlant_name)
fname = 'correlation_%s__%s_%s.png' % (whats_correlated, data_name, correlant_name)
fname = fname.replace('(', '_').replace(')', '_').replace('=', '_').replace(' ', '')

f = a.get_figure()
f.suptitle(title)
save_it(fname, ax=a)
print('saved to: %s' % fname)
return out, a, a1, a2