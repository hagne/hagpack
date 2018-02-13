from atmPy.general import timeseries
import numpy as np
import pandas as pd

def save_figur(fname, fold_name, ax):
    a = ax
    f = a.get_figure()
    if fold_name:
        f.tight_layout()
        f.savefig(fold_name+ '/' + fname, dpi = 300, transparent= True)




def do_correlate(fold_name, data, correlant, data_name, correlant_name, whats_correlated,
                 data_column=False,
                 correlant_column=False,
                 xmax=False,
                 ymax=False,
                 zero_intersect=False,
                 fit_res=(0.1, 0.9),
                 do_save = True):
    # out = data.correlate_to(correlant)
    def save_corr_res(output_folder_data):
        res = out.orthogonla_distance_regression['output']
        res_out = {'m': res.beta[0], 'c': res.beta[1]}
        res_out['std'] = np.sqrt(res.res_var)
        res_out['r'] = out.pearson_r[0]
        res_out

        df = pd.DataFrame(res_out, index=['corr_res'])

        fname = output_folder_data + 'corr_res.csv'
        df.to_csv(fname)
        return df

    out = timeseries.correlate(data, correlant, data_column=data_column, correlant_column=correlant_column)
    out._x_label_correlation = data_name
    out._y_label_correlation = correlant_name
    out._y_label_orig_correlant = correlant_name
    out._y_label_orig_data = data_name
    a, a1, a2 = out.plot_pearsonANDoriginal_data(zero_intersect=zero_intersect, xlim=xmax,
                                                 ylim=ymax, corr_kwargs={'fit_res': fit_res})

    title = 'Correlation %s  %s/%s' % (whats_correlated, data_name, correlant_name)
    fname = 'correlation_%s__%s_%s.png' % (whats_correlated, data_name, correlant_name)
    fname = fname.replace('(', '_').replace(')', '_').replace('=', '_').replace(' ', '')

    f = a.get_figure()
    f.suptitle(title)
    if fold_name:
        if do_save:
            save_figur(fname, fold_name, a)
            print('saved to: %s' % fname)
    out.save_corr_res = save_corr_res
    return out, a, a1, a2