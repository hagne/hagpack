import pandas as pd




def _get_time(file_obj):
    bt = file_obj.variables['base_time']
    toff = file_obj.variables['time_offset']
    time = pd.to_datetime(0) + pd.to_timedelta(bt[:].flatten()[0], unit = 's') + pd.to_timedelta(toff[:], unit = 's')
    return time

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