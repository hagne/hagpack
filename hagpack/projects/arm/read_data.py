from netCDF4 import Dataset
import os
from hagpack.projects.arm import _tdmasize,_tdmaapssize,_tdmahyg

arm_products = {'tdmasize':   {'read': _tdmasize._parse_netCDF,    'concat': _tdmasize._concat_rules},
                'tdmaapssize':{'read': _tdmaapssize._parse_netCDF, 'concat': _tdmaapssize._concat_rules},
                'tdmahyg':    {'read': _tdmahyg._parse_netCDF,     'concat': _tdmahyg._concat_rules}
              }








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