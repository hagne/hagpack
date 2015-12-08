import os

def convert_image_format(fname, from_format, to_format):
    """Converts between image formats using image_magic (convert)

    Parameters
    ----------
    fname: str or list
        file name or list of file names
    from_format: str
        format of images to be converted, e.g. '.png'
    to_format: str
        format of images to be converted, e.g. '.png'"""

    if isinstance(fname, (str,list)):
        if isinstance(fname, (str)):
            all_files = [fname]
        else:
            all_files = fname
    else:
        txt = 'fnama has to be of type str or list'
        raise TypeError(txt)
    for f in all_files:
        base,ending = os.path.splitext(f)
    #     print(ending)
        if ending.lower() == from_format.lower():
            newfn = base + to_format
            print('converting %s to %s'%(f,newfn))
            os.system('convert %s %s'%(f,newfn))
    print('done')