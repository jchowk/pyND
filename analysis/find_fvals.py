from numpy import *
from linetools.lists.linelist import LineList

def find_fvals(input_line,log_lambdaf=False,wavelength=False):
    """
    output = find_fvals(input_line,log_lambdaf=False,wavelength=False)
    
    :param input_line: linetools-style input, either symbol ('CIV 1548') or wavelength (1548.1)  
    :param log_lambdaf: If True, return log lam*f instead of f
    :param wavelength:  If True, return also the wavelength.
    :return: 
    """
    # Load in the linetools line lists.
    line_list = LineList('ISM', verbose=False,closest=True)
    user_line = line_list[input_line]

    # Check that we've actually got a line:
    bad = False
    if not user_line:
        ## BAIL!!
        print('No line information.')
        bad = True

    # Determine f-value for transition:
    fval = user_line['f']
    # Determine wavelength for transition:
    wave_out = user_line['wrest'].value
    # Log lambda*f
    lf = user_line['log(w*f)']

    def _ret():
        if log_lambdaf:
            #_rv = [round(lf,3)]
            _rv = round(lf, 3)
        else:
            _rv = round(fval,4)

        if wavelength:
            _rv = _rv,wave_out
            #_rv.append(wave_out)
        #return tuple(_rv)
        return _rv

    return _ret()