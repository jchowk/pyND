from numpy import *

def sum_components(Col, Err, return_output=True, print_output=True):
    """
    output = sum_components(ColArray, ErrArray, return_output=True, print_output=True)
    
    :param Col: Array of input column density.
    :param Err: Array of input errors.
    :return: Output column, output error
    
    	11/10/97 -- Created by jch.
    	04/09/17 -- Translated to python.
        
    """

    n_params = 2
    sz1 = size(Col)

    err_scale = log10(exp(1))
    sum_col = 0
    sum_err = 0
    
    for i in arange(sz1):
        sum_col += 10.** (Col[i])
        sum_err += (10.**Col[i] * Err[i]/err_scale)**2
    
    sum_err = (sum_err) ** (0.5)
    
    log_sumCol = round(log10(sum_col),3)
    #log_sumErr = log10(1 + (sum_err / sum_col))
    log_sumErr = round((sum_err / sum_col)*err_scale,3)
    
    ## log_sumCol = alog10(total(10.0D^Col))
    ## log_sumErr = sqrt( total((10.0D^Col*2.3*Err)^2) )/total(10D^Col)/2.3
    

    if print_output:
        print("\t log N = {0} +/- {1}".format(log_sumCol, log_sumErr))

    def _ret():  return (log_sumCol, log_sumErr)

    if return_output:
        return _ret()

