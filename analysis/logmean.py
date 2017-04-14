def logmean(log_data, log_err, return_straight=False):
    """Calculate the weighted average (and errors) of a series of values given in log-space.
    	log_data   = Array holding the log-space data.
    	log_err = Array holding the log-space errors.
    
    	11/01/97 -- Created: jch
    	04/10/17 -- Transferred from IDL to Python
    """

    import numpy as np

    def _ret():
        if return_straight == True:
            return (np.log10(best), (best_err/best)/np.log(10.))
        else:
            return (log_mean, log_mean_err)

    # Convert inputs to numpy arrays
    log_data = np.array(log_data)
    log_err = np.array(log_err)

    # Convert to linear data, errors
    linear_data = 10. ** log_data
    linear_err = linear_data*(log_err * np.log(10))

    sum_data = (linear_data/linear_err**2).sum()
    sum_err = (1./linear_err**2).sum()

    mean = sum_data / sum_err
    mean_err = (1 / sum_err) ** (0.5)
    
    log_mean = np.log10(mean)
    log_mean_err = (mean_err / mean) / np.log(10.)
    
    print("Weighted Mean:")
    print("{0:0.3g} +/- {1:0.4g}".format(mean,mean_err))
    print("{0:0.3f} +/- {1:0.3f}".format(log_mean,log_mean_err))

    best = linear_data.mean()
    best_err = np.sqrt((linear_err**2).sum())

    print("Straight Mean:")
    print("{0:0.4g} +/- {1:0.4g}".format(best, best_err))
    print("{0:0.3f} +/- {1:0.3f}".format(np.log10(best), (best_err/best)/np.log(10.)))

    return _ret()
