from __future__ import print_function, absolute_import, division, unicode_literals

def coveringfactor(num_hits, num_trials, 
                   quantiles = [0.025, 0.16, 0.5, 0.84, 0.975],
                   confidenceRange = None):
    """
    Calculates the covering factor and its confidence intervals using Bayesian
    inference assuming binomial statistics. This follows the formulation by Cameron, E. 2011, PASA, 28, 128.

    Parameters
    ----------
    num_hits : int
        Number of successes (hits) observed.
    num_trials : int
        Total number of trials.
    quantiles : list of float, optional
        List of quantiles to compute for the covering factor posterior
        distribution. Default is [0.025, 0.16, 0.5, 0.84, 0.975], corresponding
        to 95% and 68% confidence intervals and the median.
    confidenceRange : float or list of float, optional
        If provided, specifies the confidence interval(s) (e.g., 0.68 for 68%,
        0.95 for 95%). The function will compute the median and the lower/upper
        bounds for each confidence interval instead of using the quantiles
        argument.

    Returns
    -------
    coveringFactors : list of float
        List of the covering factors at the requested quantiles or confidence
        intervals.    

    Notes
    -----
    This function uses the beta distribution as the posterior for the covering
    factor given binomial data, following a Bayesian approach. If
    confidenceRange is specified, the function computes the median and the
    lower/upper bounds for each confidence interval. Otherwise, it returns the
    covering factors at the specified quantiles.

    Examples
    --------
    1) Calculate the median and 68% confidence interval for 2 detections out of 50 samples:

    >>> print(coveringfactor(2, 50, quantiles=[0.16,0.5,0.84]))
    [0.02708991 0.05208753 0.08837482]

    Thus the covering factor is 0.052 (+0.036, -0.025) at 60% confidence.

    2) Calculate the 90% confidence upper limit for 0 detections in 48 samples (we should get <0.046):

    >>> print(coveringfactor(0,48,confidenceRange=[0.8]))
    [0.0021479  0.01404628 0.04590452]

    OR

    >>> print(coveringfactor(0, 48, quantiles=[0.9]))
    [0.04590452]

    """
    import scipy.stats.distributions as dist
    
    if confidenceRange != None:
        # confidenceRange can be a single float or a list of floats
        if not isinstance(confidenceRange, (list, tuple)):
            confidenceRange = [confidenceRange]
        quantiles = [0.5]  # always include the median
        for cr in confidenceRange:
            lower = (1. - cr) / 2.
            upper = 1. - lower
            quantiles.extend([lower, upper])
        quantiles = sorted(set(quantiles))

    coveringFactors = dist.beta.ppf(quantiles, 
                                    num_hits+1, num_trials-num_hits+1)

    return coveringFactors
