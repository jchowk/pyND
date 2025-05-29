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
