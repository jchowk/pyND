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


def asymmetric_err_distribution(p16, p50, p84, num_samples=1000, x0=None, verbose=False):
    """Create a skew normal distribution based on 16th, 50th, and 84th percentiles.
    
    :param p16 - 16th percentile value
    :param p50 - 50th percentile value (median)
    :param p84 - 84th percentile value
    :param num_samples(=1000) - Number of random samples to generate
    :param x0(=None) - Initial parameter guess [mu, sigma, alpha], default: [p50, (p84-p16)/2, 1]
    :param verbose(=False) - Print optimization results
    :return: tuple with (random_samples, parameters), where parameters is (mu, sigma, alpha)
    """
    import numpy as np
    from scipy.stats import skewnorm
    from scipy.optimize import minimize
    
    # Set default initial guess if not provided
    if x0 is None:
        x0 = np.array([p50, (p84-p16)/2, 1])
    
    # Constraint to ensure positive scale parameter
    cons = {'type': 'ineq', 'fun': lambda x: x[1]}
    
    # Define chi-square function to minimize
    def chisqfunc(params):
        mu, sigma, alpha = params
        chisq = (skewnorm.ppf(0.16, alpha, loc=mu, scale=sigma) - p16)**2 + \
                (skewnorm.ppf(0.50, alpha, loc=mu, scale=sigma) - p50)**2 + \
                (skewnorm.ppf(0.84, alpha, loc=mu, scale=sigma) - p84)**2
        return chisq
    
    # Perform optimization
    result = minimize(chisqfunc, x0, constraints=cons)
    
    # Extract optimized parameters
    mu, sigma, alpha = result.x
    
    if verbose:
        print(f"Optimization results: {result}")
        print(f"Optimized parameters: mu={mu:.4f}, sigma={sigma:.4f}, alpha={alpha:.4f}")
        print(f"Percentiles from fitted distribution:")
        print(f"P16: {skewnorm.ppf(0.16, alpha, loc=mu, scale=sigma):.4f} (target: {p16:.4f})")
        print(f"P50: {skewnorm.ppf(0.50, alpha, loc=mu, scale=sigma):.4f} (target: {p50:.4f})")
        print(f"P84: {skewnorm.ppf(0.84, alpha, loc=mu, scale=sigma):.4f} (target: {p84:.4f})")
    
    # Generate random samples
    samples = skewnorm.rvs(alpha, loc=mu, scale=sigma, size=num_samples)
    
    return samples, (mu, sigma, alpha)