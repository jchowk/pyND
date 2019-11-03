from __future__ import print_function, absolute_import, division, unicode_literals

def coveringfactor(input_data, quantiles = []):
    """
    eps=solarabundance(input,error=False,photo=False,meteor=False)

    :param input: Input elements (by symbol or Z)
    :param error: If True, return the error in abundance [default: False]
    :param photo: If True, force photospheric abundances [default: False]
    :param meteor: If True, force meteoritic abundances [default: False]
    :return: abundance or abundance, error
    """

    import numpy as np
    from astropy.io import fits
    import os
    import scipy.stats.distributions as dist

    def _ret():
        if error:
            return bestabundance, besterr
        else:
            return bestabundance
    confidenceRangeOuter = 0.95
    confidenceRangeInner = 0.68

    intervals = [(1.-confidenceRangeOuter)/2.,
                (1.-confidenceRangeInner)/2.,
                0.5,
                1.-(1.-confidenceRangeInner)/2.,
                1.-(1.-confidenceRangeOuter)/2.]

    coveringFactors = dist.beta.ppf(intervals, num_hits+1, num_trials-num_hits+1)

    return _ret()
