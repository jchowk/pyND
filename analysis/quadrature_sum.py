def quadrature_sum(errors):
    import numpy as np

    variance = np.array(errors)** 2.
    finalerror = np.sqrt(variance.sum())
    
    return finalerror