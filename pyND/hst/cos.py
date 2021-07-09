def exp_time_func():
        import imp
        import numpy as np
        from numpy.polynomial import polynomial as P

        import pyND.hst as hst

        from astropy.table import Table


        # data_dir = hst.__path__[0]+'/'
        # exp = ascii.read(data_dir+'cos_g130m_g160m_snr10.txt',
        #                 header_start=0, data_start=1,format='csv')


        data_dir = imp.find_module('pyND')[1]+'/hst/'
        exp = Table.read(data_dir+'cos_g130m_g160m_snr10.csv')

        fuvmag = np.array(exp['FUV'])
        tg130m = np.array(exp['tg130m'])
        tg160m = np.array(exp['tg160m'])
        ti = np.arange(15.,20,0.1)
        # g130m fit
        coef1 = P.polyfit(fuvmag,tg130m,5)
        coef1 = np.ndarray.tolist(coef1)
        coef1.reverse()
        func1 = np.poly1d(coef1)
        tout_g130m = func1(ti)
        # g160m fit
        coef2 = P.polyfit(fuvmag,tg160m,5)
        coef2 = np.ndarray.tolist(coef2)
        coef2.reverse()
        func2 = np.poly1d(coef2)
        tout_g160m = func2(ti)
        return ti, tout_g130m,tout_g160m,func1,func2


def calc_orbits(fuvmag,torbit=50,snr=10):
        import numpy as np

        ti, tout_g130m,tout_g160m,func1,func2 = exp_time_func()

        t1 = func1(fuvmag) *(snr/10.)**2
        t2 = func2(fuvmag) *(snr/10.)**2
        torbit = torbit

        no1 = t1/torbit/60.
        no2 = t2/torbit/60.
        # ; first orbit
        over = (6+10+5)/torbit/60.
        # ; second
        over1 = (4 + 2)/torbit/60.
        no1 = no1 + over1
        if no1.astype(int) > 1:
            no1 = np.around(no1 + 4/torbit/60. * no1.astype(int),1)
        no2 = np.around(no2 + over1* no2.astype(int) +1./torbit/60.,1)

        return no1,no2,np.around(t1,1),np.around(t2,1)
