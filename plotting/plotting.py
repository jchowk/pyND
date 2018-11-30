def plotaxes(pltwindow=None,**kwargs):
    """Mark zero and one for normalized spectra."""
    import matplotlib.pyplot as plt

    if 'linestyle' not in kwargs.keys():
        kwargs['linestyle'] = '--'
    if 'color' not in kwargs.keys():
        kwargs['color'] = 'k'
    if 'linewidth' not in kwargs.keys():
        kwargs['linewidth'] = 1.
    if 'zorder' not in kwargs.keys():
        kwargs['zorder'] = 0

    if pltwindow == None:
        plt.axhline(0,**kwargs)
        plt.axhline(1,**kwargs)
    else:
        pltwindow.axhline(0,**kwargs)
        pltwindow.axhline(1,**kwargs)

def plotzero(pltwindow=None,**kwargs):
    """Mark the zero point of a plot."""
    import matplotlib.pyplot as plt

    if 'linestyle' not in kwargs.keys():
        kwargs['linestyle'] = '--'
    if 'color' not in kwargs.keys():
        kwargs['color'] = 'k'
    if 'linewidth' not in kwargs.keys():
        kwargs['linewidth'] = 1.
    if 'zorder' not in kwargs.keys():
        kwargs['zorder'] = 0

    if pltwindow == None:
        plt.axhline(0,**kwargs)

        # plt.axhline(0,color='k',linestyle='--',
        #             linewidth=linewidth,zorder=zorder,kwargs)
    else:
        pltwindow.axhline(0,**kwargs)

def cosmictimeaxis(ax_in, cosmo_in='',
                    ages_in=[13, 10, 8, 6, 4, 3, 2, 1],
                    top_label = 'Time since Big Bang (Gyr)',
                    label_size = 'medium',ticklabel_size='small',
                    **kwargs):
    """Mark cosmic time atop a plot vs. redshift."""

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u
    from astropy.cosmology import z_at_value

    # Test that the ages have units associated with them.
    try:
        ages_in = np.array(ages_in)
        uuu=ages_in.unit
    except:
        ages_in=np.array(ages_in)*u.Gyr

    # Test that we've defined a cosmology!!
    try:
         cosmo_in.H0
    except:
        # Default to 737 cosmology
        from astropy.cosmology import FlatLambdaCDM
        cosmo_in = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K,
                       Om0=0.3, Ob0=0.045)


    ageticks = [z_at_value(cosmo_in.age, age) for age in ages_in]

    ax_cosmo = ax_in.twiny()
    ax_cosmo.set_xticks(ageticks)
    ax_cosmo.set_xticklabels(['{:g}'.format(age) for age in ages_in.value],
                                fontsize=ticklabel_size)
    ax_cosmo.set_xlim(ax_in.get_xlim())
    ax_cosmo.minorticks_off()
    #
    # zmin, zmax = 0.0, 5.9
    # ax.set_xlim(zmin, zmax)
    # ax_cosmo.set_xlim(zmin, zmax)

    ax_cosmo.set_xlabel(top_label,fontsize=label_size)

    return ax_cosmo
