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
    """Mark cosmic time atop a plot vs. redshift.

    *** WARNING: If using log axis scales, apply set_yscale('log') after running this code if you want to have minor tick marks.
    """

    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.units as u
    from astropy.cosmology import z_at_value

    from matplotlib.ticker import NullLocator, NullFormatter


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

    # Turn off the minor ticks
    ax_cosmo.xaxis.set_minor_locator(NullLocator())
    ax_cosmo.set_xlabel(top_label,fontsize=label_size)

    return ax_cosmo


def lookbacktimeaxis(ax_in, cosmo_in='',
                    ages_in=[2, 4, 6, 8, 9, 10, 11, 12, 13],
                    top_label = 'Lookback time (Gyr)',
                    label_size = 'medium',ticklabel_size='small',
                    **kwargs):
    """Mark cosmic time atop a plot vs. redshift.

    *** WARNING: If using log axis scales, apply set_yscale('log') after running this code if you want to have minor tick marks.
    """

    import numpy as np
    import astropy.units as u
    from astropy.cosmology import z_at_value

    from matplotlib.ticker import NullLocator


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


    ageticks = [z_at_value(cosmo_in.lookback_time, age) \
                    for age in ages_in]

    ax_cosmo = ax_in.twiny()
    ax_cosmo.set_xticks(ageticks)
    ax_cosmo.set_xticklabels(['{:g}'.format(age) for age in ages_in.value],
                                fontsize=ticklabel_size)
    ax_cosmo.set_xlim(ax_in.get_xlim())

    # Turn off the minor ticks
    ax_cosmo.xaxis.set_minor_locator(NullLocator())
    ax_cosmo.set_xlabel(top_label,fontsize=label_size)

    return ax_cosmo


def error_boxes(ax, xdata, ydata, xerror, yerror,
                    boxfacecolor='r',boxedgecolor='None', boxalpha=0.5,
                    **kwargs):
    """
    Following https://matplotlib.org/gallery/statistics/errorbars_and_boxes.html

    make_error_boxes(ax, xdata, ydata, xerror, yerror,
                        boxfacecolor='r',boxedgecolor='None', boxalpha=0.5,
                        **kwargs):

    ax == The axes object for the plot.

    """

    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle

    # Create list for all the error patches
    errorboxes = []

    # Loop over data points; create box from errors at each point
    for x, y, xe, ye in zip(xdata, ydata, xerror.T, yerror.T):
        rect = Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum())
        errorboxes.append(rect)

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=boxfacecolor, alpha=boxalpha,
                         edgecolor=boxedgecolor)

    # Add collection to axes
    ax.add_collection(pc)

    # Plot errorbars
    artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror,
                          **kwargs)

    return artists
