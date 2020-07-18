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


def skyplot(longitude, latitude,
            origin=0,
            zscale=None,
            zminmax=None,
            zlabel='Use zlabel=...',
            zcmap='viridis',
            system='Galactic',
            title=None,
            projection='aitoff'):
    '''From Eduardo Martín Calleja
    [http://balbuceosastropy.blogspot.com/2013/09/the-mollweide-projection.html]

    - longitude
    - latitude
    - zscale = None: variable for scaling point colors
    - zminmax = None: Optional min/max list for color scaling

    - origin = 0: The center of the plot in the longitude coordinate. For system='RADec', assumed in hours.
    - title = None: Plot title.

    - system = 'Galactic':
    - system = 'RADec':     longitude = RA in decimal degrees, latitude = Declination in deg

    - projection = 'aitoff': Projection type: 'mollweide', 'aitoff', 'hammer', 'lambert'
    '''

    from astropy.visualization import astropy_mpl_style
    from astropy.coordinates import frame_transform_graph

    import astropy.coordinates as coord
    import astropy.units as u

    import matplotlib.patheffects as path_effects
    import numpy as np

    import matplotlib.pyplot as plt


    # Default point color
    zcolor = 'seagreen'

    if (system == "RADec") & (origin <= 24.):
        origin = np.int(origin) * 15
    else:
        origin = np.int(origin)        

    # Shift longitude values
    x = np.remainder(longitude + 360 - origin, 360)
    ind = x > 180
    x[ind] -= 360  # scale conversion to [-180, 180]
    x = -x  # reverse the scale: East to the left


    # Do the plots
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection=projection)

    # Add the gridlines
    ax.axhline(0., color='k', linestyle='--', linewidth=0.5, alpha=0.9, zorder=1)
    ax.grid(True, color='0.7', linestyle=':', linewidth=0.5, alpha=0.8, zorder=0)

    # Set the title
    ax.set_title(title, pad=20)
    ax.title.set_fontsize(16)

    # Check which coordinate system we're using.
    if system == "Galactic":
        ax.set_xlabel("Galactic Longitude", fontsize=14)
        ax.set_ylabel("Galactic Latitude", fontsize=14)

        tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
        tick_labels = np.remainder(tick_labels + 360 + origin, 360)

    elif system == "RADec":
        ax.set_xlabel("Right Ascension", fontsize=14)
        ax.set_ylabel("Declination", fontsize=14)

        tick_labels = np.array([10, 8, 6, 4, 2, 0, 20, 18, 16, 14, 12])
        num_tick_labels = np.remainder(tick_labels + 24 + np.int(origin/15), 24)

        tick_labels = []
        for j in np.arange(len(num_tick_labels)):
            new_tick = np.str(num_tick_labels[j]) + '$^h$'
            tick_labels.append(new_tick)

    # Print the tick labels
    ax.set_xticklabels(tick_labels, color='k', fontsize=14,
                           path_effects=[path_effects.Stroke(linewidth=2,
                                                             foreground='w'),
                                         path_effects.Normal()])


    # Do the scatter plots
    # Check zscale information:
    if zscale is None:
        # No z coloring
        sky = ax.scatter(np.radians(x), np.radians(latitude),
                     s=50, marker='o', alpha=0.85,
                     edgecolors='w', color=zcolor)
    else:
        if zminmax is None:
            # Autoscale the z coloring
            sky = ax.scatter(np.radians(x), np.radians(latitude),
                     c=zscale,
                     s=50, marker='o', alpha=0.85,
                     edgecolors='w', cmap=zcmap)
        else:
            sky = ax.scatter(np.radians(x), np.radians(latitude),
                     c=zscale, vmin=zminmax[0], vmax=zminmax[1],
                     s=50, marker='o', alpha=0.85,
                     edgecolors='w', cmap=zcmap)
        # Plot the colorbar
        cb = plt.colorbar(sky, shrink=0.6)
        cb.set_label(label=zlabel,fontsize=14)

    plt.tight_layout(pad=0.5, w_pad=0.9, h_pad=0.01)
