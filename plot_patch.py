def plot_patch(patch, LatLims, LonLims, CM2, LonRng, colorscale, axis,
               cbarplot=True, cbar_title="Test", cbar_reverse=False,
               vn=0.10, vx=0.20, n=6, alpha=1.0, cbarvis=True,
               anchor='W'):
    """
    Purpose:
        To plot a map patch with appropriate latitude and longitude scales,
        and, optionally, a color bar

    Parameters
    ----------
    patch : 2D array
        Data to plot.
    LatLims : list or array
        [min, max] latitude limits.
    LonLims : list or array
        [min, max] longitude limits.
    CM2 : unused
        Kept for backward compatibility.
    LonRng : unused
        Kept for backward compatibility.
    colorscale : str or colormap
        Colormap to use for imshow.
    axis : matplotlib axis
        Axis to draw the patch on.
    cbarplot : bool, optional
        Whether to draw a colorbar (default True).
    cbar_title : str, optional
        Colorbar label (default "Test").
    cbar_reverse : bool, optional
        Reverse colorbar axis if True (default False).
    vn : float, optional
        Minimum value for color scale (default 0.10).
    vx : float, optional
        Maximum value for color scale (default 0.20).
    n : int, optional
        Number of ticks on colorbar (default 6).
    alpha : float, optional
        Alpha transparency (default 1.0).
    cbarvis : bool, optional
        Show colorbar if True (default True).
    anchor : str, optional
        Anchor main axis inside its grid cell ('W', 'E', 'C', etc.; default 'W').

    Returns
    -------
    patch, vn, vx, tx, cbar
    """

    import numpy as np
    import pylab as pl
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.ticker as mticker


    # --- anchor main axis ---
    axis.set_anchor(anchor)

    # --- clean NaNs and infs ---
    np.nan_to_num(patch, copy=False, nan=0.0, posinf=0.0, neginf=0.0)

    # --- colorbar ticks ---
    tx = np.linspace(vn, vx, n, endpoint=True)
    #tx = np.array(np.linspace(vn, vx, n, endpoint=True)).astype(int)
    #tx=[f"{x:.{4}g}" for x in txtemp]

    # --- plot patch ---
    show = axis.imshow(
        patch, colorscale, origin='upper', vmin=vn, vmax=vx,
        extent=[360-LonLims[0], 360-LonLims[1], 90-LatLims[1], 90-LatLims[0]],
        aspect="equal", alpha=alpha
    )

    # --- add colorbar tied to axis ---
    im_ratio = patch.shape[0] / patch.shape[1]
    if cbarplot:
        # use make_axes_locatable to tie colorbar to axis
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right", size=0.15*im_ratio, pad=0.15)
        cbar = axis.figure.colorbar(
            show, cax=cax, orientation='vertical', ticks=tx, cmap='gist_heat',
            format=mticker.FormatStrFormatter('%d')
        )
        formatter = mticker.PercentFormatter(xmax=100.0, decimals=0)
        cbar.ax.yaxis.set_major_formatter(formatter)
        #cbar.ax.set_yticklabels(np.around(tx, 3))
        cbar.ax.set_yticklabels(np.around(tx, 3))
        cbar.ax.tick_params(labelsize=7, color="k")
        cbar.ax.set_ylabel(cbar_title, size=7)
        cbar.ax.yaxis.set_label_coords(-1.0, 0.5)
        if cbar_reverse:
            cbar.ax.invert_yaxis()
        cbar.ax.set_visible(cbarvis)
    else:
        cbar = False

    return patch, vn, vx, tx, cbar