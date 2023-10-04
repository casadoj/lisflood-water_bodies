import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cf
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Union, Dict, List, Tuple
from utils import Decomposition
from metrics import KGEmod
        
        
def create_cmap(cmap: str, bounds: List, name: str = '', specify_color: Tuple = None):
    """Given the name of a colour map and the boundaries, it creates a discrete colour ramp for future plots
    
    Inputs:
    ------
    cmap:          string. Matplotlib's name of a colourmap. E.g. 'coolwarm', 'Blues'...
    bounds:        list. Values that define the limits of the discrete colour ramp
    name:          string. Optional. Name given to the colour ramp
    specify_color: tuple (position, color). It defines a specific color for a specific position in the colour scale. Position must be an integer, and color must be either a colour name or a tuple of 4 floats (red, gren, blue, transparency)
    
    Outputs:
    --------
    cmap:   List of colours
    norm:   List of boundaries
    """
    
    cmap = plt.get_cmap(cmap)
    cmaplist = [cmap(i) for i in range(cmap.N)]
    if specify_color is not None:
        cmaplist[specify_color[0]] = specify_color[1]
    cmap = mpl.colors.LinearSegmentedColormap.from_list(name, cmaplist, cmap.N)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    return cmap, norm



def plot_resops(storage: pd.Series = None, elevation: pd.Series = None, inflow: pd.Series = None, outflow: pd.Series = None, capacity: Union[List[float], float] = None,
                   save: Union[str, Path] = None, **kwargs):
    """It creates a plot with two graphs that shows the reservoir time series. The first graph is the storage-elevation curve of the reservoir (if both storage and elevation time series are available). The second graph is the time series of storage, inflow and outflow.
    
    Parameters:
    -----------
    storage: pd.Series
        Time series of reservoir storage (hm3)
    elevation: pd.Series
        Time series of reservoir level (masl)
    inflow: pd.Series
        Time series of reservoir inflow (m3/s)
    outflow: pd.Series
        Time series of reservoir outflow (m3/s)
    capacity: Union[List[float], float]
        Values of storage capacity that will be plotted as horizontal lines in the time series plot
    save: Union[str, Path]
        File name where the plot will be saved. By default is None and the plot is not saved.
    kwargs:
        figsize: tuple
            Size of the plot
        xlim: tuple
            Limits of the X axis (time) in the time series plot
        ylim: tuple
            Limites of the Y axis (storage) in both plots. They share the Y axis
        title: str
            If given, title of the plot
    """
    
    # Create the figure and define the grid specification
    fig = plt.figure(figsize=kwargs.get('figsize', (20, 6)))
    gs = gridspec.GridSpec(1, 3)

    # Create the first graph in the left-most third
    ax1 = plt.subplot(gs[0])
    if isinstance(elevation, pd.Series) and isinstance(storage, pd.Series):
        ax1.scatter(elevation, storage, s=1, c= df.index, cmap='Greys')
    ax1.set(xlabel='elevation (m)',
            ylabel='storage (hm3)')

    # Create the second graph in the remaining area
    ax2 = plt.subplot(gs[1:], sharey=ax1)
    if isinstance(storage, pd.Series):
        ax2.fill_between(storage.index, storage, color='lightgray', alpha=.5, label='storage')
    if capacity is not None:
        if isinstance(capacity, float):
            capacity = [capacity]
        for ls, value in zip(['-', '--', ':'], capacity):
            ax2.axhline(value, lw=.8, c='k', ls=ls)
    ax2.set(xlim=kwargs.get('xlim', (storage.index[0], storage.index[-1])),
            ylim=(0, None))

    ax2_ = ax2.twinx()
    if isinstance(inflow, pd.Series):
        ax2_.plot(inflow, lw='.8', c='steelblue', label='inflow')
    if isinstance(outflow, pd.Series):
        ax2_.plot(outflow, lw='.8', c='indianred', label='outflow')
    ax2_.set(ylim=(0, None),
             ylabel='flow (m3/s)');
    
    fig.legend(ncol=3, loc=8, frameon=False, bbox_to_anchor=[.5, .0, .3, .1])
    fig.text(.5, .925, kwargs.get('title', None), fontsize=15, horizontalalignment='center')
    
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        
        
        
def plot_decomposition(obs: Decomposition, sim: Decomposition, id: int, lims: List[float] = [.1, .67, .97], save: Union[str, Path] = None, **kwargs):
    """It creates a figure that compares the decomposition of the observed and simulated time series. The figure is composed of 4 plots: original time series, trend, seasonality and residuals. Each plot includes the performance in terms of modified KGE and its components.
    
    Parameters:
    -----------
    obs:       Decomposition
        The result of decomposing the observed timeseries with the function utils.decompose_timeseries
    sim:       Decomposition
        The result of decomposing the simulated timeseries with the function utils.decompose_timeseries
    id:        int
        Identification number of the station to be plotted
    lims:      List[float]
        Values of the conservative (clim), normal (nlim) and flood (flim) limits in the LISFLOOD parameterization. If not provided, it takes the default values in GloFAS
    save:      bool
        Whether to save of not the figure. By default is None and the figure is note saved
        
    kwargs:
    -------
    figsize:   Tuple (2,)
        Size of each of the individual plots
    lw:        float
        Line width
    title:     str
        If provided, the title of the figure
    """
    
    # kwargs
    figsize = kwargs.get('figsize', (3, 6))
    lw = kwargs.get('lw', 1)
    
    # setup
    bbox_props = dict(boxstyle='round, pad=0.05', facecolor='w', edgecolor='none', alpha=.666)
    ncols = 4
    fig, axes = plt.subplots(figsize=(figsize[0] * ncols, figsize[1]), ncols=ncols, tight_layout=True, sharey=True)#, sharex=True)
    
    components = [attr for attr in dir(obs) if not attr.startswith('_')]
    for ax, comp in zip(axes, ['original', 'trend', 'seasonal', 'residual']):
        # extract series
        obs_comp = getattr(obs, comp)[id]
        sim_comp = getattr(sim, comp)[id]
        
        # line plots
        ax.plot(obs_comp, obs_comp.index, lw=lw, label='obs')
        ax.plot(sim_comp, sim_comp.index, lw=lw, label='sim')
        
        # add performance as text
        performance = KGEmod(obs_comp, sim_comp)
        performance = ['-∞' if x < -10 else '∞' if x > 10 else np.round(x, 2) for x in performance]
        ax.text(.02, .99, "KGE' = {0}\nα = {1}\nβ = {2}\nρ = {3}".format(*performance),
                va='top', transform=ax.transAxes, bbox=bbox_props)
        
        # settings
        if comp in ['original', 'trend']:
            for lim in lims:
                ax.axvline(lim, ls=':', c='k', lw=.5)
            ax.set(xlim=(-.02, 1.02),
                   ylim=(obs_comp.index.max(), obs_comp.index.min()))
        else:
            ax.set(xlim=(-.51, .51))
        ax.set(xlabel='Snorm',
               title=comp);
        
    fig.legend(*ax.get_legend_handles_labels(), ncol=2, loc=8, frameon=False, bbox_to_anchor=[.4, -.025, .2, .1])
    if 'title' in kwargs:
        fig.text(.5, 1, kwargs['title'], fontsize=11, ha='center');

    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight');