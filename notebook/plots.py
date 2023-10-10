import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
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
        
        
        
# def plot_decomposition(obs: Decomposition, sim: Decomposition, id: int, lims: List[float] = [.1, .67, .97], save: Union[str, Path] = None, **kwargs):
#     """It creates a figure that compares the decomposition of the observed and simulated time series. The figure is composed of 4 plots: original time series, trend, seasonality and residuals. Each plot includes the performance in terms of modified KGE and its components.
    
#     Parameters:
#     -----------
#     obs:       Decomposition
#         The result of decomposing the observed timeseries with the function utils.decompose_timeseries
#     sim:       Decomposition
#         The result of decomposing the simulated timeseries with the function utils.decompose_timeseries
#     id:        int
#         Identification number of the station to be plotted
#     lims:      List[float]
#         Values of the conservative (clim), normal (nlim) and flood (flim) limits in the LISFLOOD parameterization. If not provided, it takes the default values in GloFAS
#     save:      bool
#         Whether to save of not the figure. By default is None and the figure is note saved
        
#     kwargs:
#     -------
#     figsize:   Tuple (2,)
#         Size of each of the individual plots
#     lw:        float
#         Line width
#     title:     str
#         If provided, the title of the figure
#     """
    
#     # kwargs
#     figsize = kwargs.get('figsize', (3, 6))
#     lw = kwargs.get('lw', 1)
    
#     # setup
#     bbox_props = dict(boxstyle='round, pad=0.05', facecolor='w', edgecolor='none', alpha=.666)
#     ncols = 4
#     fig, axes = plt.subplots(figsize=(figsize[0] * ncols, figsize[1]), ncols=ncols, tight_layout=True, sharey=True)#, sharex=True)
    
#     components = [attr for attr in dir(obs) if not attr.startswith('_')]
#     for ax, comp in zip(axes, ['original', 'trend', 'seasonal', 'residual']):
#         # extract series
#         obs_comp = getattr(obs, comp)[id]
#         sim_comp = getattr(sim, comp)[id]
        
#         # line plots
#         ax.plot(obs_comp, obs_comp.index, lw=lw, label='obs')
#         ax.plot(sim_comp, sim_comp.index, lw=lw, label='sim')
        
#         # add performance as text
#         performance = KGEmod(obs_comp, sim_comp)
#         performance = ['-∞' if x < -10 else '∞' if x > 10 else np.round(x, 2) for x in performance]
#         ax.text(.02, .99, "KGE' = {0}\nα = {1}\nβ = {2}\nρ = {3}".format(*performance),
#                 va='top', transform=ax.transAxes, bbox=bbox_props)
        
#         # settings
#         if comp in ['original', 'trend']:
#             for lim in lims:
#                 ax.axvline(lim, ls=':', c='k', lw=.5)
#             ax.set(xlim=(-.02, 1.02),
#                    ylim=(obs_comp.index.max(), obs_comp.index.min()))
#         else:
#             ax.set(xlim=(-.51, .51))
#         ax.set(xlabel='Snorm',
#                title=comp);
        
#     fig.legend(*ax.get_legend_handles_labels(), ncol=2, loc=8, frameon=False, bbox_to_anchor=[.4, -.025, .2, .1])
#     if 'title' in kwargs:
#         fig.text(.5, 1, kwargs['title'], fontsize=11, ha='center');

#     if save is not None:
#         plt.savefig(save, dpi=300, bbox_inches='tight');



def plot_decomposition(sim: pd.DataFrame, obs: pd.DataFrame = None, lims: List[float] = None, save: Union[str, Path] = None, **kwargs):
    """It creates a figure that compares the decomposition of the observed and simulated time series. The figure is composed of 4 plots: original time series, trend, seasonality and residuals. Each plot includes the performance in terms of modified KGE and its components.
    
    Parameters:
    -----------
    sim:       Decomposition
        The result of decomposing the simulated timeseries with the function utils.decompose_timeseries
    obs:       Decomposition
        The result of decomposing the observed timeseries with the function utils.decompose_timeseries
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
    xlim = kwargs.get('xlim', [(-.02, 1.02), (-.51, .51)])
    xlabel = kwargs.get('xlabel', 'Snorm')
    
    # setup
    bbox_props = dict(boxstyle='round, pad=0.05', facecolor='w', edgecolor='none', alpha=.666)
    ncols = 4
    fig, axes = plt.subplots(figsize=(figsize[0] * ncols, figsize[1]), ncols=ncols, tight_layout=True, sharey=True)#, sharex=True)
    
    components = [attr for attr in dir(obs) if not attr.startswith('_')]
    for ax, comp in zip(axes, ['original', 'trend', 'seasonal', 'residual']):
        if obs is not None:
            # observed time series
            obs_comp = obs[comp] #getattr(obs, comp)[id]
            ax.plot(obs_comp, obs_comp.index, lw=lw, label='obs')
            
        # simulated time series
        sim_comp = sim[comp] #getattr(sim, comp)[id]
        ax.plot(sim_comp, sim_comp.index, lw=lw, label='sim')
            
        if obs is not None:        
            # add performance as text
            performance = KGEmod(obs_comp, sim_comp)
            performance = ['-∞' if x < -10 else '∞' if x > 10 else np.round(x, 2) for x in performance]
            ax.text(.02, .99, "KGE' = {0}\nα = {1}\nβ = {2}\nρ = {3}".format(*performance),
                    va='top', transform=ax.transAxes, bbox=bbox_props)
        
        # settings
        if comp in ['original', 'trend']:
            if lims is not None:
                for lim in lims:
                    ax.axvline(lim, ls=':', c='k', lw=.5)
            ax.set(xlim=xlim[0],
                   ylim=(sim_comp.index.max(), sim_comp.index.min()))
        else:
            ax.set(xlim=xlim[1])
        ax.set(xlabel=xlabel,
               title=comp);
        
    fig.legend(*ax.get_legend_handles_labels(), ncol=2, loc=8, frameon=False, bbox_to_anchor=[.4, -.025, .2, .1])
    if 'title' in kwargs:
        fig.text(.5, 1, kwargs['title'], fontsize=11, ha='center');

    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight');
        
        

def storage_outflow(storage: pd.Series, outflow: pd.Series, storage2: pd.Series = None, outflow2: pd.Series = None, s_lims: List = None, q_lims: List = None, save: Union[Path, str] = None, **kwargs):
    """It creates a figure that compares the storage and outflow time series. The figure is composed of three plots. In the center, a scatter plot of storage versus outflow; if the storage and outflow limits are provided, a line
    represents the reference LISFLOOD routine. On top, a plot shows the density function (kernel density estimation) of storage. On the right, a plot shows the density function (kernel density estimation) of outflow.
    
    Parameters:
    -----------
    storage:   pd.Series
        Series of reservoir storage. By default, it should be storage relative to the total reservoir capacity. It is supposed to be simulated storage
    outflow:   pd.Series
        Series of reservoir outflow. By default, it should be relative to the non-damaging outflow. It is supposed to be simulated outflow
    storage2:  pd.Series
        Series of reservoir storage. By default, it should be storage relative to the total reservoir capacity. It is supposed to be observed storage
    outflow2:  pd.Series
        Series of reservoir outflow. By default, it should be relative to the non-damaging outflow. It is supposed to be observed outflow
    s_lims:    List
        Storage limits (conservative, 2 times conservation, normal, adjusted normal, flood) used in the LISFLOOD reservoir routine
    q_lims:    List
        Outflow limits (minimum, minimum, normal adjusted, normal adjusted, non-damaging) used in the LISFLOOD reservoir routine
    save:      Union[str, Path]
        Path where to save the figure
    
    Keyword arguments:
    ------------------
    alpha:     float
        Transparency in the scatter plot
    color:     list(2,)
        Colours to be used in the simulated and observed data
    size:      float
        Point size in the scatter plot
    title:     str
        Title of the figure
    xlabel:    str
        Label of the X axis in the scatter plot
    ylabel:    str
        Label of the Y axis in the scatter plot
    ymin:      float
        Minimum value of the Y axis in the scatter plot
    """
    
    # extract kwargs
    s = kwargs.get('size', .5)
    a = kwargs.get('alpha', .05)
    c = kwargs.get('color', ['C1', 'C0'])
    ymin = kwargs.get('ymin', -.1)
    xlabel = kwargs.get('xlabel', 'relative storage (-)')
    ylabel = kwargs.get('ylabel', 'relative outflow (-)')
    
    if storage2 is not None:
        if storage2.isnull().all():
            storage2 = None
    if outflow2 is not None:
        if outflow2.isnull().all():
            outflow2 = None
    
    if (s_lims is not None) & (q_lims is not None):
        assert len(s_lims) == len(q_lims), 'The length of "s_lims" and "q_lims" must be the same.'
    
    # Create the figure and set the size
    fig = plt.figure(figsize=(5, 5))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 4], width_ratios=[4, 1])
    if 'title' in kwargs:
        fig.text(.95, .95, kwargs['title'], ha='right', va='top')  

    # scatter plot outflow vs storage
    ax1 = plt.subplot(gs[1, 0])
    ax1.scatter(storage, outflow, c=c[0], s=s, alpha=a)
    ax1.scatter(-1, -1, c=c[0], s=s, label=kwargs.get('label1', 'sim'))
    if (storage2 is not None) & (outflow2 is not None):
        ax1.scatter(storage2, outflow2, c=c[1], s=s, alpha=a)
        ax1.scatter(-1, -1, s=s, c=c[1], label=kwargs.get('label2', 'obs'))
    if (s_lims is not None) & (q_lims is not None):
        ax1.plot(s_lims, q_lims, c='k', lw=1, zorder=0, label='routine');
        for s, q in zip(s_lims, q_lims):
            ax1.hlines(q, xmin=-.02, xmax=s, color='k', ls=':', lw=.5, zorder=10)
            ax1.vlines(s, ymin=ymin, ymax=q, color='k', ls=':', lw=.5, zorder=10)
    ax1.set(ylim= (ymin, None),
            xlim=(-.02, 1.02), 
            xlabel=xlabel,
            ylabel=ylabel)
    ax1.spines[['top', 'right']].set_visible(False)
    ax1.legend(frameon=False, loc=2)
    
    # densidy distribution: storage
    ax2 = plt.subplot(gs[0, 0])
    sns.kdeplot(storage, color=c[0], fill=True, ax=ax2)
    if storage2 is not None:
        sns.kdeplot(storage2, color=c[1], fill=True, ax=ax2)
        kge = KGEmod(storage, storage2)[0]
        ax2.text(.5, -.2, f"KGE' = {kge:.2f}", fontsize=9, va='center', ha='center', transform=ax2.transAxes)
    if s_lims is not None:
        for s in s_lims:
            ax2.axvline(s, color='k', ls=':', lw=.5, zorder=10)
    ax2.spines[['top', 'left', 'right']].set_visible(False)
    ax2.set(xlim=(-.02, 1.02),
            ylabel=None,
            xlabel=None)
    ax2.set_yticks([])
    ax2.set_xticklabels([])

    # density distribution: outflow
    ax3 = plt.subplot(gs[1, 1])
    sns.kdeplot(y=outflow, color=c[0], fill=True, ax=ax3)
    if outflow2 is not None:
        sns.kdeplot(y=outflow2, color=c[1], fill=True, ax=ax3)
        kge = KGEmod(outflow, outflow2)[0]
        ax3.text(-0.2, .55, f"KGE' = {kge:.2f}", fontsize=9, va='center', ha='center', rotation=90, transform=ax3.transAxes)
    if q_lims is not None:
        for q in q_lims:
            ax3.axhline(q, color='k', ls=':', lw=.5, zorder=10)
    ax3.spines[['top', 'bottom', 'right']].set_visible(False)
    ax3.set(ylim=(ymin, None),
            xlabel=None,
            ylabel=None)
    ax3.set_xticks([])
    ax3.set_yticklabels([])

    # Adjust the spacing between subplots
    fig.tight_layout();
    
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')