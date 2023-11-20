import matplotlib as mpl
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import geopandas as gpd
import seaborn as sns
import cartopy.crs as ccrs
import cartopy.feature as cf
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Union, Dict, List, Tuple

# from utils import Decomposition
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
    c = kwargs.get('color', ['C1', 'C0'])
    
    # setup
    bbox_props = dict(boxstyle='round, pad=0.05', facecolor='w', edgecolor='none', alpha=.666)
    ncols = 4
    fig, axes = plt.subplots(figsize=(figsize[0] * ncols, figsize[1]), ncols=ncols, tight_layout=True, sharey=True)#, sharex=True)
    
    components = [attr for attr in dir(obs) if not attr.startswith('_')]
    for ax, comp in zip(axes, sim.columns):
        if obs is not None:
            # observed time series
            obs_comp = obs[comp] #getattr(obs, comp)[id]
            ax.plot(obs_comp, obs_comp.index, c=c[1], lw=lw, label='obs')
            
        # simulated time series
        sim_comp = sim[comp] #getattr(sim, comp)[id]
        ax.plot(sim_comp, sim_comp.index, c=c[0], lw=lw, label='sim')
            
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
        
        

# def storage_outflow(storage: pd.Series, outflow: pd.Series, storage2: pd.Series = None, outflow2: pd.Series = None, s_lims: List = None, q_lims: List = None, save: Union[Path, str] = None, **kwargs):
#     """It creates a figure that compares the storage and outflow time series. The figure is composed of three plots. In the center, a scatter plot of storage versus outflow; if the storage and outflow limits are provided, a line
#     represents the reference LISFLOOD routine. On top, a plot shows the density function (kernel density estimation) of storage. On the right, a plot shows the density function (kernel density estimation) of outflow.
    
#     Parameters:
#     -----------
#     storage:   pd.Series
#         Series of reservoir storage. By default, it should be storage relative to the total reservoir capacity. It is supposed to be simulated storage
#     outflow:   pd.Series
#         Series of reservoir outflow. By default, it should be relative to the non-damaging outflow. It is supposed to be simulated outflow
#     storage2:  pd.Series
#         Series of reservoir storage. By default, it should be storage relative to the total reservoir capacity. It is supposed to be observed storage
#     outflow2:  pd.Series
#         Series of reservoir outflow. By default, it should be relative to the non-damaging outflow. It is supposed to be observed outflow
#     s_lims:    List
#         Storage limits (conservative, 2 times conservation, normal, adjusted normal, flood) used in the LISFLOOD reservoir routine
#     q_lims:    List
#         Outflow limits (minimum, minimum, normal adjusted, normal adjusted, non-damaging) used in the LISFLOOD reservoir routine
#     save:      Union[str, Path]
#         Path where to save the figure
    
#     Keyword arguments:
#     ------------------
#     alpha:     float
#         Transparency in the scatter plot
#     color:     list(2,)
#         Colours to be used in the simulated and observed data
#     size:      float
#         Point size in the scatter plot
#     title:     str
#         Title of the figure
#     xlabel:    str
#         Label of the X axis in the scatter plot
#     ylabel:    str
#         Label of the Y axis in the scatter plot
#     ymin:      float
#         Minimum value of the Y axis in the scatter plot
#     """
    
#     # extract kwargs
#     s = kwargs.get('size', .5)
#     a = kwargs.get('alpha', .05)
#     c = kwargs.get('color', ['C1', 'C0'])
#     ymin = kwargs.get('ymin', -.1)
#     xlabel = kwargs.get('xlabel', 'relative storage (-)')
#     ylabel = kwargs.get('ylabel', 'relative outflow (-)')
    
#     if storage2 is not None:
#         if storage2.isnull().all():
#             storage2 = None
#     if outflow2 is not None:
#         if outflow2.isnull().all():
#             outflow2 = None
    
#     if (s_lims is not None) & (q_lims is not None):
#         assert len(s_lims) == len(q_lims), 'The length of "s_lims" and "q_lims" must be the same.'
    
#     # Create the figure and set the size
#     fig = plt.figure(figsize=(5, 5))
#     gs = gridspec.GridSpec(2, 2, height_ratios=[1, 4], width_ratios=[4, 1])
#     if 'title' in kwargs:
#         fig.text(.95, .95, kwargs['title'], ha='right', va='top')  

#     # scatter plot outflow vs storage
#     ax1 = plt.subplot(gs[1, 0])
#     ax1.scatter(storage, outflow, c=c[0], s=s, alpha=a)
#     ax1.scatter(-1, -1, c=c[0], s=s, label=kwargs.get('label1', 'sim'))
#     if (storage2 is not None) & (outflow2 is not None):
#         ax1.scatter(storage2, outflow2, c=c[1], s=s, alpha=a)
#         ax1.scatter(-1, -1, s=s, c=c[1], label=kwargs.get('label2', 'obs'))
#     if (s_lims is not None) & (q_lims is not None):
#         ax1.plot(s_lims, q_lims, c='k', lw=1, zorder=0, label='routine');
#         for s, q in zip(s_lims, q_lims):
#             ax1.hlines(q, xmin=-.02, xmax=s, color='k', ls=':', lw=.5, zorder=10)
#             ax1.vlines(s, ymin=ymin, ymax=q, color='k', ls=':', lw=.5, zorder=10)
#     ax1.set(ylim= (ymin, None),
#             xlim=(-.02, 1.02), 
#             xlabel=xlabel,
#             ylabel=ylabel)
#     ax1.spines[['top', 'right']].set_visible(False)
#     ax1.legend(frameon=False, loc=2)
    
#     # densidy distribution: storage
#     ax2 = plt.subplot(gs[0, 0])
#     sns.kdeplot(storage, color=c[0], fill=True, ax=ax2)
#     if storage2 is not None:
#         sns.kdeplot(storage2, color=c[1], fill=True, ax=ax2)
#         kge = KGEmod(storage, storage2)[0]
#         ax2.text(.5, -.2, f"KGE' = {kge:.2f}", fontsize=9, va='center', ha='center', transform=ax2.transAxes)
#     if s_lims is not None:
#         for s in s_lims:
#             ax2.axvline(s, color='k', ls=':', lw=.5, zorder=10)
#     ax2.spines[['top', 'left', 'right']].set_visible(False)
#     ax2.set(xlim=(-.02, 1.02),
#             ylabel=None,
#             xlabel=None)
#     ax2.set_yticks([])
#     ax2.set_xticklabels([])

#     # density distribution: outflow
#     ax3 = plt.subplot(gs[1, 1])
#     sns.kdeplot(y=outflow, color=c[0], fill=True, ax=ax3)
#     if outflow2 is not None:
#         sns.kdeplot(y=outflow2, color=c[1], fill=True, ax=ax3)
#         kge = KGEmod(outflow, outflow2)[0]
#         ax3.text(-0.2, .55, f"KGE' = {kge:.2f}", fontsize=9, va='center', ha='center', rotation=90, transform=ax3.transAxes)
#     if q_lims is not None:
#         for q in q_lims:
#             ax3.axhline(q, color='k', ls=':', lw=.5, zorder=10)
#     ax3.spines[['top', 'bottom', 'right']].set_visible(False)
#     ax3.set(ylim=(ymin, None),
#             xlabel=None,
#             ylabel=None)
#     ax3.set_xticks([])
#     ax3.set_yticklabels([])

#     # Adjust the spacing between subplots
#     fig.tight_layout();
    
#     if save is not None:
#         plt.savefig(save, dpi=300, bbox_inches='tight')
        
        
        
def reservoir_scatter(sim: pd.DataFrame, x: str, y: str, obs: pd.DataFrame = None, x_thr: List = None, y_thr: List = None, legend: bool = True, ax: Axes = None, **kwargs):
    """It creates a figure that compares the storage and outflow time series. The figure is composed of three plots. In the center, a scatter plot of storage versus outflow; if the storage and outflow limits are provided, a line
    represents the reference LISFLOOD routine. On top, a plot shows the density function (kernel density estimation) of storage. On the right, a plot shows the density function (kernel density estimation) of outflow.
    
    Parameters:
    -----------
    sim:   pd.DataFrame
        Simulated time series of reservoir behaviour. It should contain colums "x" and "y"
    x:     str
        Column of "sim" (and "obs") to be used in the X axis
    y:     str
        Column of "sim" (and "obs") to be used in the Y axis
    obs:   pd.DataFrame
        Oberved time series of reservoir behaviour. It should contain colums "x" and "y"
    x_thr:    List
        Thresholds in the LISFLOOD reservoir routine to be used in the X axis
    y_thr:    List
        Thresholds in the LISFLOOD reservoir routine to be used in the Y axis
    legend:   bool
        Whether to plot the legend or not
    ax:       Axes
        Matplotlib axes in which to insert the plot
    
    Keyword arguments:
    ------------------
    alpha:     float
        Transparency in the scatter plot
    color:     List(2,)
        Colours to be used in the simulated and observed data
    figsize    Tuple(2,)
        Size of the figure
    labels:    List(2,)
        Labels of the datasets "sim" and "obs"
    size:      float
        Point size in the scatter plot
    xlabel:    str
        Label of the X axis in the scatter plot
    xlim:      Tuple(2,)
        Limits of the X axis
    xticklabels: bool
        Whether to include values in the X ticks or not
    ylabel:    str
        Label of the Y axis in the scatter plot
    ylim:      Tuple(2,)
        Limits of the Y axis
    yticklabels: bool
        Whether to include values in the Y ticks or not
    """
    
    # extract kwargs
    a = kwargs.get('alpha', .05)
    c = kwargs.get('color', ['C1', 'C0'])
    figsize = kwargs.get('figsize', (4, 4))
    labels = kwargs.get('labels', ['sim', 'obs'])
    s = kwargs.get('size', .5)
    xlabel = kwargs.get('xlabel', x)
    xlim = kwargs.get('xlim', (-.1, None))
    xticklabels = kwargs.get('xticklabels', True)
    ylabel = kwargs.get('ylabel', y)
    ylim = kwargs.get('ylim', (-.1, None))
    yticklabels = kwargs.get('yticklabels', True)
    
    if (x_thr is not None) & (y_thr is not None):
        assert len(x_thr) == len(y_thr), 'The length of "x_thr" and "y_thr" must be the same.'
    
    # Create the figure and set the size
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # scatter plot outflow vs storage
    ax.scatter(sim[x], sim[y], c=c[0], s=s, alpha=a, zorder=2)
    ax.scatter(-1, -1, c=c[0], s=s, label=labels[0])
    if obs is not None:
        ax.scatter(obs[x], obs[y], c=c[1], s=s, alpha=a, zorder=1)
        ax.scatter(-1, -1, s=s, c=c[1], label=labels[1])
    if (x_thr is not None) & (y_thr is not None):
        #ax.plot(x_thr, y_thr, c='k', lw=1, zorder=0, label='routine');
        for s, q in zip(x_thr, y_thr):
            ax.hlines(q, xmin=-.02, xmax=s, color='k', ls=':', lw=.5, zorder=10)
            ax.vlines(s, ymin=ylim[0], ymax=q, color='k', ls=':', lw=.5, zorder=10)
    ax.set(xlim=xlim,
           ylim=ylim,
           xlabel=xlabel,
           ylabel=ylabel)
    if xticklabels is False:
        ax.set_xticklabels([])
    if yticklabels is False:
        ax.set_yticklabels([])
    ax.spines[['top', 'right']].set_visible(False)
    
    if legend:
        ax.legend(frameon=False, loc=2)

        
        
def reservoir_kde(sim: pd.DataFrame, obs: pd.DataFrame = None, x: str = None, y: str = None, thr: List = None, ax: Axes = None, **kwargs):
    """It creates a figure that compares the storage and outflow time series. The figure is composed of three plots. In the center, a scatter plot of storage versus outflow; if the storage and outflow limits are provided, a line
    represents the reference LISFLOOD routine. On top, a plot shows the density function (kernel density estimation) of storage. On the right, a plot shows the density function (kernel density estimation) of outflow.
    
    Parameters:
    -----------
    sim:   pd.DataFrame
        Simulated time series of reservoir behaviour. It should contain colums "x" and "y"
    obs:   pd.DataFrame
        Oberved time series of reservoir behaviour. It should contain colums "x" and "y"
    x:     str
        Column of "sim" (and "obs") for which the distribution will be computed. The plot will be horizontal
    y:     str
        Column of "sim" (and "obs") for which the distribution will be computed. The plot will be vertical
    thr:    List
        Thresholds in the LISFLOOD reservoir routine
    ax:       Axes
        Matplotlib axes in which to insert the plot
    
    Keyword arguments:
    ------------------
    color:     List(2,)
        Colours to be used in the simulated and observed data
    figsize    Tuple(2,)
        Size of the figure
    xlabel:    str
        Label of the X axis in the scatter plot
    xlim:      Tuple(2,)
        Limits of the X axis
    xticklabels: bool
        Whether to include values in the X ticks or not
    ylabel:    str
        Label of the Y axis in the scatter plot
    ylim:      Tuple(2,)
        Limits of the Y axis
    yticklabels: bool
        Whether to include values in the Y ticks or not
    """
    
    assert (x != y) & ((x is not None) | (y is not None)), 'Either "x" or "y" must indicate a column in "sim".'
    
    # extract kwargs
    c = kwargs.get('color', ['C1', 'C0'])
    figsize = kwargs.get('figsize', (4, 4))
    xlim = kwargs.get('xlim', (-.1, None))
    xlabel = kwargs.get('xlabel', None)
    xticklabels = kwargs.get('xticklabels', True)
    ylabel = kwargs.get('ylabel', None)
    ylim = kwargs.get('ylim', (-.1, None))
    yticklabels = kwargs.get('yticklabels', True)
    
    # Create the figure and set the size
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)  
    
    # density distribution
    if x is not None:
        sns.kdeplot(sim[x], color=c[0], fill=True, ax=ax)
        if obs is not None:
            if not obs[x].isnull().all():
                sns.kdeplot(obs[x], color=c[1], fill=True, ax=ax)
                kge = KGEmod(obs[x], sim[x])[0]
                ax.text(.5, -.2, f"KGE' = {kge:.2f}", fontsize=9, va='center', ha='center', transform=ax.transAxes)
        if thr is not None:
            for x in thr:
                ax.axvline(x, color='k', ls=':', lw=.5, zorder=10)
        ax.spines[['top', 'left', 'right']].set_visible(False)
        ax.set(xlim=xlim,
                 ylabel=None,
                 xlabel=None)
        ax.set_yticks([])
        if xticklabels is False:
            ax.set_xticklabels([])
    elif y is not None:
        sns.kdeplot(y=sim[y], color=c[0], fill=True, ax=ax)
        if obs is not None:
            if not obs[y].isnull().all():
                sns.kdeplot(y=obs[y], color=c[1], fill=True, ax=ax)
                kge = KGEmod(obs[y], sim[y])[0]           
                ax.text(-0.2, .55, f"KGE' = {kge:.2f}", fontsize=9, va='center', ha='center', rotation=90, transform=ax.transAxes)
        if thr is not None:
            for y in thr:
                ax.axhline(y, color='k', ls=':', lw=.5, zorder=10)
        ax.spines[['top', 'bottom', 'right']].set_visible(False)
        ax.set(ylim=ylim,
               xlabel=None,
               ylabel=None)
        ax.set_xticks([])
        if yticklabels is False:
            ax.set_yticklabels([])

            
            
def reservoir_analysis(sim: pd.DataFrame, obs: pd.DataFrame = None, x1: str = 'storage', x2: str = 'inflow', y: str = 'outflow', x_thr: List = None, y_thr: List = None, save: Union[Path, str] = None, **kwargs):
    """It creates a figure that compares the storage and outflow time series. The figure is composed of three plots. In the center, a scatter plot of storage versus outflow; if the storage and outflow limits are provided, a line
    represents the reference LISFLOOD routine. On top, a plot shows the density function (kernel density estimation) of storage. On the right, a plot shows the density function (kernel density estimation) of outflow.
    
    Parameters:
    -----------
    sim:   pd.DataFrame
        Simulated time series of reservoir behaviour. It should contain colums "x" and "y"
    obs:   pd.DataFrame
        Oberved time series of reservoir behaviour. It should contain colums "x" and "y"
    x1:     str
        Column of "sim" (and "obs") that will be used in the X axis of the first scatter plot
    x2:     str
        Column of "sim" (and "obs") that will be used in the X axis of the second scatter plot
    y:     str
        Column of "sim" (and "obs") that will be used in the Y axis of both scatter plots
    x_thr:    List
        Thresholds in the LISFLOOD reservoir routine to be used in the "x1" variable
    y_thr:    List
        Thresholds in the LISFLOOD reservoir routine to be used in the "y" axis
    save:      Union[str, Path]
        Path where to save the figure
    
    Keyword arguments:
    ------------------
    alpha:     float
        Transparency in the scatter plot
    color:     list(2,)
        Colours to be used in the simulated and observed data
    figsize:   Tuple(2,)
        Size of the figure
    labels:    List(2,)
        Labels of the datasets "sim" and "obs"
    size:      float
        Point size in the scatter plot
    title:     str
        Title of the figure
    """
    
    # extract kwargs
    a = kwargs.get('alpha', .05)
    c = kwargs.get('color', ['C1', 'C0'])
    figsize = kwargs.get('figsize', (9, 5))
    labels = kwargs.get('labels', ['sim', 'obs'])
    s = kwargs.get('size', .5)
    x1_lim = kwargs.get('x1lim', (-.02, 1.02))
    
    if (x_thr is not None) & (y_thr is not None):
        assert len(x_thr) == len(y_thr), 'The length of "x_thr" and "y_thr" must be the same.'
    
    # Create the figure and set the size
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 3, height_ratios=[1, 4], width_ratios=[4, 4, 1])
    if 'title' in kwargs:
        fig.text(.95, .95, kwargs['title'], ha='right', va='top')  

    # scatter plot: x1 vs y
    ax10 = plt.subplot(gs[1, 0])
    reservoir_scatter(sim, x1, y, obs, x_thr=x_thr, y_thr=y_thr, xlim=x1_lim, ax=ax10, legend=False,
                      size=s, alpha=a, color=c, labels=labels)
    
    # scatter plot: x2 vs y
    ax11 = plt.subplot(gs[1, 1])
    reservoir_scatter(sim, x2, y, obs, x_thr=y_thr, y_thr=y_thr, ax=ax11, legend=False, ylim=ax10.get_ylim(), ylabel='', yticklabels=False,
                      size=s, alpha=a, color=c, labels=labels)
    ax11.plot(ax10.get_ylim(), ax10.get_ylim(), c='k', lw=.5, ls=':', zorder=0)
    
    # density distribution: x1
    ax00 = plt.subplot(gs[0, 0])
    reservoir_kde(sim, obs, x=x1, thr=x_thr, ax=ax00, xlim=ax10.get_xlim(), xticklabels=False,
                  color=c)
    
    # density distribution: x2
    ax01 = plt.subplot(gs[0, 1])
    reservoir_kde(sim, obs, x=x2, thr=y_thr, ax=ax01, xlim=ax11.get_xlim(), xticklabels=False,
                  color=c)

    # density distribution: y
    ax12 = plt.subplot(gs[1, 2])
    reservoir_kde(sim, obs, y=y, thr=y_thr, ax=ax12, ylim=ax10.get_ylim(), yticklabels=False)
    
    fig.legend(*ax10.get_legend_handles_labels(), frameon=False, ncol=3, loc=8, bbox_to_anchor=[.25, -.04, .5, .05])
    
    # Adjust the spacing between subplots
    fig.tight_layout();
    
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        
        
        
def maps_performance(x: pd.Series, y: pd.Series, performance: pd.DataFrame, s: Union[pd.Series, int] = None, polygons: gpd.GeoDataFrame = None, save: Union[Path, str] = None, **kwargs):
    """It creates a figure that contains 4 maps with the KGE and its 3 components.
    
    Inputs:
    -------
    x: pd.Series
        X coordinate of the points 
    y: pd.Series
        Y coordinate of the points
    performance: pd.DataFrame
        Performance of each of the points. It must contain at least 4 columns: 'KGE', 'r', 'alpha', 'beta'
    s: Union[pd.Series, int]
        Size of each point. In case of an integer all points will have the same size. In case of a pd.Series every point will have an specific size
    save: Union[Path, str]
        File in which the plot will be saved
        
    kwargs:
    -------
    figsize: Tuple
        Size of the figure
    proj: 
        Projection of the map
    extent: List
        Extension of the map [xmin, xmax, ymin, ymax]
    alpha: float
        Transparency of the points
    title: str
        Title of the figure
    """
    
    figsize = kwargs.get('figsize', (15, 7))
    proj = kwargs.get('proj', ccrs.PlateCarree())
    extent = kwargs.get('extent', None)#[-125, -67.5, 24, 51]
    alpha = kwargs.get('alpha', .8)
    if s is None:
        s = 5
    
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=figsize, tight_layout=True, subplot_kw={'projection': proj})

    for i, metric in enumerate(['KGE', 'r', 'alpha', 'beta']):

        if metric == 'KGE':
            cmap, norm = create_cmap('Spectral', [-100, -1, -.75, -.5, -.25 ,0, .25, .5, .75, 1])
        elif metric in ['alpha', 'beta']:
            cmap, norm = create_cmap('RdBu', [1e-6, 1/16, 1/8, 1/4, 1/2, 1/1.2, 1.2, 2, 4, 8, 16, 1e6])
        elif metric == 'r':
            cmap, norm = create_cmap('Spectral', [-1, -.75, -.5, -.25, 0, .25, .5, .75, 1])

        # background map
        ax = axes[int(i / 2), i % 2]
        ax.add_feature(cf.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='lightgray'), alpha=.5, zorder=0)
        if polygons is not None:
            polygons.plot(facecolor='none', edgecolor='white', ax=ax)
        if extent is not None:
            ax.set_extent(extent)
        ax.axis('off')

        # scatter plot
        sct = ax.scatter(x, y, c=performance[metric], cmap=cmap, norm=norm, edgecolor='w', lw=1, 
                          s=s, alpha=alpha)
        # # setup: color bar, title...
        cbar = plt.colorbar(sct, ax=ax, shrink=.66)#, orientation='horizontal')
        cbar.set_label(metric, rotation=90)
    
    if 'title' in kwargs:
        fig.text(.5, 1.0, kwargs['title'], ha='center', va='bottom', fontsize=12);
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight')
        
        
        
def plot_iterations(iters: pd.DataFrame, pareto: pd.DataFrame, best_iter: int, cols: List = ['like1', 'like2'], save: Union[str, Path] = None, **kwargs):
    """It creates a scatter plot that shows the performance of the iterations in the calibration. On top of the scatter plot a line depicts the Pareto front, from which the best iteration is taken.
    
    Inputs:
    -------
    iters: pd.DataFrame
        A table that contains the two objective functions ("cols") to be plotted
    pareto: pd.DataFrame
        A table that contains the Pareto front from "iters[cols]"
    best_iter: int
        The index of "iters" that contains the best iteration
    save: Union[str, Path]
        If provided, where to save the plot
    """
    
    figsize = kwargs.get('figsize', (4, 4))
    vmax = kwargs.get('vmax', None)
    xlabel = kwargs.get('xlabel', r'$L_{outflow}$')
    ylabel = kwargs.get('ylabel', r'$L_{storage}$')
    title = kwargs.get('title', None)
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(iters[cols[0]], iters[cols[1]], s=1, c='gray', label='iteration', alpha=.2)
    ax.scatter(*iters.loc[best_iter, cols], c='steelblue', label='optimum', s=4)
    ax.plot(pareto[cols[0]], pareto[cols[1]], c='k', lw=1, ls=':', label='pareto front', zorder=0)
    ax.set(xlim=(-.025, vmax),
           xlabel=xlabel,
           ylim=(-.025, vmax),
           ylabel=ylabel)
    if 'title' in kwargs:
        ax.set_title(kwargs['title'])
    fig.legend(frameon=False, loc=1, bbox_to_anchor=[1.175, .7, .1, .2])
    if save is not None:
        plt.savefig(save, dpi=300, bbox_inches='tight');