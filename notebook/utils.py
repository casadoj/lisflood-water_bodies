import os
os.environ['USE_PYGEOS'] = '0'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from typing import Union, List, Dict, Tuple
from pathlib import Path



def select_reservoirs(df: pd.DataFrame, sort: str, storage: str, target: float, plot: bool = True, **kwargs) -> pd.DataFrame:
    """It selects reservoirs that fulfil a target total storage capacity by prioritizing based on another characteristic
    
    Inputs:
    -------
    df:    pandas.DataFrame. Table of reservoirs
    sort:  string. Name of the field in 'df' that will be use to sort (prioritize) the selection
    storage: string. Name of the field in 'df' that contains the reservoir storage capacity
    plot:    boolean. If True, a map of the selected reservoirs will be plotted. The size of the dots represents the reservoir storage capacity and the colours the sorting field.
    
    Outputs:
    --------
    df_sel: pandas.DataFrame. A subset of 'df' with the selection of reservoirs.
    """
    
    mask = df.sort_values(sort, ascending=False)[storage].cumsum() <= target
    df_sel = df.loc[mask]
    volume = df_sel[storage].sum()
    
    if plot:
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (20, 5)), subplot_kw=dict(projection=ccrs.PlateCarree()))
        ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m', edgecolor='face', facecolor='lightgray'), alpha=.5, zorder=0)
        if 'c' in kwargs:
            if isinstance(kwargs['c'], str):
                c = kwargs['c']
            elif isinstance(kwargs['c'], pd.Series):
                c = kwargs['c'].loc[mask]
        else:
            c = df_sel[sort]
        scatter = ax.scatter(df_sel.geometry.x, df_sel.geometry.y, s=df_sel[storage] / 1000, cmap=kwargs.get('cmap', 'coolwarm'), c=c, alpha=kwargs.get('alpha', .5))
        if 'title' in kwargs:
            ax.text(.5, 1.07, kwargs['title'], horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes, fontsize=12)
        text = '{0} reservoirs   {1:.0f} km³'.format(mask.sum(), volume / 1000)
        ax.text(.5, 1.02, text, horizontalalignment='center', verticalalignment='bottom', transform=ax.transAxes)
        ax.axis('off');
        # if 'c' in kwargs:
        #     if isinstance(kwargs['c'], pd.Series):
        #         legend1 = ax.legend(*scatter.legend_elements(prop='colors', num=4, alpha=.5), title=kwargs.get('legend_title', ''), bbox_to_anchor=[1.025, .65, .09, .25], frameon=False)
        #         ax.add_artist(legend1)
        legend2 = ax.legend(*scatter.legend_elements(prop='sizes', num=4, alpha=.5), title='storage (km³)', bbox_to_anchor=[1.025, .35, .1, .25], frameon=False)
        ax.add_artist(legend2);
    
    return df_sel



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
        
        

def decomposition(data: Union[pd.DataFrame, pd.Series]) -> Tuple[Union[pd.DataFrame, pd.Series], Union[pd.DataFrame, pd.Series], Union[pd.DataFrame, pd.Series]]:
    """It decomposes the timeseries in three components: annual average, seasonality and residuals.
    
    Parameters:
    -----------
    data: Union[pd.DataFrame, pd.Series]
        Time series to be decomposed
        
    Returns:
    --------
    annual: Union[pd.DataFrame, pd.Series]
        Timeseries of mean annual values. The length of the timeseries will be the number of years in "data"
    seasonality: Union[pd.DataFrame, pd.Series]
        Timeseries of montly mean values after removal of the annual trend. The length of the timeseries will be alwas 12, the months
    residuals: Union[pd.DataFrame, pd.Series]
        Timeseries of residuals, that is, the difference of the original "data" and the annual and seosonal timeseris. The length of the timeseries is the same as the input "data"
    """
    
    # annual average storage
    annual = data.resample('Y').mean()
    annual.index = annual.index.year

    # seasonal variability
    detrended = data - data.resample('Y').transform('mean')
    seasonality = detrended.groupby(detrended.index.month).mean()

    # residual
    residuals = detrended - detrended.groupby(detrended.index.month).transform('mean')
    
    return annual, seasonality, residuals