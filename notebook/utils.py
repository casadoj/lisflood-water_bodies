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
import xml.etree.ElementTree as ET



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
        
        

# def decomposition(data: Union[pd.DataFrame, pd.Series]) -> Tuple[Union[pd.DataFrame, pd.Series], Union[pd.DataFrame, pd.Series], Union[pd.DataFrame, pd.Series]]:
#     """It decomposes the timeseries in three components: annual average, seasonality and residuals.
    
#     Parameters:
#     -----------
#     data: Union[pd.DataFrame, pd.Series]
#         Time series to be decomposed
        
#     Returns:
#     --------
#     annual: Union[pd.DataFrame, pd.Series]
#         Timeseries of mean annual values. The length of the timeseries will be the number of years in "data"
#     seasonality: Union[pd.DataFrame, pd.Series]
#         Timeseries of montly mean values after removal of the annual trend. The length of the timeseries will be alwas 12, the months
#     residuals: Union[pd.DataFrame, pd.Series]
#         Timeseries of residuals, that is, the difference of the original "data" and the annual and seosonal timeseris. The length of the timeseries is the same as the input "data"
#     """
    
#     # annual average storage
#     annual = data.resample('Y').mean()
#     annual.index = annual.index.year

#     # seasonal variability
#     detrended = data - data.resample('Y').transform('mean')
#     seasonality = detrended.groupby(detrended.index.month).mean()

#     # residual
#     residuals = detrended - detrended.groupby(detrended.index.month).transform('mean')
    
#     return annual, seasonality, residuals



# class Decomposition:
#     def __init__(self, original: Tuple[Union[pd.DataFrame, pd.Series]], trend: Tuple[Union[pd.DataFrame, pd.Series]], seasonal: Tuple[Union[pd.DataFrame, pd.Series]], residuals: Tuple[Union[pd.DataFrame, pd.Series]]):
#         self.original = original
#         self.trend = trend
#         self.seasonal = seasonal
#         self.residual = residuals

        
        
# def decompose_timeseries(data: Union[pd.DataFrame, pd.Series], window: int = 365, center: bool = True) -> Decomposition:
#     """It decomposes the timeseries in three components: trend, seasonality and residuals.
    
#     Parameters:
#     -----------
#     data: Union[pd.DataFrame, pd.Series]
#         Time series to be decomposed
        
#     Returns:
#     --------
#     DecompositionResult:
#         Object with three methods: trend(), seasonal(), residuals()
#     """ 

#     # trend as the 365 rolling mean
#     trend = data.rolling(window=365, min_periods=180, center=center).mean()

#     # seasonality
#     detrended = data - trend
#     seasonal = detrended.groupby(detrended.index.month).transform('mean')

#     # residuals
#     residual = detrended - seasonal

#     return Decomposition(data, trend, seasonal, residual)



def decompose_timeseries(data: pd.Series, window: int = 365, center: bool = True) -> pd.DataFrame:
    """It decomposes the timeseries in three components: trend, seasonality and residuals.
    
    Parameters:
    -----------
    data: pd.Series
        Time series to be decomposed
        
    Returns:
    --------
    DecompositionResult:
        Object with three methods: trend(), seasonal(), residuals()
    """ 
    
    assert isinstance(data, pd.Series), '"data" must be a pandas.Series'
    
    # trend as the 365 rolling mean
    trend = data.rolling(window=365, min_periods=180, center=center).mean()

    # seasonality
    detrended = data - trend
    seasonal = detrended.groupby(detrended.index.month).transform('mean')

    # residuals
    residual = detrended - seasonal
    
    decomposition = pd.concat((data, trend, seasonal, residual), axis=1)
    decomposition.columns = ['original', 'trend', 'seasonal', 'residual']
    
    return decomposition



def xml_parameters(xml: Union[str, Path], pars: Union[str, List[str]] = None) -> Dict:
    """It extracts the temporal information from the settings XML file.
    
    Input:
    ------
    xml:         Union[str, Path] 
        A XML settings file (path, filename and extension)
    pars:        Union[str, List[str]]
        Name(s) of the parameters to be extracted
        
    Output:
    -------
    parameters:  Dict
        Keys are parameter names and values the calibrated parameter value
    """
    
    # extract temporal info from the XML
    tree = ET.parse(xml)
    root = tree.getroot()
    
    if pars is None:
        pars = ['b_Xinanjiang', 'UpperZoneTimeConstant', 'LowerZoneTimeConstant', 'LZThreshold',
                'GwPercValue', 'GwLoss', 'PowerPrefFlow', 'SnowMeltCoef',
                'AvWaterRateThreshold' , 'LakeMultiplier', 'adjust_Normal_Flood', 'ReservoirRnormqMult', 
                'QSplitMult', 'CalChanMan', 'CalChanMan2', 'ChanBottomWMult', 'ChanDepthTMult', 'ChanSMult']
    
    parameters = {par: float(root.find(f'.//textvar[@name="{par}"]').attrib['value']) for par in pars}
        
    return parameters