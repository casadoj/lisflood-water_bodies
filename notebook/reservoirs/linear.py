import os
os.environ['USE_PYGEOS'] = '0'
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from tqdm.notebook import tqdm
from typing import Union, List, Tuple, Dict
from pathlib import Path

from plots import reservoir_analysis
from metrics import KGEmod

class Linear:
    """Representation of a reservoir according to Hanazaki, Yamazaki & Yoshimura (2021)."""
    
    def __init__(self, Vmin: float, Vtot: float, Qmin: float, T: int, At: int = 86400):
        """        
        Parameters:
        -----------
        Vmin: float
            Volume (m3) associated to the conservative storage
        Vtot: float
            Total reservoir storage capacity (m3)
        Qmin: float
            Minimum outflow (m3/s)
        T: int
            Residence time in days. The coefficient of the linear reservoir is the inverse of T (1/T)
        At: int
            Simulation time step in seconds.
        """
        
        # storage limits
        self.Vmin = Vmin
        self.Vtot = Vtot
        
        # outflow limits
        self.Qmin = Qmin
        
        # release coefficient
        self.k = 1 / (T * 24 * 3600)
        
        # time step duration in seconds
        self.At = At
        
    def timestep(self, I: float, V: float, verbose: bool = False) -> List[float]:
        """Given an inflow and an initial storage values, it computes the corresponding outflow
        
        Parameters:
        -----------
        I: float
            Inflow (m3/s)
        V: float
            Volume stored in the reservoir (m3)
        verbose: bool
            Whether to show on screen the evolution
            
        Returns:
        --------
        Q, V: List[float]
            Outflow (m3/s) and updated storage (m3)
        """
        
        # update reservoir storage with the inflow volume
        V += I * self.At
        
        # ouflow depending on the inflow and storage level
        Qmin = min(self.Qmin, (V - self.Vmin) / self.At)
        Q = max(Qmin, V * self.k, (V - self.Vtot) / self.At)
        
        # update reservoir storage with the outflow volume
        V -= Q * self.At
        
        return Q, V
        
    def simulate(self, inflow: pd.Series, Vo: float = None):
        """Given a inflow time series (m3/s) and an initial storage (m3), it computes the time series of outflow (m3/s) and storage (m3)
        
        Parameters:
        -----------
        inflow: pd.Series
            Time series of flow coming into the reservoir (m3/s)
        Vo: float
            Initial value of reservoir storage (m3). If not provided, it is assumed that the normal storage is the initial condition
            
        Returns:
        --------
        pd.DataFrame
            A table that concatenates the storage, inflow and outflow time series.
        """
        
        if Vo is None:
            Vo = self.Qtot * .5
        
        storage = pd.Series(index=inflow.index, dtype=float, name='storage')
        outflow = pd.Series(index=inflow.index, dtype=float, name='outflow')
        for ts in tqdm(inflow.index):
            # compute outflow and new storage
            Q, V = self.timestep(inflow[ts], Vo)
            storage[ts] = V
            outflow[ts] = Q
            # update current storage
            Vo = V

        return pd.concat((storage, inflow, outflow), axis=1)
    
    def routine(self, V: pd.Series, I: Union[float, pd.Series], modified: bool = True):
        """Given a time series of reservoir storage (m3) and a value or a time series of inflow (m3/s), it computes the ouflow (m3/s). This function is only meant for explanatory purposes; since the volume time series is given, the computed outflow does not update the reservoir storage. If the intention is to simulate the behaviour of the reservoir, refer to the function "simulate"
        
        Parameters:
        -----------
        V: pd.Series
            Time series of reservoir storage (m3)
        I: Union[float, pd.Series]
            Reservor inflow (m3/s)
        modified: bool
            Whether to use the modified (default) of the orinigal Hanazaki's routine. The modified routine avoids the breaks in the outflow function at Vmin, Vf and Ve.
            
        Returns:
        --------
        O: pd.Series
            Time series of reservoir outflow (m3/s)
        """

        if isinstance(I, float) or isinstance(I, int):
            assert I >= 0, '"I" must be a positive value'
            I = pd.Series(I, index=V.index)

        maskI = I < self.Qf
        maskV1 = V < self.Vmin
        maskV2 = (self.Vmin <= V) & (V < self.Vf)
        maskV3 = (self.Vf <= V) & (V < self.Ve)
        maskV4 = self.Ve <= V

        O = pd.Series(index=V.index, dtype=float)

        # if in the lower storage level
        O[maskV1] = self.Qn * V[maskV1] / self.Vf

        # if input below flood... 
        # ... and storage below emergency level
        mask = (maskI & maskV2) | (maskI & maskV3)
        if np.sum(mask) > 0:
            if modified:
                O[mask] = self.Vmin / self.Vf * self.Qn + ((V[mask] - self.Vmin) / (self.Ve - self.Vmin))**2 * (self.Qf - self.Vmin / self.Vf * self.Qn)
            else:
                O[mask] = .5 * self.Qn + ((V[mask] - self.Vmin) / (self.Ve - self.Vmin))**2 * (self.Qf - self.Qn)
        # ... and storage above emergency level
        O[maskI & maskV4] = self.Qf

        # if inflow over flood...
        # ... and storage in the normal level
        mask = ~maskI & maskV2
        if np.sum(mask) > 0:
            if modified:
                O[mask] = self.Vmin / self.Vf * self.Qn + (V[mask] - self.Vmin) / (self.Vf - self.Vmin) * (self.Qf - self.Vmin / self.Vf * self.Qn)
            else:
                O[mask] = .5 * self.Qn + (V[mask] - self.Vmin) / (self.Vf - self.Vmin) * (self.Qf - self.Qn)
        # ... and storage in flood level
        mask = ~maskI & maskV3
        if np.sum(mask) > 0:
            O[mask] = self.Qf + self.k * (V[mask] - self.Vf) / (self.Ve - self.Vf) * (I[mask] - self.Qf)
        # ... and storage in emergency level
        O[~maskI & maskV4] = I[~maskI & maskV4]

        return O
           
#    def normalize_timeseries(self, timeseries: pd.DataFrame) -> pd.DataFrame:
#        """It normalizes the timeseries using the total reservoir capacity and the outflow associated to the flood limit. In this way, the storage time series ranges between 0 and 1, and the inflow and outflow time series are in the order of units.
        
#        Parameters:
#        -----------
#        timeseries: pd.DataFrame
#            A table with three columns ('storage', 'inflow', 'outflow') with the time series of a reservoir
            
#        Returns:
#        --------
#        ts_norm: pd.DataFrame
#            Table similar to the original but with normalized values
#        """

#        ts_norm = timeseries.copy()
#        ts_norm.storage /= self.Vtot
#        ts_norm[['inflow', 'outflow']] /= self.Qf

#        return ts_norm
    
    def get_params(self):
        """It generates a dictionary with the reservoir paramenters in the Hanazaki model."""

        params = {'Vmin': self.Vmin,
                  'Vtot': self.Vtot,
                  'Qmin': self.Qmin,
                  'T': 1 / (self.k * 24 * 3600)}

        return params
        
    def scatter(self, series1: pd.DataFrame, series2: pd.DataFrame = None, norm: bool = True, save: Union[Path, str] = None, **kwargs):
        """It compares two reservoir timeseries (inflow, outflow and storage) using the function 'reservoir_analysis'. If only 1 time series is given, the plot will simply show the reservoir behaviour of that set of time series.
        
        Inputs:
        -------
        series1: pd.DataFrame
            A table with the time series of 'inflow', 'outflow' and 'storage'
        series2: pd.DataFrame
            A second table with the time series of 'inflow', 'outflow' and 'storage'
        norm: bool
            Whether to normalize or not the time series by the total reservoir capacity (storage) and the non-damaging flow (outflow and inflow)
        save: Union[Path, str]
            Directory and file where the figure will be saved
                    
        kwargs:
        -------
        title: str
            If provided, title of the figure
        labels: List[str]
            A list of 2 strings to be used as labels for each set of time series
        alpha: float
            The transparency of the scatter plot
        """
        
        # storage limits
        #Vlims = np.array([2 * self.Vmin, self.Vn, self.Vn_adj, self.Vf])
        
        # outflow limits
        #Qlims = np.array([self.Qmin, self.Qn, self.Qn, self.Qnd])
        
        if norm:
            series1_ = self.normalize_timeseries(series1)
            if series2 is not None:
                series2_ = self.normalize_timeseries(series2)
            Vlims /= self.Vtot
            Qlims /= self.Qnd
            x1lim = (-.02, 1.02)
        else:
            series1_ = series1
            if series2 is not None:
                series2_ = series2
            x1lim = (0, None)
        reservoir_analysis(series1_, series2_,
                           x_thr=None, y_thr=None,
                           title=kwargs.get('title', None),
                           labels=kwargs.get('labels', ['sim', 'obs']),
                           alpha=kwargs.get('alpha', .05),
                           x1lim=x1lim,
                           save=save)
        
    def lineplot(self, sim: Dict[str, pd.DataFrame], obs: pd.DataFrame = None, save: Union[Path, str] = None, **kwargs):
        """It plots the simulated time series of outflow and storage. If the observed time series is provided, it is plotted and the modified KGE shown.

        Input:
        ------
        sim: Dict[str, pd.DataFrame]
            A dictionary that contains the name and simulated time series in a pandas.DataFrame format. This DataFrame must have at least the columns 'outflow' and 'storage'
        obs: pd.DataFrame
            The oberved time series. This DataFrame must have at least the columns 'outflow' and 'storage'
        save: Union[Path, str]
            Directory and file where the figure will be saved
        """
    
        figsize = kwargs.get('figsize', (12, 6))
        lw = kwargs.get('lw', 1)
        
        fig, axes = plt.subplots(nrows=2, figsize=figsize, sharex=True)

        variables = {'outflow': {'unit': 'm3',
                                 'factor': 1,
                                 'thresholds': None}, #[self.Qmin, self.Qn, self.Qnd]},
                     'storage': {'unit': 'hm3',
                                 'factor': 1e-6,
                                 'thresholds': None}} #[self.Vmin, self.Vn, self.Vn_adj, self.Vf, self.Vtot]}}

        for ax, (var, dct) in zip(axes, variables.items()):
            f = dct['factor']
            if obs is not None:
                ax.plot(obs[var] * f, lw=.5 * lw, c='k', label='obs')
            for i, (label, serie) in enumerate(sim.items()):
                ax.plot(serie[var] * f, lw=lw, label=label)
                if obs is not None:
                    kge, alpha, beta, corr = KGEmod(obs[var], serie[var])
                    text = f'KGE={kge:.2f}  α={alpha:.2f}  β={beta:.2f}  ρ={corr:.2f}'
                    if var == 'outflow':
                        y = .97 - .08 * i
                        ha = 'top'
                    elif var == 'storage':
                        y = .03 + .08 * i
                        ha = 'bottom'
                    ax.text(0.01, y, text, ha='left', va=ha,
                            color=f'C{i}', transform=ax.transAxes, fontsize=10,
                            bbox=dict(facecolor='white', edgecolor='none', alpha=0.5))
            #for y in dct['thresholds']:
            #    ax.axhline(y * f, c='gray', lw=.5, ls=':', zorder=0)
            ax.set(title=var,
                   ylabel=dct['unit'],
                   xlim=(serie.index.min(), serie.index.max()))
            ax.spines[['top', 'right']].set_visible(False)
        
        fig.legend(*ax.get_legend_handles_labels(), loc=8, ncol=1 + len(sim), frameon=False)
        
        if save is not None:
            plt.savefig(save, dpi=300, bbox_inches='tight')