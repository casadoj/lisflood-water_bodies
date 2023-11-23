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

       
        
class Hanazaki:
    """Representation of a reservoir according to Hanazaki, Yamazaki & Yoshimura (2021)."""
    
    def __init__(self, Vc: float, Vf:float, Ve: float, Vtot: float, Qn: float, Qf: float, A: int, At: int = 86400):
        """        
        Parameters:
        -----------
        Vc: float
            Volume (m3) associated to the conservative storage
        Vf: float
            Volume (m3) associated to the flood storage
        Ve: float
            Volume (m3) associated with an emergency situation
        Vtot: float
            Total reservoir storage capacity (m3)
        Qn: float
            Normal outflow (m3/s)
        Qf: float
            Outflow (m3/s) in case of flood
        A: float
            Area (m2) of the reservoir catchment
        At: int
            Simulation time step in seconds.
        """
        
        # storage limits
        self.Vc = Vc
        self.Vf = Vf
        self.Ve = Ve
        self.Vtot = Vtot
        
        # outflow limits
        self.Qn = Qn
        self.Qf = Qf
        
        # release coefficient
        self.k = max(1 - 5 * (Vtot - Vf) / A, 0)
        
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
        if V < self.Vc:
            Q = V * self.Qn / self.Vf
        elif V < self.Vf:
            if I < self.Qf:
                Q = self.Vc / self.Vf * self.Qn + ((V - self.Vc) / (self.Ve - self.Vc))**2 * (self.Qf - self.Vc / self.Vf * self.Qn)
            elif I >= self.Qf:
                Q = self.Vc / self.Vf * self.Qn + (V - self.Vc) / (self.Vf - self.Vc) * (self.Qf - self.Vc / self.Vf * self.Qn)
        elif V < self.Ve:
            if I < self.Qf:
                Q = self.Vc / self.Vf * self.Qn + ((V - self.Vc) / (self.Ve - self.Vc))**2 * (self.Qf - self.Vc / self.Vf * self.Qn)
            elif I >= self.Qf:
                Q = self.Qf + self.k * (V - self.Vf) / (self.Ve - self.Vf) * (I - self.Qf)
        elif self.Ve <= V:
            if I < self.Qf:
                Q = self.Qf
            elif I >= self.Qf:
                Q = I
            if verbose:
                if V > self.Vtot:
                    print(f'{V} m3 is greater than the reservoir capacity of {self.Vtot} m3')
            
        # 
        
        # update reservoir storage with the outflow volume
        AV = np.min([Q * self.At, V])
        AV = np.max([AV, V - self.Vtot])
        V -= AV
        
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
            Vo = self.Qn
        
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
            Whether to use the modified (default) of the orinigal Hanazaki's routine. The modified routine avoids the breaks in the outflow function at Vc, Vf and Ve.
            
        Returns:
        --------
        O: pd.Series
            Time series of reservoir outflow (m3/s)
        """

        if isinstance(I, float) or isinstance(I, int):
            assert I >= 0, '"I" must be a positive value'
            I = pd.Series(I, index=V.index)

        maskI = I < self.Qf
        maskV1 = V < self.Vc
        maskV2 = (self.Vc <= V) & (V < self.Vf)
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
                O[mask] = self.Vc / self.Vf * self.Qn + ((V[mask] - self.Vc) / (self.Ve - self.Vc))**2 * (self.Qf - self.Vc / self.Vf * self.Qn)
            else:
                O[mask] = .5 * self.Qn + ((V[mask] - self.Vc) / (self.Ve - self.Vc))**2 * (self.Qf - self.Qn)
        # ... and storage above emergency level
        O[maskI & maskV4] = self.Qf

        # if inflow over flood...
        # ... and storage in the normal level
        mask = ~maskI & maskV2
        if np.sum(mask) > 0:
            if modified:
                O[mask] = self.Vc / self.Vf * self.Qn + (V[mask] - self.Vc) / (self.Vf - self.Vc) * (self.Qf - self.Vc / self.Vf * self.Qn)
            else:
                O[mask] = .5 * self.Qn + (V[mask] - self.Vc) / (self.Vf - self.Vc) * (self.Qf - self.Qn)
        # ... and storage in flood level
        mask = ~maskI & maskV3
        if np.sum(mask) > 0:
            O[mask] = self.Qf + self.k * (V[mask] - self.Vf) / (self.Ve - self.Vf) * (I[mask] - self.Qf)
        # ... and storage in emergency level
        O[~maskI & maskV4] = I[~maskI & maskV4]

        return O
    
    def plot_routine(self, modified: bool = True, ax: Axes = None, **kwargs):
        """It creates a plot that explains the reservoir routine.
        
        Parameters:
        -----------
        modified: bool
            Whether to use the modified (default) of the orinigal Hanazaki's routine. The modified routine avoids the breaks in the outflow function at Vc, Vf and Ve.
        ax: Axes
            If provided, the plot will be added to the given axes
        """

        # dummy storage time series
        V = pd.Series(np.linspace(0, self.Vtot + .01, 1000))

        # create scatter plot
        if ax is None:
            fig, ax = plt.subplots(figsize=kwargs.get('figsize', (4, 4)))

        # outflow when inflow is lower than the flood outflow
        outflow1 = self.routine(V, .9 * self.Qf, modified=modified)
        ax.scatter(V, outflow1, s=.05, c='C0', label=r'$I < Q_f$')

        # outflow when inflow is larger than the flood outflow
        outflow2 = self.routine(V, 1.2 * self.Qf, modified=modified)
        ax.scatter(V, outflow2, s=.05, c='C1', label=r'$I \geq Q_f$')

        # reference storages and outflows
        vs = [self.Vc, self.Vf, self.Ve]
        if modified:
            qs = [self.Vc / self.Vf * self.Qn, self.Qf, self.Qf]
        else:
            qs = [.5 * self.Qn, self.Qf, self.Qf]
        
        for v, q in zip(vs, qs):
            ax.vlines(v, 0, q, color='k', ls=':', lw=.5, zorder=0)
            ax.hlines(q, 0, v, color='k', ls=':', lw=.5, zorder=0)
        
        # labels
        if modified:
            ax.text(0, qs[0], r'$\frac{V_c}{V_f} Q_n$', ha='left', va='bottom')
        else:
            ax.text(0, qs[0], r'$0.5 Q_n$', ha='left', va='bottom')
        ax.text(0, qs[1], r'$Q_f$', ha='left', va='bottom')
        ax.text(self.Vc, 0, r'$V_c$', rotation=90, ha='left', va='bottom')
        ax.text(self.Vf, 0, r'$V_f$', rotation=90, ha='right', va='bottom')
        ax.text(self.Ve, 0, r'$V_e$', rotation=90, ha='right', va='bottom')
        
        # setup
        ax.set(xlim=(0, self.Vtot),
               xlabel='storage (hm3)',
               ylim=(0, None),
               ylabel='outflow (m3/s)')
        ax.legend(frameon=False, loc=2)
        
    def normalize_timeseries(self, timeseries: pd.DataFrame) -> pd.DataFrame:
        """It normalizes the timeseries using the total reservoir capacity and the outflow associated to the flood limit. In this way, the storage time series ranges between 0 and 1, and the inflow and outflow time series are in the order of units.
        
        Parameters:
        -----------
        timeseries: pd.DataFrame
            A table with three columns ('storage', 'inflow', 'outflow') with the time series of a reservoir
            
        Returns:
        --------
        ts_norm: pd.DataFrame
            Table similar to the original but with normalized values
        """

        ts_norm = timeseries.copy()
        ts_norm.storage /= self.Vtot
        ts_norm[['inflow', 'outflow']] /= self.Qf

        return ts_norm
    
    def get_params(self):
        """It generates a dictionary with the reservoir paramenters in the Hanazaki model."""

        params = {'Vc': self.Vc,
                  'Vf': self.Vf,
                  'Ve': self.Ve,
                  'Vtot': self.Vtot,
                  'Qn': self.Qn,
                  'Qf': self.Qf,
                  'k': self.k}

        return params
    
    def compare(self, series1: pd.DataFrame, series2: pd.DataFrame = None, save: Union[Path, str] = None, **kwargs):
        """It compares two reservoir timeseries (inflow, outflow and storage) using the function 'reservoir_analysis'. If only 1 time series is given, the plot will simply show the reservoir behaviour of that set of time series.
        
        Inputs:
        -------
        series1: pd.DataFrame
            A table with the time series of 'inflow', 'outflow' and 'storage'
        series2: pd.DataFrame
            A second table with the time series of 'inflow', 'outflow' and 'storage'
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
        Vlims = np.array([self.Vc, self.Vf, self.Ve]) / self.Vtot
        
        # outflow limits
        Qlims = np.array([self.Vc / self.Vf * self.Qn, self.Qf, self.Qf]) / self.Qf
        
        # plot analysis
        series1_norm = self.normalize_timeseries(series1)
        if series2 is not None:
            series2_norm = self.normalize_timeseries(series2)
        else:
            series2_norm = series2
        reservoir_analysis(series1_norm, series2_norm,
                           x_thr=Vlims, y_thr=Qlims,
                           title=kwargs.get('title', None),
                           labels=kwargs.get('labels', ['sim', 'obs']),
                           alpha=kwargs.get('alpha', .05),
                           save=save)