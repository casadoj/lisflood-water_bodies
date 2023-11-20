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


class Lisflood:
    """Representation of a reservoir in the LISFLOOD-OS hydrological model."""
    
    def __init__(self, Vc: float, Vn: float, Vn_adj: float, Vf: float, Vtot: float, Qmin: float, Qn: float, Qnd: float, At: int = 86400):
        """
        Parameters:
        -----------
        Vc: float
            Volume (m3) associated to the conservative storage
        Vn: float
            Volume (m3) associated to the normal storage
        Vn_adj: float
            Volume (m3) associated to the adjusted (calibrated) normal storage
        Vf: float
            Volume (m3) associated to the flood storage
        Vtot: float
            Total reservoir storage capacity (m3)
        Qmin: float
            Minimum outflow (m3/s)
        Qn: float
            Normal outflow (m3/s)
        Qnd: float
            Non-damaging outflow (m3/s)
        At: int
            Simulation time step in seconds.
        """
        
        # storage limits
        self.Vc = Vc
        self.Vn = Vn
        self.Vn_adj = Vn_adj
        self.Vf = Vf
        self.Vtot = Vtot
        
        # outflow limits
        self.Qmin = Qmin
        self.Qn = Qn
        self.Qnd = Qnd
        
        # time step duration in seconds
        self.At = At
    
    def timestep(self, I: float, V: float, limit_Q: bool = True, type: int = 1, k: float  = 1.2, verbose: bool = False) -> List[float]:
        """Given an inflow and an initial storage values, it computes the corresponding outflow
        
        Parameters:
        -----------
        I: float
            Inflow (m3/s)
        V: float
            Volume stored in the reservoir (m3)
        limit_Q: bool
            Whether to limit the outflow in the flood zone when it exceeds inflow by more than 1.2 times
        k: float
            Release coefficient. If the reservoir is in the flood zone, the outflow is limited to k times the inflow
        verbose: bool
            Whether to show on screen the evolution
            
        Returns:
        --------
        Q, V: List[float]
            Outflow (m3/s) and updated storage (m3)
        """
        
        # update reservoir storage with the inflow volume
        V += I * self.At
        
        # ouflow depending on the storage level
        if V < 2 * self.Vc:
            Q = np.min((self.Qmin, V / self.At))
        elif V < self.Vn:
            Q = self.Qmin + (self.Qn - self.Qmin) * (V - 2 * self.Vc) / (self.Vn - 2 * self.Vc)
        elif V < self.Vn_adj:
            Q = self.Qn
        elif V < self.Vf:
            Q = self.Qn + (self.Qnd - self.Qn) * (V - self.Vn_adj) / (self.Vf - self.Vn_adj)
        elif V > self.Vf:
            Q = np.max([(V - self.Vf - .01 * self.Vtot) / self.At, np.min([self.Qnd, np.max([1.2 * I, self.Qn])])])
            if verbose:
                if V > self.Vtot:
                    print(f'{V} m3 is greater than the reservoir capacity of {self.Vtot} m3')
        
        # limit outflow depending on the inflow
        if limit_Q:
            if (Q > 1.2 * I) & (Q > self.Qn) & (V < self.Vf): # isn't this the same as: (Q > 1.2 * I) & (V > self.Vn_adj) & (V < self.Vf)
                if type == 1:
                    Q = np.min([Q, np.max([k * I, self.Qn])])
                elif type == 2:
                    Q = np.min([self.Qnd, np.max([k * I, self.Qn])])

        # update reservoir storage with the outflow volume
        AV = np.min([Q * self.At, V])
        AV = np.max([AV, V - self.Vtot])
        V -= AV
        
        return Q, V

    def simulate(self, inflow: pd.Series, Vo: float = None, limit_Q: bool = True, k: float = 1):
        """Given a inflow time series (m3/s) and an initial storage (m3), it computes the time series of outflow (m3/s) and storage (m3)
        
        Parameters:
        -----------
        inflow: pd.Series
            Time series of flow coming into the reservoir (m3/s)
        Vo: float
            Initial value of reservoir storage (m3). If not provided, it is assumed that the normal storage is the initial condition
        limit_Q: bool
            Whether to limit the outflow in the flood zone when it exceeds inflow by more than 1.2 times
        k: float
            Release coefficient. If the reservoir is in the flood zone, the outflow is limited to k times the inflow
            
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
            Q, V = self.timestep(inflow[ts], Vo, limit_Q=limit_Q, k=k)
            storage[ts] = V
            outflow[ts] = Q
            # update current storage
            Vo = V

        return pd.concat((storage, inflow, outflow), axis=1)
        
        
    def routine(self, V: pd.Series, I: Union[float, pd.Series]):
        """Given a time series of reservoir storage (m3) and a value or a time series of inflow (m3/s), it computes the ouflow (m3/s). This function is only meant for explanatory purposes; since the volume time series is given, the computed outflow does not update the reservoir storage. If the intention is to simulate the behaviour of the reservoir, refer to the function "simulate"
        
        Parameters:
        -----------
        V: pd.Series
            Time series of reservoir storage (m3)
        I: Union[float, pd.Series]
            Reservor inflow (m3/s)
            
        Returns:
        --------
        O: pd.Series
            Time series of reservoir outflow (m3/s)
        """
        
        if isinstance(I, float) or isinstance(I, int):
            assert I >= 0, '"I" must be a positive value'
            I = pd.Series(I, index=V.index)
        
        O1 = V / self.At 
        O1[O1 > self.Qmin] = self.Qmin
        O = O1.copy()
        
        O2 = self.Qmin + (self.Qn - self.Qmin) * (V - 2 * self.Vc) / (self.Vn - 2 * self.Vc)
        maskV2 = (2 * self.Vc <= V) & (V < self.Vn)
        O[maskV2] = O2[maskV2]
        
        O3 = pd.Series(self.Qn, index=V.index)
        maskV3 = (self.Vn <= V) & (V < self.Vn_adj)
        O[maskV3] = O3[maskV3]
        
        O4 = self.Qn + (self.Qnd - self.Qn) * (V - self.Vn_adj) / (self.Vf - self.Vn_adj)
        maskV4 = (self.Vn_adj <= V) & (V < self.Vf)
        O[maskV4] = O4[maskV4]
        
        Omax = 1.2 * I
        Omax[Omax < self.Qn] = self.Qn
        Omax[Omax > self.Qnd] = self.Qnd
        O5 = pd.concat(((V - self.Vf - .01 * self.Vtot) / self.At, Omax), axis=1).max(axis=1)
        maskV5 = self.Vf <= V
        O[maskV5] = O5[maskV5]
        
        Oreg = I
        Oreg[Oreg < self.Qn] = self.Qn
        Oreg = pd.concat((O, Oreg), axis=1).min(axis=1)
        maskO = (O > 1.2 * I) & (O > self.Qn) & (V < self.Vf)
        O[maskO] = Oreg[maskO]
        
        temp = pd.concat((O1, O2, O3, O4, O5, Omax, Oreg), axis=1)
        temp.columns = ['O1', 'O2', 'O3', 'O4', 'O5', 'Omax', 'Oreg']
        self.O = temp
        
        return O
    
    def plot_routine(self, ax: Axes = None, **kwargs):
        """It creates a plot that explains the reservoir routine.
        
        Parameters:
        -----------
        ax: Axes
            If provided, the plot will be added to the given axes
        """

        # dummy storage time series
        V = pd.Series(np.linspace(0, self.Vtot + .01, 1000))

        # create scatter plot
        if ax is None:
            fig, ax = plt.subplots(figsize=kwargs.get('figsize', (5, 5)))

        # outflow
        outflow = self.routine(V, I=self.Qnd)
        ax.plot(V, outflow, lw=1, c='C0')

        # reference storages and outflows
        vs = [self.Vc, 2 * self.Vc, self.Vn, self.Vn_adj, self.Vf]
        qs = [self.Qmin, self.Qmin, self.Qn, self.Qn, self.Qnd]
        for v, q in zip(vs, qs):
            ax.vlines(v, 0, q, color='k', ls=':', lw=.5, zorder=0)
            ax.hlines(q, 0, v, color='k', ls=':', lw=.5, zorder=0)
        
        # labels
        ax.text(0, self.Qmin, r'$Q_{min}$', ha='left', va='bottom')
        ax.text(0, self.Qn, r'$Q_{n,adj}$', ha='left', va='bottom')
        ax.text(0, self.Qnd, r'$Q_{nd}$', ha='left', va='bottom')
        ax.text(self.Vn, 0, r'$V_n$', rotation=90, ha='right', va='bottom')
        ax.text(self.Vn_adj, 0, r'$V_{n,adj}$', rotation=90, ha='right', va='bottom')
        ax.text(self.Vf, 0, r'$V_f$', rotation=90, ha='right', va='bottom')
        
        # setup
        ax.set(xlim=(0, self.Vtot),
               xlabel='storage (hm3)',
               ylim=(0, None),
               ylabel='outflow (m3/s)')
        ax.set_title('LISFLOOD reservoir routine')
        
    def normalize_timeseries(self, timeseries: pd.DataFrame) -> pd.DataFrame:
        """It normalizes the timeseries using the total reservoir capacity and the non-damaging outflow. In this way, the storage time series ranges between 0 and 1, and the inflow and outflow time series are in the order of units.
        
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
        ts_norm[['inflow', 'outflow']] /= self.Qnd

        return ts_norm
        

        
class Lisflood:
    """Representation of a reservoir in the LISFLOOD-OS hydrological model."""
    
    def __init__(self, Vc: float, Vn: float, Vn_adj: float, Vf: float, Vtot: float, Qmin: float, Qn: float, Qnd: float, At: int = 86400):
        """
        Parameters:
        -----------
        Vc: float
            Volume (m3) associated to the conservative storage
        Vn: float
            Volume (m3) associated to the normal storage
        Vn_adj: float
            Volume (m3) associated to the adjusted (calibrated) normal storage
        Vf: float
            Volume (m3) associated to the flood storage
        Vtot: float
            Total reservoir storage capacity (m3)
        Qmin: float
            Minimum outflow (m3/s)
        Qn: float
            Normal outflow (m3/s)
        Qnd: float
            Non-damaging outflow (m3/s)
        At: int
            Simulation time step in seconds.
        """
        
        # storage limits
        self.Vc = Vc
        self.Vn = Vn
        self.Vn_adj = Vn_adj
        self.Vf = Vf
        self.Vtot = Vtot
        
        # outflow limits
        self.Qmin = Qmin
        self.Qn = Qn
        self.Qnd = Qnd
        
        # time step duration in seconds
        self.At = At
    
    def timestep(self, I: float, V: float, limit_Q: bool = True, k: float  = 1.2) -> List[float]:
        """Given an inflow and an initial storage values, it computes the corresponding outflow
        
        Parameters:
        -----------
        I: float
            Inflow (m3/s)
        V: float
            Volume stored in the reservoir (m3)
        limit_Q: bool
            Whether to limit the outflow in the flood zone when it exceeds inflow by more than 'k' times
        k: float
            Release coefficient. If the reservoir is in the flood zone, the outflow is limited to k times the inflow
            
        Returns:
        --------
        Q, V: List[float]
            Outflow (m3/s) and updated storage (m3)
        """
        
        # update reservoir storage with the inflow volume
        V += I * self.At
        
        # ouflow depending on the storage level
        if V < 2 * self.Vc:
            Q = self.Qmin
        elif V < self.Vn:
            Q = self.Qmin + (self.Qn - self.Qmin) * (V - 2 * self.Vc) / (self.Vn - 2 * self.Vc)
        elif V < self.Vn_adj:
            Q = self.Qn
        elif V < self.Vf:
            Q = self.Qn + (self.Qnd - self.Qn) * (V - self.Vn_adj) / (self.Vf - self.Vn_adj)
            if limit_Q:
                if Q > k * I:
                    Q = np.max([k * I, self.Qn])
                    # # Q <= Qnd at this storage zone, so this second approach (from the documentation) makes no sense
                    # Q = np.min([self.Qnd, np.max([k * I, self.Qn])])
        elif V > self.Vf:
            Q = np.max([(V - self.Vf) / self.At, np.min([self.Qnd, np.max([k * I, self.Qn])])])
        
        # limit outflow so the final storage is between 0 and 1
        Q = np.max([np.min([Q, (V - self.Vc) / self.At]), (V - self.Vtot) / self.At])

        # update reservoir storage with the outflow volume
        V -= Q * self.At
        
        assert 0 <= V, 'The volume at the end of the timestep is negative.'
        assert V <= self.Vtot, 'The volume at the end of the timestep is larger than the total reservoir capacity.'
        
        return Q, V
    
    def timestep2(self, I: float, V: float, k: float = 1) -> List[float]:
        """Given an inflow and an initial storage values, it computes the corresponding outflow
        
        Parameters:
        -----------
        I: float
            Inflow (m3/s)
        V: float
            Volume stored in the reservoir (m3)
        limit_Q: bool
            Whether to limit the outflow in the flood zone when it exceeds inflow by more than 1.2 times
        k: float
            Release coefficient. If the reservoir is in the flood zone, the outflow is limited to k times the inflow
        verbose: bool
            Whether to show on screen the evolution
            
        Returns:
        --------
        Q, V: List[float]
            Outflow (m3/s) and updated storage (m3)
        """
        
        # update reservoir storage with the inflow volume
        V += I * self.At
        
        # ouflow depending on the storage level
        if V < 2 * self.Vc:
            Q = self.Qmin
        elif V < self.Vn:
            Q = self.Qmin + (self.Qn - self.Qmin) * (V - 2 * self.Vc) / (self.Vn - 2 * self.Vc)
        elif V < self.Vn_adj:
            Q = self.Qn
        elif V < self.Vf:
            Q = self.Qn + (self.Qnd - self.Qn) * (V - self.Vn_adj) / (self.Vf - self.Vn_adj)
        elif V > self.Vf:
            Q = np.min([(V - self.Vf) / self.At, np.max([self.Qnd, k * I])])
            
        # limit outflow so the final storage is between 0 and 1
        Q = np.max([np.min([Q, (V - self.Vc) / self.At]), (V - self.Vtot) / self.At])

        # update reservoir storage with the outflow volume
        # AV = np.min([Q * self.At, V])
        # AV = np.max([AV, V - self.Vtot])
        V -= Q * self.At
        
        assert 0 <= V, 'The volume at the end of the timestep is negative.'
        assert V <= self.Vtot, 'The volume at the end of the timestep is larger than the total reservoir capacity.'
        
        return Q, V
    
    def timestep3(self, I: float, V: float, limit_Q: bool = True, k: float  = 1.2) -> List[float]:
        """Given an inflow and an initial storage values, it computes the corresponding outflow
        
        Parameters:
        -----------
        I: float
            Inflow (m3/s)
        V: float
            Volume stored in the reservoir (m3)
        limit_Q: bool
            Whether to limit the outflow in the flood zone when it exceeds inflow by more than 'k' times
        k: float
            Release coefficient. If the reservoir is in the flood zone, the outflow is limited to k times the inflow
            
        Returns:
        --------
        Q, V: List[float]
            Outflow (m3/s) and updated storage (m3)
        """
        
        # update reservoir storage with the inflow volume
        V += I * self.At
        
        # ouflow depending on the storage level
        if V < 2 * self.Vc:
            Q = self.Qmin
        elif V < self.Vn:
            Q = self.Qmin + (self.Qn - self.Qmin) * (V - 2 * self.Vc) / (self.Vn - 2 * self.Vc)
        elif V < self.Vn_adj:
            Q = self.Qn
        elif V < self.Vf:
            Q = self.Qn + (self.Qnd - self.Qn) * (V - self.Vn_adj) / (self.Vf - self.Vn_adj)
            if limit_Q:
                if Q > k * I:
                    Q = np.max([k * I, self.Qn])
                    # # Q <= Qnd at this storage zone, so this second approach (from the documentation) makes no sense
                    # Q = np.min([self.Qnd, np.max([k * I, self.Qn])])
        elif V > self.Vf:
            Q = np.min([self.Qnd, np.max([k * I, self.Qn])])
        
        # limit outflow so the final storage is between 0 and 1
        Q = np.max([np.min([Q, (V - self.Vc) / self.At]), (V - self.Vtot) / self.At])

        # update reservoir storage with the outflow volume
        V -= Q * self.At
        
        assert 0 <= V, 'The volume at the end of the timestep is negative.'
        assert V <= self.Vtot, 'The volume at the end of the timestep is larger than the total reservoir capacity.'
        
        return Q, V
    
    def timestep4(self, I: float, V: float, limit_Q: bool = True, k: float  = 1.2, p: float = 3.333) -> List[float]:
        """Given an inflow and an initial storage values, it computes the corresponding outflow
        
        Parameters:
        -----------
        I: float
            Inflow (m3/s)
        V: float
            Volume stored in the reservoir (m3)
        limit_Q: bool
            Whether to limit the outflow in the flood zone when it exceeds inflow by more than 'k' times
        k: float
            Release coefficient. If the reservoir is in the flood zone, the outflow is limited to k times the inflow
            
        Returns:
        --------
        Q, V: List[float]
            Outflow (m3/s) and updated storage (m3)
        """
        
        # update reservoir storage with the inflow volume
        V += I * self.At
        
        # ouflow depending on the storage level
        if V < 2 * self.Vc:
            Q = self.Qmin
        elif V < self.Vn:
            Q = self.Qmin + (self.Qn - self.Qmin) * (V - 2 * self.Vc) / (self.Vn - 2 * self.Vc)
        elif V < self.Vn_adj:
            Q = self.Qn
        elif V < self.Vf:
            Q = self.Qn + (self.Qnd - self.Qn) * (V - self.Vn_adj) / (self.Vf - self.Vn_adj)
            if limit_Q:
                if Q > k * I:
                    Q = np.max([k * I, self.Qn])
                    # # Q <= Qnd at this storage zone, so this second approach (from the documentation) makes no sense
                    # Q = np.min([self.Qnd, np.max([k * I, self.Qn])])
        elif V > self.Vf:
            Q = np.max([np.min([(V - self.Vf) / self.At, p * self.Qnd]), np.min([self.Qnd, np.max([k * I, self.Qn])])])
        
        # limit outflow so the final storage is between 0 and 1
        Q = np.max([np.min([Q, (V - self.Vc) / self.At]), (V - self.Vtot) / self.At])

        # update reservoir storage with the outflow volume
        V -= Q * self.At
        
        assert 0 <= V, 'The volume at the end of the timestep is negative.'
        assert V <= self.Vtot, 'The volume at the end of the timestep is larger than the total reservoir capacity.'
        
        return Q, V
    
    def timestep5(self, I: float, V: float, limit_Q: bool = True, k: float = 1.2, tol: float = 1e-6) -> List[float]:
        """Given an inflow and an initial storage values, it computes the corresponding outflow
        
        Parameters:
        -----------
        I: float
            Inflow (m3/s)
        V: float
            Volume stored in the reservoir (m3)
        limit_Q: bool
            Whether to limit the outflow in the flood zone when it exceeds inflow by more than 1.2 times
        k: float
            Release coefficient. If the reservoir is in the flood zone, the outflow is limited to k times the inflow
            
        Returns:
        --------
        Q, V: List[float]
            Outflow (m3/s) and updated storage (m3)
        """
        
        # update reservoir storage with the inflow volume
        V += I * self.At
        
        # ouflow depending on the storage level
        if V < 2 * self.Vc:
            Q = np.min([self.Qmin, (V - self.Vc) / self.At])
        elif V < self.Vn:
            Q = self.Qmin + (self.Qn - self.Qmin) * (V - 2 * self.Vc) / (self.Vn - 2 * self.Vc)
        elif V < self.Vn_adj:
            Q = self.Qn
        elif V < self.Vf:
            Q = self.Qn + (self.Qnd - self.Qn) * (V - self.Vn_adj) / (self.Vf - self.Vn_adj)
            if limit_Q:
                if Q > k * I:
                    Q = np.max([k * I, self.Qn])
        elif V > self.Vf:
            # Q = np.max([(V - self.Vf - tol * self.Vtot) / self.At, np.min([self.Qnd, np.max([1.2 * I, self.Qn])])])
            Q = np.max([self.Qnd, I])
            
        # limit outflow so the final storage is between 0 and 1
        Q = np.max([np.min([Q, (V - self.Vc) / self.At]), (V - self.Vtot) / self.At])
        # Q = np.max([np.min([Q, (V - self.Vc) / self.At]), I])

        # update reservoir storage with the outflow volume
        # AV = np.min([Q * self.At, V])
        # AV = np.max([AV, V - self.Vtot])
        V -= Q * self.At
        
        assert 0 <= V, 'The volume at the end of the timestep is negative.'
        assert V <= self.Vtot, 'The volume at the end of the timestep is larger than the total reservoir capacity.'
        
        return Q, V
    
    def timestep6(self, I: float, Io: float, V: float, limit_Q: bool = True, k: float = 1.2) -> List[float]:
        """Given an inflow and an initial storage values, it computes the corresponding outflow
        
        Parameters:
        -----------
        I: float
            Inflow (m3/s)
        V: float
            Volume stored in the reservoir (m3)
        limit_Q: bool
            Whether to limit the outflow in the flood zone when it exceeds inflow by more than 1.2 times
        k: float
            Release coefficient. If the reservoir is in the flood zone, the outflow is limited to k times the inflow
            
        Returns:
        --------
        Q, V: List[float]
            Outflow (m3/s) and updated storage (m3)
        """
        
        # update reservoir storage with the inflow volume
        V += I * self.At
        
        # ouflow depending on the storage level
        if V < 2 * self.Vc:
            Q = self.Qmin
        elif V < self.Vn:
            Q = self.Qmin + (self.Qn - self.Qmin) * (V - 2 * self.Vc) / (self.Vn - 2 * self.Vc)
        elif V < self.Vn_adj:
            Q = self.Qn
        elif V < self.Vf:
            Q = self.Qn + (self.Qnd - self.Qn) * (V - self.Vn_adj) / (self.Vf - self.Vn_adj)
            if limit_Q:
                if Q > k * I:
                    Q = np.max([k * I, self.Qn])
        elif V > self.Vf:
            if I > Io:
                Q = np.max([self.Qnd, I / k])
            else:
                Q = np.max([self.Qnd, k * I])
            
        # limit outflow so the final storage is between 0 and 1
        Q = np.max([np.min([Q, (V - self.Vc) / self.At]), (V - self.Vtot) / self.At])

        # update reservoir storage with the outflow volume
        V -= Q * self.At
        
        assert 0 <= V, 'The volume at the end of the timestep is negative.'
        assert V <= self.Vtot, 'The volume at the end of the timestep is larger than the total reservoir capacity.'
        
        return Q, V
    
    def simulate(self, inflow: pd.Series, Vo: float = None, limit_Q: bool = True, routine: int = 1, k: float = 1):
        """Given a inflow time series (m3/s) and an initial storage (m3), it computes the time series of outflow (m3/s) and storage (m3)
        
        Parameters:
        -----------
        inflow: pd.Series
            Time series of flow coming into the reservoir (m3/s)
        Vo: float
            Initial value of reservoir storage (m3). If not provided, it is assumed that the normal storage is the initial condition
        limit_Q: bool
            Whether to limit the outflow in the flood zone when it exceeds inflow by more than 1.2 times
        k: float
            Release coefficient. If the reservoir is in the flood zone, the outflow is limited to k times the inflow
            
        Returns:
        --------
        pd.DataFrame
            A table that concatenates the storage, inflow and outflow time series.
        """
        
        if Vo is None:
            Vo = self.Qn
            
        routines = {1: self.timestep,
                    2: self.timestep2,
                    3: self.timestep3,
                    4: self.timestep4,
                    5: self.timestep5,
                    6: self.timestep6}
        
        storage = pd.Series(index=inflow.index, dtype=float, name='storage')
        outflow = pd.Series(index=inflow.index, dtype=float, name='outflow')
        # for ts in tqdm(inflow.index):
        for ts in inflow.index:
            try:
                # compute outflow and new storage
                if routine == 2:
                    Q, V = routines[routine](inflow[ts], Vo, k=k)
                elif routine == 6:
                    try:
                        Q, V = routines[routine](inflow[ts], inflow[ts - timedelta(seconds=self.At)], Vo, limit_Q=limit_Q, k=k)
                    except:
                        Q, V = routines[routine](inflow[ts], inflow[ts], Vo, limit_Q=limit_Q, k=k)
                else:
                    Q, V = routines[routine](inflow[ts], Vo, limit_Q=limit_Q, k=k)
            except:
                print(ts)
                return pd.concat((storage, inflow, outflow), axis=1)
            storage[ts] = V
            outflow[ts] = Q
            # update current storage
            Vo = V

        return pd.concat((storage, inflow, outflow), axis=1)
        
    def routine(self, V: pd.Series, I: Union[float, pd.Series]):
        """Given a time series of reservoir storage (m3) and a value or a time series of inflow (m3/s), it computes the ouflow (m3/s). This function is only meant for explanatory purposes; since the volume time series is given, the computed outflow does not update the reservoir storage. If the intention is to simulate the behaviour of the reservoir, refer to the function "simulate"
        
        Parameters:
        -----------
        V: pd.Series
            Time series of reservoir storage (m3)
        I: Union[float, pd.Series]
            Reservor inflow (m3/s)
            
        Returns:
        --------
        O: pd.Series
            Time series of reservoir outflow (m3/s)
        """
        
        if isinstance(I, float) or isinstance(I, int):
            assert I >= 0, '"I" must be a positive value'
            I = pd.Series(I, index=V.index)
        
        O1 = V / self.At 
        O1[O1 > self.Qmin] = self.Qmin
        O = O1.copy()
        
        O2 = self.Qmin + (self.Qn - self.Qmin) * (V - 2 * self.Vc) / (self.Vn - 2 * self.Vc)
        maskV2 = (2 * self.Vc <= V) & (V < self.Vn)
        O[maskV2] = O2[maskV2]
        
        O3 = pd.Series(self.Qn, index=V.index)
        maskV3 = (self.Vn <= V) & (V < self.Vn_adj)
        O[maskV3] = O3[maskV3]
        
        O4 = self.Qn + (self.Qnd - self.Qn) * (V - self.Vn_adj) / (self.Vf - self.Vn_adj)
        maskV4 = (self.Vn_adj <= V) & (V < self.Vf)
        O[maskV4] = O4[maskV4]
        
        Omax = 1.2 * I
        Omax[Omax < self.Qn] = self.Qn
        Omax[Omax > self.Qnd] = self.Qnd
        O5 = pd.concat(((V - self.Vf - .01 * self.Vtot) / self.At, Omax), axis=1).max(axis=1)
        maskV5 = self.Vf <= V
        O[maskV5] = O5[maskV5]
        
        Oreg = I
        Oreg[Oreg < self.Qn] = self.Qn
        Oreg = pd.concat((O, Oreg), axis=1).min(axis=1)
        maskO = (O > 1.2 * I) & (O > self.Qn) & (V < self.Vf)
        O[maskO] = Oreg[maskO]
        
        temp = pd.concat((O1, O2, O3, O4, O5, Omax, Oreg), axis=1)
        temp.columns = ['O1', 'O2', 'O3', 'O4', 'O5', 'Omax', 'Oreg']
        self.O = temp
        
        return O
       
    def plot_routine(self, ax: Axes = None, **kwargs):
        """It creates a plot that explains the reservoir routine.
        
        Parameters:
        -----------
        ax: Axes
            If provided, the plot will be added to the given axes
        """

        # dummy storage time series
        V = pd.Series(np.linspace(0, self.Vtot + .01, 1000))

        # create scatter plot
        if ax is None:
            fig, ax = plt.subplots(figsize=kwargs.get('figsize', (5, 5)))

        # outflow
        outflow = self.routine(V, I=self.Qnd)
        ax.plot(V, outflow, lw=1, c='C0')

        # reference storages and outflows
        vs = [self.Vc, 2 * self.Vc, self.Vn, self.Vn_adj, self.Vf]
        qs = [self.Qmin, self.Qmin, self.Qn, self.Qn, self.Qnd]
        for v, q in zip(vs, qs):
            ax.vlines(v, 0, q, color='k', ls=':', lw=.5, zorder=0)
            ax.hlines(q, 0, v, color='k', ls=':', lw=.5, zorder=0)
        
        # labels
        ax.text(0, self.Qmin, r'$Q_{min}$', ha='left', va='bottom')
        ax.text(0, self.Qn, r'$Q_{n,adj}$', ha='left', va='bottom')
        ax.text(0, self.Qnd, r'$Q_nd$', ha='left', va='bottom')
        ax.text(self.Vn, 0, r'$V_n$', rotation=90, ha='right', va='bottom')
        ax.text(self.Vn_adj, 0, r'$V_{n,adj}$', rotation=90, ha='right', va='bottom')
        ax.text(self.Vf, 0, r'$V_f$', rotation=90, ha='right', va='bottom')
        
        # setup
        ax.set(xlim=(0, self.Vtot),
               xlabel='storage (hm3)',
               ylim=(0, None),
               ylabel='outflow (m3/s)')
        ax.set_title('LISFLOOD reservoir routine')
        
    def normalize_timeseries(self, timeseries: pd.DataFrame) -> pd.DataFrame:
        """It normalizes the timeseries using the total reservoir capacity and the non-damaging outflow. In this way, the storage time series ranges between 0 and 1, and the inflow and outflow time series are in the order of units.
        
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
        ts_norm[['inflow', 'outflow']] /= self.Qnd

        return ts_norm
    
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
        Vlims = np.array([2 * self.Vc, self.Vn, self.Vn_adj, self.Vf])
        #Vlims = np.array([2 * self.Vc, self.Vn, self.Vn_adj, self.Vf]) / self.Vtot
        
        # outflow limits
        Qlims = np.array([self.Qmin, self.Qn, self.Qn, self.Qnd])
        #Qlims = np.array([self.Qmin, self.Qn, self.Qn, self.Qnd]) / self.Qnd
        
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
                           x_thr=Vlims, y_thr=Qlims,
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
                                 'thresholds': [self.Qmin, self.Qn, self.Qnd]},
                     'storage': {'unit': 'hm3',
                                 'factor': 1e-6,
                                 'thresholds': [self.Vc, self.Vn, self.Vn_adj, self.Vf, self.Vtot]}}

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
            for y in dct['thresholds']:
                ax.axhline(y * f, c='gray', lw=.5, ls=':', zorder=0)
            ax.set(title=var,
                   ylabel=dct['unit'],
                   xlim=(serie.index.min(), serie.index.max()))
            ax.spines[['top', 'right']].set_visible(False)
        
        fig.legend(*ax.get_legend_handles_labels(), loc=8, ncol=1 + len(sim), frameon=False)
        
        if save is not None:
            plt.savefig(save, dpi=300, bbox_inches='tight')
        
    def get_params(self) -> Dict:
        """It generates a dictionary with the reservoir parameters
        
        Returns:
        --------
        params: Dict
            A dictionary with the name and value of the reservoir parameters
        """

        params = {'Vc': self.Vc,
                  'Vn': self.Vn,
                  'Vn_adj': self.Vn_adj,
                  'Vf': self.Vf,
                  'Vtot': self.Vtot,
                  'Qmin': self.Qmin,
                  'Qn': self.Qn,
                  'Qnd': self.Qnd}

        return params
        
        
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