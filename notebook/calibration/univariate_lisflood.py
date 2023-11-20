"""
Copyright 2023 by Jesús Casado Rodríguez
This file is part of Statistical Parameter Estimation Tool (SPOTPY).

:author: Jesús Casado Rodríguez

This example implements the python version of the reservoir routine in the hydrological model LISFLOOD-OS into SPOTPY.
"""

import os
import numpy as np
import pandas as pd
from spotpy.objectivefunctions import kge
from spotpy.parameter import Uniform
from typing import List
from reservoirs import Lisflood


class univariate_3pars(object):
    """This class allows for calibrating 3 parameters in the LISFLOOD reservoir routine.
    
    alpha: a value between 0 and 1 that defines the limit between the normal and flood zones
    beta: a factor that adjusts the initial estimate of the normal outflow
    k: release coefficient. A factor of the inflow that limits the outflow
    """
    
    alpha = Uniform(name='alpha', low=0.01, high=0.99)#, optguess=0.01)
    beta = Uniform(name='beta', low=0.25, high=2.0)#, optguess=1.59126)
    k = Uniform(name='k', low=1.0, high=5.0)

    def __init__(self, inflow: pd.Series, storage_init: float, outflow: pd.Series, Vc: float, Vn: float, Vf: float, Vtot: float, Qmin: float, Qn: float, Qnd: float, obj_func=kge):
        """
        Parameters:
        -----------
        inflow: pd.Series
            Inflow time seris used to force the model
        storage_init: float
            Initial reservoir storage
        outflow: pd.Series
            Observed outflow time series that will be the target in the calibration
        Vc: float
            storage (m3) associated to the conservative storage
        Vn: float
            Volume (m3) associated to the normal storage
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
        obj_func=None
        """
        
        # time series
        self.inflow = inflow
        self.outflow = outflow
        self.storage_init = storage_init
        
        # reservoir limits
        # volume
        self.Vc, self.Vn, self.Vf, self.Vtot = Vc, Vn, Vf, Vtot
        # outflow
        self.Qmin, self.Qn, self.Qnd = Qmin, Qn, Qnd
        
        # Just a way to keep this example flexible and applicable to various examples
        self.obj_func = obj_func       

    def simulation(self, pars: List[float], inflow: pd.Series = None, storage_init: float = None, spinup: int = None):
        """Given a parameter set, it declares the reservoir and runs the simulation.
        
        Inputs:
        -------
        pars: List
            A pair of values of the parameters 'alpha' and 'beta'
        inflow: pd.Series
            Inflow time seris used to force the model. If not given, the 'inflow' stored in the class will be used
        storage_init: float
            Initial reservoir storage. If not provided, the 'storage_init' stored in the class will be used
        spinup: int
            Numer or time steps to use to warm up the model. This initial time steps will not be taken into account in the computation of model performance
        """
        
        # forcings
        if inflow is None:
            inflow = self.inflow
        if storage_init is None:
            storage_init = self.storage_init
            
        # declare the reservoir with the effect of the parameters in 'x'
        Vn_adj = self.Vn + pars[0] * (self.Vf - self.Vn)
        Qn = max(self.Qmin, min(self.Qn * pars[1], self.Qnd))
        res = Lisflood(self.Vc, self.Vn, Vn_adj, self.Vf, self.Vtot, self.Qmin, Qn, self.Qnd)
        
        # simulate
        sim = res.simulate(inflow, storage_init, limit_Q=True, k=pars[2])
        if spinup is not None:
            sim = sim.iloc[spinup:]

        return sim.outflow.round(2)

    def evaluation(self, spinup: int = None):
        """It simply extracts the observed outflow from the class and removes (if necessary) the spinup time
        
        Inputs:
        -------
        spinup: int
            Numer or time steps to use to warm up the model. This initial time steps will not be taken into account in the computation of model performance
            """
        
        obs = self.outflow
        if spinup is not None:
            obs = obs.iloc[spinup:]
            
        return obs

    def objectivefunction(self, simulation, evaluation):
        """It computes the objective function defined in the class from the results of the simulation and the target series defined in the class
        
        Inputs:
        -------
        simulation: pd.Series
            Simulated time series
        evaluation: pd.Series
            Target time series
        """
        
        # compute the objective function
        of = 1 - self.obj_func(evaluation, simulation)
        
        return of
    
    
    
class univariate_6pars(object):
    """This class allows for calibrating 6 parameters in the LISFLOOD reservoir routine, 3 related to the storage limits, 2 to the outflow limits and the last one to the relation between inflow and outflow.
    
    FFn: fraction filled normal. The proportion of reservoir capacity that defines the lower limit of the normal storage zone
    FFf: fraction filled flood. The proportion of reservoir capacity that defines the upper limit of the flood zone
    alpha: a value between 0 and 1 that defines the limit between the normal and flood zones
    QQn: quantile outflow normal. The quantile of the inflow records that defines the normal outflow
    QQf: quantile outflow flood. The quantile of the inflow records that defines the flood outflow
    k: release coefficient. A factor of the inflow that limits the outflow
    """
    
    # FFn = Uniform(name='FFn', low=0.201, high=0.666)
    # FFf = Uniform(name='FFf', low=0.67, high=0.97)
    # alpha = Uniform(name='alpha', low=0.001, high=0.999)
    # QQn = Uniform(name='QQn', low=0.001, high=0.5)
    # QQf = Uniform(name='QQf', low=0.501, high=0.99)
    # k = Uniform(name='k', low=1.0, high=5.0)
    
    FFf = Uniform(name='FFf', low=0.20, high=0.99)
    alpha = Uniform(name='alpha', low=0.001, high=0.999)
    beta = Uniform(name='beta', low=0.001, high=0.999)
    QQf = Uniform(name='QQf', low=0.1, high=0.99)
    gamma = Uniform(name='gamma', low=0.001, high=0.999)
    k = Uniform(name='k', low=1.0, high=5.0)

    def __init__(self, inflow: pd.Series, storage: pd.Series, outflow: pd.Series, Vc: float, Vtot: float, Qmin: float, target: str = 'outflow', obj_func=kge):
        """
        Parameters:
        -----------
        inflow: pd.Series
            Inflow time seris used to force the model
        storage: pd.Series
            Time series of reservoir storage
        outflow: pd.Series
            Observed outflow time series that will be the target in the calibration
        Vc: float
            Volume (m3) associated to the conservative storage
        Vtot: float
            Total reservoir storage capacity (m3)
        Qmin: float
            Minimum outflow (m3/s)
        target: str
            Variable targeted in the calibration
        obj_func=None
        """
        
        # time series
        self.inflow = inflow
        self.outflow = outflow
        self.storage = storage
        
        # reservoir limits
        # volume
        self.Vc, self.Vtot = Vc, Vtot
        # outflow
        self.Qmin = Qmin
        
        # target variable and objective function
        self.target = target
        self.obj_func = obj_func       

    def simulation(self, pars: List[float], inflow: pd.Series = None, storage_init: float = None, spinup: int = None):
        """Given a parameter set, it declares the reservoir and runs the simulation.
        
        Inputs:
        -------
        pars: List
            A pair of values of the parameters 'alpha' and 'beta'
        inflow: pd.Series
            Inflow time seris used to force the model. If not given, the 'inflow' stored in the class will be used
        storage: float
            Initial reservoir storage. If not provided, the first value of the method 'storage' stored in the class will be used
        spinup: int
            Numer or time steps to use to warm up the model. This initial time steps will not be taken into account in the computation of model performance
        """
        
        # forcings
        if inflow is None:
            inflow = self.inflow
        if storage_init is None:
            storage_init = self.storage[0]
        
        # volume and outflow limits
        # Vn = self.Vtot * pars[0]
        # Vf = self.Vtot * pars[1]
        # Vn_adj = Vn + pars[2] * (Vf - Vn)
        # Qn = self.inflow.quantile(pars[3])
        # Qf = self.inflow.quantile(pars[4])
        
        # volume and outflow limits
        Vf = pars[0] * self.Vtot 
        Vn = self.Vc + pars[1] * (Vf - self.Vc)
        Vn_adj = Vn + pars[2] * (Vf - Vn)
        Qf = self.inflow.quantile(pars[3])
        Qn = pars[4] * Qf
        
        # declare the reservoir with the effect of the parameters in 'x'
        res = Lisflood(self.Vc, Vn, Vn_adj, Vf, self.Vtot, self.Qmin, Qn, Qf)
        
        # simulate
        sim = res.simulate(inflow, storage_init, limit_Q=True, k=pars[5])
        if spinup is not None:
            sim = sim.iloc[spinup:]
        
        self.reservoir = res
        return sim[self.target].round(2)

    def evaluation(self, spinup: int = None):
        """It simply extracts the observed outflow from the class and removes (if necessary) the spinup time
        
        Inputs:
        -------
        spinup: int
            Numer or time steps to use to warm up the model. This initial time steps will not be taken into account in the computation of model performance
            """
        
        if self.target == 'outflow':
            obs = self.outflow
        elif self.target == 'storage':
            obs = self.storage
        if spinup is not None:
            obs = obs.iloc[spinup:]
            
        return obs

    def objectivefunction(self, simulation, evaluation):
        """It computes the objective function defined in the class from the results of the simulation and the target series defined in the class
        
        Inputs:
        -------
        simulation: pd.Series
            Simulated time series
        evaluation: pd.Series
            Target time series
        """
        
        # compute the objective function       
        of = self.obj_func(evaluation, simulation)
        if isinstance(of, tuple):
            of = 1 - of[0]
        else:
            of = 1 - of
        
        return of