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
  
    
    
class bivariate_6pars(object):
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

    def __init__(self, inflow: pd.Series, storage: pd.Series, outflow: pd.Series, Vc: float, Vtot: float, Qmin: float, obj_func=kge):
        """
        Parameters:
        -----------
        inflow: pd.Series
            Inflow time seris used to force the model
        storage: pd.Series
            Observed storage time series that will be one of the targets in the calibration. The first value will be used as initial condition
        outflow: pd.Series
            Observed outflow time series that will be another target in the calibration
        Vc: float
            Volume (m3) associated to the conservative storage
        Vtot: float
            Total reservoir storage capacity (m3)
        Qmin: float
            Minimum outflow (m3/s)
        obj_func: 
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
            storage_init = self.storage.iloc[0]
        
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
        return [sim.outflow.round(2), sim.storage.round(1)]

    def evaluation(self, spinup: int = None):
        """It simply extracts the observed outflow from the class and removes (if necessary) the spinup time
        
        Inputs:
        -------
        spinup: int
            Numer or time steps to use to warm up the model. This initial time steps will not be taken into account in the computation of model performance
        """
        
        if spinup is not None:
            return [self.outflow.iloc[spinup:], self.storage.iloc[spinup:]]
        else:
            return [self.outflow, self.storage]

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
        of = []
        for i in range(len(evaluation)):
            perf = self.obj_func(evaluation[i], simulation[i])
            if isinstance(perf, tuple):
                of.append(1 - perf[0])
            else:
                of.append(1 - perf)
            
        return of
    
    
    
class bivariate_6pars_1of(object):
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

    def __init__(self, inflow: pd.Series, storage: pd.Series, outflow: pd.Series, Vc: float, Vtot: float, Qmin: float, obj_func=kge):
        """
        Parameters:
        -----------
        inflow: pd.Series
            Inflow time seris used to force the model
        storage: pd.Series
            Observed storage time series that will be one of the targets in the calibration. The first value will be used as initial condition
        outflow: pd.Series
            Observed outflow time series that will be another target in the calibration
        Vc: float
            Volume (m3) associated to the conservative storage
        Vtot: float
            Total reservoir storage capacity (m3)
        Qmin: float
            Minimum outflow (m3/s)
        obj_func: 
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
        
        # Just a way to keep this example flexible and applicable to various examples
        self.obj_func = obj_func       

    def simulation(self, pars: List[float], inflow: pd.Series = None, storage_init: float = None, spinup: int = None) -> List[pd.Series]:
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
            
        Outputs:
        --------
        sim: List[pd.Series]
            A list with two time series: outflow and storage
        """
        
        # forcings
        if inflow is None:
            inflow = self.inflow
        if storage_init is None:
            storage_init = self.storage.iloc[0]
        
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
        return [sim.outflow.round(2), sim.storage.round(1)]

    def evaluation(self, spinup: int = None):
        """It simply extracts the observed outflow from the class and removes (if necessary) the spinup time
        
        Inputs:
        -------
        spinup: int
            Numer or time steps to use to warm up the model. This initial time steps will not be taken into account in the computation of model performance
        """
        
        if spinup is not None:
            return [self.outflow.iloc[spinup:], self.storage.iloc[spinup:]]
        else:
            return [self.outflow, self.storage]

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
        of = []
        for i in range(len(evaluation)):
            perf = self.obj_func(evaluation[i], simulation[i])
            if isinstance(perf, tuple):
                of.append(1 - perf[0])
            else:
                of.append(1 - perf)
            
        # the objective function is the Euclidean distance from the origin (0, 0)
        return np.sqrt(np.sum(np.array(of)**2))
        # return of