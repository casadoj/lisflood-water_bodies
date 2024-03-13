import numpy as np
import pandas as pd
import logging
import matplotlib.pyplot as plt
from typing import Union, List, Dict

# crear logger
FRND_logger = logging.getLogger('metricas_rendimiento')

def NSE(obs: pd.Series, sim: pd.Series) -> float:
    """It computes the Nash-Sutcliffe efficiency coefficient.
    
    Parameters:
    -----------
    obs:   pd.Series.
        Observed time series
    sim:   pd.Series
        Simulated time series
    
    Returns:
    -------
    NSE: float
        Nash-Sutcliffe efficiency coefficient
    """
    
    # Eliminar pasos sin dato
    data = pd.concat((obs, sim), axis=1)
    data.columns = ['obs', 'sim']
    data.dropna(axis=0, how='any', inplace=True)
    # Para la función si no hay datos
    assert data.shape[0] > 0, "ERROR. No indices match between the observed and the simulated series."
    # Calcular NSE
    nse = 1 - sum((data.obs - data.sim)**2) / sum((data.obs - data.obs.mean())**2)
    
    return nse


def RMSE(serie1, serie2):
    """Calcula la raíz del error medio cuadrático.
    
    Parámetros:
    -----------
    serie1:    series. Primera del par de series a comparar
    serie2:    series. Segunda del par de series a comparar"""
    
    data = pd.concat((serie1, serie2), axis=1)
    data.columns = ['obs', 'sim']
    # Eliminar pasos sin dato
    data.dropna(axis=0, how='any', inplace=True)
    # Para la función si no hay datos
    if data.shape[0] == 0:
        FRND_logger.error('No hay valores')
        return 
    # Calcular RMSE
    rmse = np.sqrt(sum((data.obs - data.sim)**2) / data.shape[0])
    
    return rmse


def sesgo(observado, simulado):
    """Calcula el sesgo del hidrograma, es decir, el porcentaje de error en el volumen simulado.
    
    sesgo = (Vsim - Vobs) / Vobs * 100
    
    Parámetros:
    -----------
    observado:   series. Serie observada
    simulado:    series. Serie simulada"""
    
    # Eliminar pasos sin dato
    data = pd.concat((observado, simulado), axis=1)
    data.columns = ['obs', 'sim']
    data.dropna(axis=0, how='any', inplace=True)
    # Para la función si no hay datos
    assert data.shape[0] > 0, "ERROR. No indices match between the observed and the simulated series."
    # Calcular el sesgo    
    return (data.sim.sum() - data.obs.sum()) / data.obs.sum() * 100

  
def KGE(obs: pd.Series, sim: pd.Series, sa: float = 1, sb: float = 1, sr: float = 1) -> List[float]:
    """It computes the Kling-Gupta efficiency coefficient.
    
    Parameters:
    -----------
    obs:   pd.Series.
        Observed time series
    sim:   pd.Series
        Simulated time series
    sa, sb, sr: float
        Scale factors of the three terms of the modified KGE: ratio of the coefficient of variation (alpha), bias (beta), and coefficient of correlation (r), respectively
    
    Returns:
    -------
    KGE: float
        Modified KGE
    alpha: float
        Ratio of the standard deviations
    beta: float
        Bias, i.e., ratio of the mean values
    r: float
        Coefficient of correlation
    """
    
    # Eliminar pasos sin dato
    data = pd.concat((obs, sim), axis=1)
    data.columns = ['obs', 'sim']
    data.dropna(axis=0, how='any', inplace=True)
    # Para la función si no hay datos
    assert data.shape[0] > 0, "ERROR. No indices match between the observed and the simulated series."
    
    # calcular cada uno de los términos del KGE
    alpha = data.sim.std() / data.obs.std()
    beta = data.sim.mean() / data.obs.mean()
    r = np.corrcoef(data.obs, data.sim)[0, 1]
    
    # Cacular KGE
    ED = np.sqrt((sr * (r - 1))**2 + (sa * (alpha - 1))**2 + (sb * (beta - 1))**2)
    KGE = 1 - ED
    
    return KGE, alpha, beta, r



def KGEmod(obs: pd.Series, sim: pd.Series, sa: float = 1, sb: float = 1, sr: float = 1) -> List[float]:
    """It computes the modified Kling-Gupta efficiency coefficient.
    
    Parameters:
    -----------
    obs:   pd.Series.
        Observed time series
    sim:   pd.Series
        Simulated time series
    sa, sb, sr: float
        Scale factors of the three terms of the modified KGE: ratio of the coefficient of variation (alpha), bias (beta), and coefficient of correlation (r), respectively
    
    Returns:
    -------
    KGE: float
        Modified KGE
    alpha: float
        Ratio of the coefficient of variation
    beta: float
        Bias, i.e., ratio of the mean values
    r: float
        Coefficient of correlation
    """
    
    # Eliminar pasos sin dato
    data = pd.concat((obs, sim), axis=1)
    data.columns = ['obs', 'sim']
    data.dropna(axis=0, how='any', inplace=True)
    # Para la función si no hay datos
    assert data.shape[0] > 0, "ERROR. No indices match between the observed and the simulated series."

    # calcular cada uno de los términos del KGE
    alpha = (data.sim.std() / data.sim.mean()) / (data.obs.std() / data.obs.mean())
    beta = data.sim.mean() / data.obs.mean()
    r = np.corrcoef(data.obs, data.sim)[0, 1]
    
    # Cacular KGE
    ED = np.sqrt((sr * (r - 1))**2 + (sa * (alpha - 1))**2 + (sb * (beta - 1))**2)
    KGEmod = 1 - ED
    
    return KGEmod, alpha, beta, r
    


def matriz_confusion(obs, sim):
    """Calcula la matriz de confunsión del acierto en la ocurrencia o ausencia de precipitación diaria.
    
    Parámetros:
    -----------
    obs:       series. Serie observada
    sim:       series. Serie simulada"""
    
    data = pd.concat((obs, sim), axis=1)
    data.columns = ['obs', 'sim']
    # convertir días de lluvia en 1
    data[data > 0] = 1
    # Eliminar pasos sin dato
    data.dropna(axis=0, how='any', inplace=True)
    # Para la función si no hay datos
    if data.shape[0] == 0:
        FRND_logger.error('Series no coincidentes')
    # días con lluvia en la observación
    data1 = data[data.obs ==1]
    # días secos en la observación
    data0 = data[data.obs == 0]
    # calcular acierto
    acierto00 = sum(data0.sim == 0) / data0.shape[0]
    acierto01 = sum(data0.sim == 1) / data0.shape[0]
    acierto10 = sum(data1.sim == 0) / data1.shape[0]
    acierto11 = sum(data1.sim == 1) / data1.shape[0]
    acierto = [[acierto00, acierto01],
               [acierto10, acierto11]]
    acierto = pd.DataFrame(acierto, index=['obs0', 'obs1'], columns=['sim0', 'sim1'])
    
    return acierto


def logBiasNSE(obs, sim):
    """Calcula el coeficiente de eficiencia penalizado con el sesgo.
    
        logBiasNSE = NSE - 5 * abs(ln(1 + bias))**2.5
    
    Parámetros:
    -----------
    obs:       series. Serie observada
    sim:       series. Serie simulada"""
    
    # calcular Nash-Sutcliffe y error en volumen
    nse = NSE(obs, sim)
    bias = sesgo(obs, sim) / 100
    # calcular la función objetivo    
    return nse - 5 * abs(np.log(1 + bias))**2.5



def ECDF(serie: pd.Series, plot: bool = True, **kwargs) -> pd.Series:
    """it computes the empirical cumulative distribution function (ECDF) of the input data. If specified, the ECDF is plotted.

    Parametres:
    -----------
    serie:     pd.Series. 
        Input data
    plot:      bool
        Whether to plot or not the ECMWF

    Output:
    -------
    ecdf:      pd.Series
        The ECDF, where the index represents the non-exceedance probability and the values are those of the input data in ascending order
        
    Keyword arguments:
    ------------------
    lw:        float
        Line width for the line plot
    c:         str
        Line colour for the line plot
    xlim:      Tuple (2,)
        Limits of the X axis in the plot
    ylim:      Tuple (2,)
        Limits of the Y axis in the plot
    ylabel:    str
        Lable of the Y axis in the plot
    title:     str
        Text to be added as the plot title.
    """

    ecdf = pd.Series(data=serie.sort_values().values,
                     index=np.arange(1, serie.shape[0] + 1) / serie.count(),
                     name='ECDF')

    if plot:
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.plot(ecdf.index, ecdf, lw=kwargs.get('lw', 1), c=kwargs.get('c', 'steelblue'))
        ax.set(ylim=kwargs.get('ylim', (-1.02, 1.02)), ylabel=kwargs.get('ylabel', serie.name),
               xlim=(-.02, 1.02), xlabel='ECDF (-)',
               title=kwargs.get('title', ''));

    return ecdf



def PDE(obs: pd.Series, sim: pd.Series) -> float:
    """It computes the peak discharge error: mean of the absolute values of the relative error in the annual maximum discharge.
    
    Parameters:
    -----------
    obs:   pd.Series.
        Observed time series
    sim:   pd.Series
        Simulated time series
    
    Returns:
    -------
    PDE: float
        Peak discharge error
    """
    
    # remove timesteps with missing values
    data = pd.concat((obs, sim), axis=1)
    data.columns = ['obs', 'sim']
    data.dropna(axis=0, how='any', inplace=True)
    
    # annual maxima
    maxima = data.groupby(data.index.year).max()
    
    # peak discharge error
    PDE = ((maxima.obs - maxima.sim).abs() / maxima.obs).mean()
    
    return PDE



def pareto_front(of1: pd.Series, of2: pd.Series) -> pd.DataFrame:
    """It computes the Pareto front of two sets of objective functions.
    
    Inputs:
    -------
    of1: pd.Series
        Values of the objective function 1
    of2: pd.Series
        Values of the objective function 2
        
    Output:
    -------
    front: pd.DataFrame
        the pairs of values in 'of1' and 'of2' that delineate the Pareto front
    """
    

    df = pd.concat((of1, of2), axis=1).dropna(axis=0, how='any')
    cols = df.columns
    df.sort_values(by=cols[0], inplace=True)  

    # Iterate through each row in the sorted DataFrame
    pareto_front = []
    for index, row in df.iterrows():
        is_dominated = False

        # Check if the current row is dominated by any row in the Pareto front list
        for p_row in pareto_front:
            if row[cols[0]] >= p_row[cols[0]] and row[cols[1]] > p_row[cols[1]]:
                is_dominated = True
                break

        # If the current row is not dominated, add it to the Pareto front list
        if not is_dominated:
            pareto_front.append(row)

    # Convert the Pareto front list to a new DataFrame
    return pd.DataFrame(pareto_front)
