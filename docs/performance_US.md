# Comparison of the GloFAS reservoir simulation with reservoir records in the US
***

## Introduction

As a first approach to the improvement of the LISFLOOD reservoir routine, I have compared the results of the LISFLOOD reservoir simulations with available observations. The idea is to identify limitations of the current routine or parameterization and come up with ideas for improving the representation of reservoirs in LISFLOOD.

## Data

The main limitation for this analysis is the availability of reservoir observations. I will use the data set **ResOpsUS** [(Steyaert et al., 2022)](https://www.nature.com/articles/s41597-022-01134-7), a compilation of reservoir operation in the conterminous US. ResOpsUS contains daily time series of reservoir operations (storage, inflow, outflow, evapotranspiration and level) for 679 major reservoirs in the US.

The ResOpsUS time series will be compared with the daily time series simulated in the **GloFASv4** long run from 1982 to 2019. GloFASv4 represents 150 reservoirs in the US, but not all of them are included in ResOpsUS. 121 reservoirs are both represented in the ResOpsUS data set and GloFASv4, but in 8 of them the observed and simulated time series do not overlap. On top of that, there are 6 reservoirs I could not find their corresponding catchment in the GloFAS calibration, so I do not have the simulated time series. To sum up, the analysis below  comprises **113 reservoirs**.

### Uncertainty on the total reservoir capacity

I've noticed that the reservoir storage capacity in GloFAS (taken from GLWD), GRanD and ResOpsUS are note coherent. [Steyaert and Condon (2023)](https://hess.copernicus.org/preprints/hess-2023-194/) also found that the maximum capacity of 100 reservoirs in GRanD were not coherent with the observed records from ResOpsUS. As the reservoir routine in LISFLOOD is based on the fraction filled (the quotient of the storage and the total capacity), this uncertainty in the total reservoir capacity may affect the comparison of the simulated and observed storage.

To limit this uncertainty, I have compared the values of total capacity from three data sets (GloFAS, GRanD and ResOpsUS) and manually select for each reservoir whether GRanD or GloFAS is the data source more coherent when comparing against the observations from ResOpsUS. The plot below compares the maximum observed storage (ResOpsUS) against the total reservoir capacity reported both in GloFAS and GRanD. The Spearman correlation coefficient ($R^2$) gives an idea of the coherence between each pair of data sets.

<img src='../notebook/GloFAS/parameters/scatter_capacity.jpg' width="600">

***Figure 1**. Scatter plots of total reservoir capacity in three different data sets: ResOpsUS, GRanD and GloFAS.*

In general, the agreement between the 3 capacity data sets is good (correlation larger than 0.95). GRanD seems to agree better with ResOpsUS than GloFAS, but still some points deviate notably from the 1:1 relationship.

## Performance

### LISFLOOD model parameters

<font color='red'>Add the sketch of the reservoir routine</font>

![Calibrated parameters](../notebook/GloFAS/parameters/parameters_calibration_maps.jpg)

***Figure 2**. Values of the LISFLOOD reservoir parameters calibrated in GloFASv4. The size of the dots corresponds to the total reservoir capacity.*

![Reservoir limits](../notebook/GloFAS/parameters/parameters_maps.jpg)

***Figure 3**. Values of the LISFLOOD reservoir limits in GloFASv4. The size of the dots corresponds to the total reservoir capacity.*

### Reservoir storage

![Performance storage](../notebook/GloFAS/storage/maps_performance.jpg)

***Figure 4**. GloFASv4 performance in simulating reservoir storage. It shows the modified Klin-Gupta efficiency coefficient (KGE) and its three components: correlation (r), ratio of the coefficient of variance (alpha) and bias (beta). The size of the dots represents total reservoir capacity.*

#### Time series decomposition

### Reservoir inflow

![Performance inflow](../notebook/GloFAS/inflow/maps_performance.jpg)

***Figure 5**. GloFASv4 performance in simulating reservoir inflow. It shows the modified Klin-Gupta efficiency coefficient (KGE) and its three components: correlation (r), ratio of the coefficient of variance (alpha) and bias (beta). The size of the dots represents total reservoir capacity.*

#### Time series decomposition

### Reservoir outflow

![Performance outflow](../notebook/GloFAS/outflow/maps_performance.jpg)

***Figure 5**. GloFASv4 performance in simulating reservoir outflow. It shows the modified Klin-Gupta efficiency coefficient (KGE) and its three components: correlation (r), ratio of the coefficient of variance (alpha) and bias (beta). The size of the dots represents total reservoir capacity.*

#### Time series decomposition


### Comparison of the 3 variables

<img src='../notebook/GloFAS/behaviour/scatter_KGE.jpg' width='600'>

## Reservoir routine

<img src='../notebook/GloFAS/behaviour/016.jpg' width='675'>
<img src='../notebook/GloFAS/behaviour/140.jpg' width='675'>
<img src='../notebook/GloFAS/behaviour/227.jpg' width='675'>



## Ideas

* [ ] Fit the storage and outflow limits from the records.
* [ ] Regionalize those parameters to the reservoirs without records.
* [ ] Remove the fix value of 2 ($2 \cdot c_{lim}$) as the upper limit for the minimum outflow, and fit that multiplier of the conservative limit.
