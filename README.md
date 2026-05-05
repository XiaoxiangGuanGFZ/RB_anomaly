# Record-Breaking (RB) Anomaly

An R toolkit to compute and analyze record-breaking (RB) anomalies in climate time series, with a focus on hydrologic extremes.

This repository implements methods for detecting changes in the frequency of record-breaking events relative to a stationary climate. 


# Overview

Extreme climate events, especially record-breaking precipitation and temperature, have increased in frequency under global warming. Instead of relying on parametric extreme value theory, the record-breaking framework offers a robust, distribution-free approach.

A value is considered record-breaking if it exceeds all previous values in a time series. The key metric implemented here is the Record-Breaking Anomaly (RB anomaly):

$$
RB_{\text{anom}} = \frac{
    R_{\text{obs}} - R_{\text{exp}}
}{R_{\text{exp}} } \times 100 \%
$$

where:

- $R_{\text{obs}}$: observed number of record-breaking events
- $R_{\text{exp}} = \sum_{n=1}^{N}{\frac{1}{n}}$: expected number under a stationary (iid) process.

This method allows direct comparison between observations and a null hypothesis of no climate trend.

# Features

- Compute record-breaking events for time series (e.g., precipitation, temperature)
- Estimate RB anomaly relative to iid expectation
- Perform Monte Carlo shuffling to generate confidence intervals
- Spatial aggregation across regions


# Reference
1. Lehmann, J., Coumou, D. and Frieler, K.  2015.  Increased record-breaking precipitation events under global warming. Climatic Change, 132(4), 501-515. doi: 10.1007/s10584-015-1434-y.
2. Robinson, A., Lehmann, J., Barriopedro, D., Rahmstorf, S. and Coumou, D.  2021.  Increasing heat and rainfall extremes now far outside the historical climate. npj Climate and Atmospheric Science, 4(1), 45. doi: 10.1038/s41612-021-00202-w.

# Contact

[Xiaoxiang Guan](guan@gfz.de)

[GFZ Helmholtz Centre for Geosciences](https://www.gfz.de)

