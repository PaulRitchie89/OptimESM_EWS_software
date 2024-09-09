# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 10:15:04 2024

@author: Paul Ritchie

Script to calculate EWS (increasing lag-1 autocorrelation and variance) for a
given time series containing an abrupt shift. EWS statistics are measured by
the Kendall tau (indicates the monotonicity of a time series) and the warning
time defined to be the period of time the signal remains continuously above
a 95% confidence interval determined from a provided control run.

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import seaborn as sns
import numpy.ma as ma
from scipy.stats import norm
from matplotlib import rc

from regimeshifts import regime_shifts as rs
from regimeshifts import ews

fontsize = 14
rc('font', **{'size' : fontsize})
rc('text', usetex=True)

def sig_test(ts,control_ts):
    '''
    Function to calculate the confidence interval for EWS statistic based on the
    control run and determine the corresponding warning time.

    Parameters
    ----------
    ts : DataFrame
        EWS time series up to abrupt shift
    control_ts : Dataframe
        EWS time series of control run

    Returns
    -------
    std_err : float64
        Reutrns standard error for 95% confidence interval
    time_ind : int
        Returns warning time EWS statistic provides

    '''
    ts = ts.loc[ts.first_valid_index():ts.last_valid_index()]
    time_ind = np.nan
    std_err = norm.ppf(0.95)*np.sqrt(control_ts.loc[:,'Time series'].var(ddof=0))
    for j in range(len(ts)):
        # Identify first instance that all remaining points of the change in the
        # EWS statistic are above the confidence interval, and record how far
        # from end of time series to determine warning time
        if all(i >= float(std_err) for i in ts.loc[ts.first_valid_index()+j:,'Time series']-ts.loc[ts.first_valid_index(),'Time series']):
            time_ind = len(ts) - j
            break
    return (std_err,time_ind)


############### Specify Gaussian bandwidth for detrending ################
bW = 60

############### Specify window length to calculate EWS over ################
wL = 100

######### Specify minimum length of time required to monitor EWS for #########
margin = 40

################# Specify variable of interest ######################
var = 'msftbarot'

################ Specify region of interest #######################
region = 'Atlantic'

################# Specify model details here ################################## 
models = 'NorESM2-LM'#,'CESM2-WACCM','MRI-ESM2-0'
variant_ids = 'r1i1p1f1'
date_ends = '2100'
control_ends = '0501'#,'0499','0701'

####################### Specify experiments ##################################
experiment_hist = 'historical'
experiment_ssps = 'ssp126'#'ssp245'

################### Initialise time interval ########################
t = np.linspace(1850,2100,251)

############## Initialise figure ####################
fig, ax = plt.subplots(4,1,sharex=True,figsize=(6.4,8.5))

# ax[0].set_ylabel('Sea surface salinity',fontsize=fontsize)
ax[0].set_ylabel('Change in barotropic\nstreamfunction (Sv)',fontsize=fontsize)
# ax[0].set_ylabel('Change in air\ntemperature ($^o$C)',fontsize=fontsize)
# ax[0].set_ylabel('Mixed layer depth (m)',fontsize=fontsize)
ax[1].set_ylabel('Detection index',fontsize=fontsize)
ax[2].set_ylabel('AR(1)',fontsize=fontsize)
ax[3].set_ylabel('Variance',fontsize=fontsize)
ax[3].set_xlabel('Year',fontsize=fontsize)

ax[0].set_xlim(t[0],t[-1])


    
   
#### File names containing data of SSP and PI control runs ####
fname = 'C:/Users/pdlr201/OneDrive - University of Exeter/OptimESM/CMIP6_data/'+var+'/'+experiment_ssps+'/Processed_data/'+region+'/'+var+'_'+models+'_'+experiment_hist+'_'+experiment_ssps+'_'+variant_ids+'_1850-'+date_ends+'_'+region+'_corrected.mat'
fname2 = 'C:/Users/pdlr201/OneDrive - University of Exeter/OptimESM/CMIP6_data/'+var+'/piControl/Processed_data/'+region+'/'+var+'_'+models+'_piControl_'+variant_ids+'_0001-'+control_ends+'_'+region+'_corrected.mat'

LONS = loadmat(fname)['LONS'][0]    # Retrieve longitudes
LATS = loadmat(fname)['LATS'][0]    # Retrieve latitudes
x_ts = loadmat(fname)['region_data']    # Retrieve time series data of SSP run for all grid cells

x_ts_control = loadmat(fname2)['region_data']   # Retrieve time series data of control run for all grid cells

weights = np.cos(np.deg2rad(LATS))  # Calculate weigths for each grid cell

#### Initialise time series for weigthed averages of SSP and control runs ####
x_ts_mean = np.zeros(251)
x_ts_control_mean = np.zeros(499)

#### Calculate weigthed averages of SSP and control runs ####
for i in range(251):
    x_ts_mean[i] = ma.average(x_ts[i,:].data,weights=weights)
for i in range(499):
    x_ts_control_mean[i] = ma.average(x_ts_control[i,:].data,weights=weights)

### Run Regime_shift to generate a Panda Series to implement a regime shift detection index ###
ts = rs.Regime_shift(x_ts_mean)
detection = ts.as_detect()  # Calculate detection time series
detection_index = detection[np.argmax(np.abs(detection))]   # Identify index of largest maximum or minimum

### Retrieve the time series only up to the largest abrupt change ###     
if detection_index < 0:
    bef_rs = ts.before_drs()
else:
    bef_rs = ts.before_urs()

### Plot the change in time series and detection time series in first 2 panels ###
ax[0].plot(t,ts-ts[0])
ax[1].plot(t,detection)

### The Ews class returns an extended Dataframe object, if we provided a series, it sets 0 for the column name ###
series = ews.Ews(bef_rs)
series = series.rename(columns={0:'Time series'})  

### Calculate EWS if time series is of sufficient length to monitor EWS ###
if (len(bef_rs) > wL + margin):
          
    ar1 = series.ar1(detrend=True,bW=bW,wL=wL)  # Compute lag-1 autocorrelation using the ar1() method
    vari = series.var(detrend=True,bW=bW,wL=wL) # Compute variance
    
    ar1_tau = ar1.kendall   # Compute Kendall tau of lag-1 autocorrelation
    var_tau = vari.kendall  # Compute Kendall tau of variance
    
    ### Create Dataframe object for control run ###
    ts2 = rs.Regime_shift(x_ts_control_mean)
    series2 = ews.Ews(ts2)
    series2 = series2.rename(columns={0:'Time series'})
    
    control_ar1 = series2.ar1(detrend=True,bW=bW,wL=wL) # Compute control lag-1 autocorrelation using the ar1() method
    control_var = series2.var(detrend=True,bW=bW,wL=wL) # Compute control variance
    
    ar1_err, ar1_time = sig_test(ar1,control_ar1)   # Compute warning time for lag-1 autocorrelation
    var_err, var_time = sig_test(vari,control_var)  # Compute warning time for variance   
    
    ####### Complete plotting of EWS ##########
    ax[2].plot(t[:len(ar1)],ar1,label='$\\tau$ = '+"{:.2f}".format(ar1_tau))
    ax[3].plot(t[:len(vari)],vari,label='$\\tau$ = '+"{:.2f}".format(var_tau))
    ax[2].fill_between([1850,2100],ar1.loc[wL-1,'Time series']-ar1_err,ar1.loc[wL-1,'Time series']+ar1_err,alpha=0.3)
    ax[3].fill_between([1850,2100],vari.loc[wL-1,'Time series']-var_err,vari.loc[wL-1,'Time series']+var_err,alpha=0.3)


ax[2].legend(frameon=False)
ax[3].legend(frameon=False)


fig.subplots_adjust(left=0.175, right=0.95, bottom=0.075, top=0.975, hspace=0.1)

sns.despine()

plt.text(-0.21,0.95,'$\\textbf{(a)}$',transform=ax[0].transAxes)
plt.text(-0.21,0.95,'$\\textbf{(b)}$',transform=ax[1].transAxes)
plt.text(-0.21,0.95,'$\\textbf{(c)}$',transform=ax[2].transAxes)
plt.text(-0.21,0.95,'$\\textbf{(d)}$',transform=ax[3].transAxes)