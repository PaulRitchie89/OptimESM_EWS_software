# OptimESM_EWS_software
Early warning signal software for use on the OptimESM project

Running the python code 'EWS_detection_box_mean.py' will compute and plot early warning signals (increasing lag-1 autocorrelation and variance) for an area weighted average of the barotropic streamfunction in the North Atlantic from the NorESM2-LM model under the historical + ssp126 scenario. Data provided in the file 'msftbarot_NorESM2-LM_historical_ssp126_r1i1p1f1_1850_2100_Atlantic_corrected.mat'. 

The code uses the previously created python packages 'regime_shifts.py' and 'ews.py' (with slight modifications) from https://github.com/BeatrizArellano/regimeshifts, where further information on the calculation of the early warning signals (EWSs) can be found. Also calculated in the code are statistics for the EWS such as the Kendall tau value (a measure of a time series monotonicity), and a metric that determines if the signal is statistically significant by calculating a 95% confidence interval for the variability of the EWS from the control run data 'msftbarot_NorESM2-LM_piControl_r1i1p1f1_0001_0501_Atlantic_corrected.mat'.
