import numpy as np
import pandas as pd
import copy
import math
import string
from sklearn import metrics
import shutil
import os
import scipy
import random
from copy import deepcopy 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def linear_growth_model(x, y0, mu):
    return y0 + mu*x

def exponential_growth_model(x, y0, y1, mu):
    return y0 + y1 * np.exp(mu * x)

def Zwietering_Gompertz_growth_model(x, y0, A, lag, mu):
    return y0 + (A-y0)*np.exp(-np.exp(mu*np.exp(1)/A*(lag-x)+1))

def Zwietering_Logistic_growth_model(x, y0, A, lag, mu):
    return y0 + (A-y0)/(1+np.exp(4*mu/A*(lag-x)+2))

def get_r2_coefficient(function, popt, x, y):
    ypred = function(x, *popt)
    _, _, r_value, _, _ = scipy.stats.linregress(ypred, y)
    r_squared = r_value ** 2
    return r_squared

def get_slope(model, popt, x):
    if model == 'EXP':
        y0, y1, mu =  popt
        gr = y1*mu*np.exp(mu*x)
        sgr = gr/exponential_growth_model(x, *popt)        
    elif model == 'GOM':
        y0, A, lag, mu  = popt
        _exp = np.exp(mu*np.exp(1)/A*(lag-x)+1)
        gr = (A-y0)*np.exp(-_exp)*_exp*mu*np.exp(1)/A
        sgr = gr/Zwietering_Gompertz_growth_model(x, *popt)
    elif model == 'LGX':
        y0, A, lag, mu  = popt
        _exp = np.exp(4*mu/A*(lag-x)+2)
        gr = (A-y0)*_exp/((1+_exp)**2)*4*mu/A
        sgr = gr/Zwietering_Logistic_growth_model(x, *popt)
    else:
        print('unknown fitting model.')
        raise
        
    gr_to_return = [abs(max(gr)), np.mean(gr)]
    sgr_to_return = [abs(max(sgr)), np.mean(sgr)]
    return gr_to_return, sgr_to_return 

def fit_model_parameters(curve, phase, model):
    if model == 'EXP':
        init_guess = [random.random(), random.random(), random.random()]
        try:
            best_fits, covar = curve_fit(exponential_growth_model,
                                         curve.index,
                                         curve.values,
                                         p0=init_guess,
                                         bounds = ([0,-math.inf,-10], [math.inf,math.inf,10]),
                                         maxfev=1000
                                        )
            r2 = get_r2_coefficient(exponential_growth_model, best_fits, curve.index, curve.values)
        except:
            best_fits = [np.nan, np.nan, np.nan]
            r2 = 0.0
    elif model == 'GOM':
        init_guess = [random.random(), random.random(), random.random(), random.random()]
        try:
            best_fits, covar = curve_fit(Zwietering_Gompertz_growth_model,
                                         curve.index,
                                         curve.values,
                                         p0=init_guess,
                                         bounds = ([0,0,0,-10], [10,10,curve.index[-1],10]),
                                         maxfev=1000
                                        )
            r2 = get_r2_coefficient(Zwietering_Gompertz_growth_model, best_fits, curve.index, curve.values)
        except:
            best_fits = [np.nan, np.nan, np.nan, np.nan]
            r2 = 0.0
    elif model == 'LGX':
        init_guess = [random.random(), random.random(), random.random(), random.random()]
        try:
            best_fits, covar = curve_fit(Zwietering_Logistic_growth_model,
                                         curve.index,
                                         curve.values,
                                         p0=init_guess,
                                         bounds = ([0,0,0,-10], [10,10,curve.index[-1],10]),
                                         maxfev=1000
                                        )
            r2 = get_r2_coefficient(Zwietering_Logistic_growth_model, best_fits, curve.index, curve.values)
        except:
            best_fits = [np.nan, np.nan, np.nan, np.nan]
            r2 = 0.0
    else:
        print('unknown fitting model.')
        raise
    
    return best_fits, r2

def extract_growth_curve_features(file_od,
                                  file_phase,
                                  is_plot=True,
                                  plot_dim=(8,12),
                                  plot_logscale=False,
                                  xlim=None,
                                  ylim=None
                                 ):

    '''
    file_od: growth curves (OD measurement)
    file_phase: start time of each growth phase (I, II, III)
    is_plot: whether to plot curve fitting
    plot_dim: number of rows and columns in a plot. by default, it is 8 x 12
    plot_logscale: whether to plot in log scale
    xlim: x-axis range if is_plot is true
    ylim: y-axis range if is_plot is true
    '''
    #--------------
    # Read OD data
    #-------------- 
    if file_od.split('.')[-1] == 'xlsx':
        df_od = pd.read_excel(file_od, index_col=0)
    elif file_od.split('.')[-1] == 'csv':
        df_od = pd.read_csv(file_od, index_col=0)
    else:
        print('file format not recognized. Support xlsx and csv only.')
        raise
        
    #----------------------
    # Read phase start time
    #----------------------
    if file_phase.split('.')[-1] == 'xlsx':
        df_phase = pd.read_excel(file_phase, index_col=0)
    elif file_phase.split('.')[-1] == 'csv':
        df_phase = pd.read_csv(file_phase, index_col=0)
    else:
        print('file format not recognized. Support xlsx and csv only.')
        raise
        
    #---------------------------------
    # Reformat growth phase start time
    #---------------------------------
    line2append = []
    for strain in df_od.columns:
        phase_start_time_for_curr_strain = df_phase.loc[strain].values
        if not np.isnan(phase_start_time_for_curr_strain).any():
            line2append.append({'strain':strain,
                                'phase':'1',
                                'start_time':phase_start_time_for_curr_strain[0],
                                'end_time':phase_start_time_for_curr_strain[1],
                                })
            line2append.append({'strain':strain,
                                'phase':'2',
                                'start_time':phase_start_time_for_curr_strain[1],
                                'end_time':phase_start_time_for_curr_strain[2],
                                })
            line2append.append({'strain':strain,
                                'phase':'3',
                                'start_time':phase_start_time_for_curr_strain[2],
                                'end_time':df_od.index[-1],
                                })
    df_growth_phase = pd.DataFrame(line2append)

    #---------------
    # Fit parameters
    #---------------
    line2append = []
    fitting_models = ['EXP','GOM','LGX']
    for row_index in df_growth_phase.index:
        strain = df_growth_phase.loc[row_index,'strain']
        curr_curve = df_od[strain]
        phase_start_time = df_growth_phase.loc[row_index,'start_time']
        phase_end_time = df_growth_phase.loc[row_index,'end_time']
        curr_curve = curr_curve[(curr_curve.index >= phase_start_time) & (curr_curve.index <= phase_end_time)]
        if len(curr_curve.index) == 1:
            area = 0
        else:
            area = metrics.auc(curr_curve.index,curr_curve.values)
        #print(phase_start_time, phase_end_time)
        
        # fit phase curve using different models
        max_r2 = 0
        ntrial = 0
        max_trial = 50
        while max_r2 < 0.90 and ntrial < max_trial:
            best_fits = []
            r2 = []
            for model in fitting_models:
                curr_best_fits, curr_r2 = fit_model_parameters(
                    curr_curve, 
                    df_growth_phase.loc[row_index,'phase'],
                    model
                )
                best_fits.append(curr_best_fits)
                r2.append(curr_r2)
            max_r2 = max(r2)
            ntrial += 1
        
        # which model has best fitting
        best_model_index = np.argmax(r2)
                      
        # get growth rate and specific growth rate
        growth_rate, specific_growth_rate = get_slope(fitting_models[best_model_index],
                                                      best_fits[best_model_index],
                                                      curr_curve.index.values
                                                     )
                      
        # add data to output
        line2append.append({'strain':df_growth_phase.loc[row_index,'strain'],
                            'phase':df_growth_phase.loc[row_index,'phase'],
                            'start_time':df_growth_phase.loc[row_index,'start_time'],
                            'end_time':df_growth_phase.loc[row_index,'end_time'],
                            'delta_time':df_growth_phase.loc[row_index,'end_time']-df_growth_phase.loc[row_index,'start_time'],
                            'start_od':curr_curve.values[0],
                            'delta_od':curr_curve.values[-1]-curr_curve.values[0],
                            'area':area,
                            'best_model':fitting_models[best_model_index],
                            'r2':r2[best_model_index],
                            'growth_rate_max':growth_rate[0],
                            'growth_rate_mean':growth_rate[1],
                            'specific_growth_rate_max':specific_growth_rate[0],
                            'specific_growth_rate_mean':specific_growth_rate[1]
                           })
        
        # print
        print('strain=%s, phase=%s, start_time=%2.2f, start_od=%2.2f, area=%2.2f, best_model=%s (r2=%2.2f), growth_rate = [%2.2f, %2.2f], specific_growth_rate = [%2.2f, %2.2f]' %
              (line2append[-1]['strain'], line2append[-1]['phase'], line2append[-1]['start_time'], line2append[-1]['start_od'], line2append[-1]['area'],
               line2append[-1]['best_model'], line2append[-1]['r2'], line2append[-1]['growth_rate_max'], line2append[-1]['growth_rate_mean'], 
               line2append[-1]['specific_growth_rate_max'], line2append[-1]['specific_growth_rate_mean']))


    df_output = pd.DataFrame(line2append)          
    return df_output
