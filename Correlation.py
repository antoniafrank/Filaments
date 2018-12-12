#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 15:54:01 2018

@author: antoniafrank
"""

import math
import numpy as np
from numpy import sqrt, pi, exp, linspace, loadtxt

import matplotlib.pyplot as plt
from lmfit import Model



from sdas.tests.LoadSdasData import LoadSdasData

from sdas.tests.StartSdas import StartSdas

client = StartSdas()

import Correlation

#%%
def gaussian(x, amp, cen, wid):
    "1-d gaussian: gaussian(x, amp, cen, wid)"
    
    return (amp) * exp(-(x-cen)**2 /(2*wid**2))
        

#%%

def analyse(x,y,time,win,res,maxl=30,dist=3.,barrierset=1/math.e,barrier=True):

#%%    
    #window has to be higher or equal to res 
    points=int(len(x)/res) #data points sind vielfaches von meiner resolution
    tp=list(range(1,points+1))
    tp=[i * res for i in tp]
    tp=[i - (res/2) for i in tp]
    tp=[int(x) for x in tp] #for pytjon 2.7 not needed
    avt=time[tp]
    
#%%    
    cc=[]
    cc_err=[]
    
    lag=[]
    lag_err=[]
    
    v=[]
    v_err=[]
    
#    x=x-np.average(x) 
#    y=y-np.average(y)
               
    for i in range(0,points):
        
        t_min=int(tp[i]-(win/2))
        t_max=int(tp[i]+(win/2))
        
        if t_min < 0:
            cc=np.append(cc,np.nan)
            cc_err=np.append(cc_err,np.nan)
            
            lag=np.append(lag,np.nan)
            lag_err=np.append(lag_err,np.nan)
            
            v=np.append(v,np.nan)
            v_err=np.append(v_err,np.nan)
            
            i=i+1
            
        elif t_max > avt[points-1]:
            
            cc=np.append(cc,np.nan)
            cc_err=np.append(cc_err,np.nan)
            
            lag=np.append(lag,np.nan)
            lag_err=np.append(lag_err,np.nan)
            
            v=np.append(v,np.nan)
            v_err=np.append(v_err,np.nan)
            
       
        else:  
#%%            
           corrx=x[t_min:t_max]
           corry=y[t_min:t_max]
           
           corrx=corrx-np.average(corrx) 
           corry=corry-np.average(corry) 
           
#%%           
           #make the crosscorrelation
           
           #autocorrelation _a
           lags_a,c_a,line_a,b_a=plt.xcorr(corrx,corrx, usevlines=False , maxlags=maxl, linestyle='solid', marker=None)
           #crosscorrelation _c
           lags_c,c_c,line_c,b_c=plt.xcorr(corry,corrx, usevlines=False , maxlags=maxl, linestyle='solid', marker=None)
           #plt.close()
           
           #get lag and corrcoef
           cc1=np.max(c_c)
#%%           
           lag1=np.where(c_c==cc1)
           lag1=lag1[0]
           lag1=lag1[0] 
           lag2=lags_c[lag1]
           #in 0.5 mus steps
#           lag2=lag2/2.
#           
#           v1=(dist/lag2)


##%%
           
           fit_c=c_c[(lag1-10):(lag1+10)]
           fit_lags=lags_c[(lag1-10):(lag1+10)]
           fit_lags=fit_lags/2.
           #print(fit_c)
           #print(fit_lags)
           
#%%           
           # make a model that is a Gaussian + a constant:
           try:
               gmod = Model(gaussian)
             
               result = gmod.fit(fit_c, x=fit_lags, amp=cc1, cen=lag2/2., wid=1)
            
               #print(result.fit_report())
            
               plt.plot(fit_lags, fit_c,         'bo')
               plt.plot(fit_lags, result.init_fit, 'k--')
               plt.plot(fit_lags, result.best_fit, 'r-')
               plt.show()
               plt.close()
            
            
               cc1=result.best_values['amp']
               cc1_err=abs(result.params['amp'].stderr)
               
               lag2=result.best_values['cen']
               lag2_err=abs(result.params['cen'].stderr)
               
               v1=(dist/lag2)
               v1_err=abs((-(dist)/(lag2**2))*lag2_err)
               
               
           except TypeError:
               
               cc1=np.nan
               cc1_err=np.nan
               
               lag2=np.nan
               lag2_err=np.nan
               
               v1=np.nan
               v1_err=np.nan
               
#%%           
           
           
           
           if barrier==False:
           
#%%           
               #version1
               cc=np.append(cc,cc1)
               cc_err=np.append(cc_err,cc1_err)
               
               lag=np.append(lag,lag2)
               lag_err=np.append(lag_err,lag2_err)
               
               v=np.append(v,v1)
               v_err=np.append(v_err,v1_err)
               
#%%               
               i=i+1
              
 #%%         
            
           if barrier==True:
               
               
               #version2
                   
               if cc1 < barrierset:
                   
                   lag=np.append(lag,np.nan)
                   lag_err=np.append(lag_err,np.nan)
                   
                   v=np.append(v,np.nan)
                   v_err=np.append(v_err,np.nan)
                   
                   cc=np.append(cc,np.nan)
                   cc_err=np.append(cc_err,np.nan)
                               
                   i=i+1
                              
               else:    
                   
#                   lag1=np.where(c_c==cc1)
#                   lag1=lag1[0]
#                   lag1=lag1[0] 
#                   lag2=lags_c[lag1]
#                   #in 0.5 mus steps
#                   lag2=lag2/2.
#                   lag=np.append(lag,lag2)
#                       
#                   v1=(dist/lag2)
#                   v=np.append(v,v1)
#                   
#                   cc=np.append(cc,cc1)
#                   
#                   i=i+1
                   
                   cc=np.append(cc,cc1)
                   cc_err=np.append(cc_err,cc1_err)
               
                   lag=np.append(lag,lag2)
                   lag_err=np.append(lag_err,lag2_err)
               
                   v=np.append(v,v1)
                   v_err=np.append(v_err,v1_err)
               
                   i=i+1


    return cc, cc_err, lag, lag_err, v, v_err, avt


#%%
    
def evolution(dist,shot,win,res,maxl=30,barrierset=1/math.e,barrier=True):
    
    start=24
    end=17
    vel=[]
    d=int(dist/3.)
    
    for x in range(start,int(end+d-1),-1):
        
        I1, time = LoadSdasData(
                client, 'PCIE_ATCA_ADC_TEST.BOARD_1.CHANNEL_0'+str(x), 26494)
        I2, time = LoadSdasData(
                client, 'PCIE_ATCA_ADC_TEST.BOARD_1.CHANNEL_0'+str(x-d), 26494)
        
        cc,lag,v,avt=Correlation.analyse(I1,I2,time,win,res,maxl=maxl,dist=dist,
                                     barrierset=barrierset, barrier=barrier)
        
        vel=np.append(vel,v)
        
        x=x+1
    
    return vel



