# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 14:23:16 2019

@author: yx2957
"""

import geostatspy.GSLIB as GSLIB          # GSLIB utilies, visualization and wrapper
import geostatspy.geostats as geostats    # GSLIB methods convert to Python    
import os                                                   # to set current working directory 
import numpy as np                                          # arrays and matrix math
import pandas as pd                                         # DataFrames
from datetime import datetime
from datetime import timedelta
import matplotlib.pyplot as plt                             # plotting
import math as math
from mpl_toolkits.mplot3d import Axes3D

#%%
# Set the model parameters 

inj_min = 0.0; inj_max = 71000000.0; cmap = plt.cm.plasma                        # color min and max and using the plasma color map
xmin = 1733000; xmax = 2250000; ymin = 152300; ymax = 718900
mag_min = 2.5; mag_max = 3.0

#%%
# Read in all the data
os.chdir("C:/Users/yx2957/Documents/Research/Induced Seismicity/Jack_Induced Seismicity_USGS") 
well_loc = pd.read_csv("WellLocations.csv")  
well_inj = pd.read_csv("MichaelPyrcz_MonthlyInjVolumes_AfterQC_06222018.csv")
pressure = pd.read_csv("HighetP_All_TXNC.csv")
events = pd.read_csv("EarthQuakes.csv") 
nevent = len(events)
print(nevent)
events.head()


# For event, convert the date to a date formate
events['DATE'] = pd.to_datetime(events['DATE'])
print(events.head())
events.describe().transpose()

# Combine the year and month columns into a single date and remove the old columns
well_inj['DATE'] = well_inj.apply(lambda x:datetime.strptime("{0} {1} 00:00:00".format(x['YEAR'],x['MO']), "%Y %m %H:%M:%S"),axis=1)
well_inj = well_inj.drop('YEAR',axis=1)                                    # remove the zero column
well_inj = well_inj.drop('MO',axis=1)                                      # remove the zero column
well_inj.head()

# Assign the spatial extents 
xmin = 1733000; xmax = 2193000; ymin = 152300; ymax = 718900
# Set the UWI as the index / record key for the well locations and injections
well_loc = well_loc.set_index('UWI')
well_inj = well_inj.set_index('UWI')

# Add well coordinates to injection data
well_merge = well_inj.join(well_loc, how='outer')
cols = well_merge.columns.tolist()
cols = ['SURFX','SURFY','DATE','INJ']
well_merge = well_merge[cols] 
well_merge.head()

# Calculate the total injected at each well
well_summary = well_merge.groupby(['UWI'])
well_summary_total = well_summary.agg({'SURFX' : np.average, 'SURFY' : np.average, 'INJ' : np.sum,'DATE' : [np.min,np.max]})
well_summary_total.head()
print(len(well_summary_total.index))
print(well_summary_total)
print(well_summary_total.columns)
well_summary_total.to_csv('well_summary_total.csv')

#%%
# Make a grid of information on injection and events

nx = 10; ny = 10
total_injected = np.zeros((nx,ny))
max_event = np.zeros((nx,ny))
num_event = np.zeros((nx,ny))

#ninj = len(well_merge.index)
ninj = len(well_summary_total.index)

#ninj = 1000
nevent = len(events.index)


xsiz = (xmax - xmin)/nx; ysiz = (ymax - ymin)/ny
xcoords = np.arange(xmin,xmax, xsiz)
ycoords = np.arange(ymin,ymax, ysiz)

for iy in range(0,ny):
    for ix in range(0,nx):
        for inj in range(0,ninj):
            wx = float(well_summary_total.iloc[inj]['SURFX'])
            wy = float(well_summary_total.iloc[inj]['SURFY'])             
            if (wy >= ycoords[iy]) and (wy <= (ycoords[iy] + ysiz)):
                if (wx >= xcoords[ix]) and (wx <= (xcoords[ix] + xsiz)):
                    total_injected[ny-iy-1,ix] = total_injected[ny-iy-1,ix] + float(well_summary_total.iloc[inj]['INJ'])
                    if ix == 8 and iy == 2:
                        print(well_merge.iloc[inj]['INJ'],wx,wy) 
                        print(wx,wy)
#                        quit()

        for ievent in range(0,nevent): 
            ex = events.iloc[ievent]['SURFX']
            ey = events.iloc[ievent]['SURFY'] 
            smag = events.iloc[ievent]['MAG']
            if (ey >= ycoords[iy]) and (ey <= (ycoords[iy] + ysiz)):
                if ex >= xcoords[ix] and (ex <= (xcoords[ix] + xsiz)):
                    num_event[ny-iy-1,ix] = num_event[ny-iy-1,ix]+ 1
                    if max_event[ny-iy-1,ix] < smag: 
                        max_event[ny-iy-1,ix] = smag


#%%
# Make a grid of pressure information

nx = 10; ny = 10
max_pressure = np.zeros((nx,ny))
num_pres = np.zeros((nx,ny))
npres = len(pressure.index)

xsiz = (xmax - xmin)/nx; ysiz = (ymax - ymin)/ny
xcoords = np.arange(xmin,xmax, xsiz)
ycoords = np.arange(ymin,ymax, ysiz)

for iy in range(0,ny):
    print('Completed ' + str(iy) + ' of ' + str(ny))
    for ix in range(0,nx):
        for ipres in range(0,npres):
            wx = float(pressure.iloc[ipres]['SURF_X'])
            wy = float(pressure.iloc[ipres]['SURF_Y'])             
            if (wy >= ycoords[iy]) and (wy <= (ycoords[iy] + ysiz)):
                if (wx >= xcoords[ix]) and (wx <= (xcoords[ix] + xsiz)):
                    max_pressure[ny-iy-1,ix] = max_pressure[ny-iy-1,ix] + float(pressure.iloc[ipres]['HighestP'])
                    num_pres[ny-iy-1,ix] = num_pres[ny-iy-1,ix]+ 1

for iy in range(0,ny):
    for ix in range(0,nx): 
        if num_pres[ny-iy-1,ix] > 0.0:
            max_pressure[ny-iy-1,ix] = max_pressure[ny-iy-1,ix] / num_pres[ny-iy-1,ix]           

#%%

#combine pressure data with num_events and total_injected
            
GSLIB.ndarray2GSLIB(max_pressure,'max_pres_grid.out','max_pres_grid')

df_grid = pd.DataFrame({'max_pres': max_pressure.flatten(), 'num_event': num_event.flatten(), 'max_event': max_event.flatten(),'total_injected':total_injected.flatten()})
df_grid.loc[df_grid['num_event'] <= 0.0,'num_event'] = np.nan
df_grid.loc[df_grid['max_pres'] <= 0.0,'max_pres'] = np.nan
df_grid.loc[df_grid['max_event'] <= 0.0,'max_event'] = np.nan
df_grid.loc[df_grid['total_injected'] <= 0.0,'total_injected'] = np.nan
num_event_null = df_grid['num_event']
num_event_null = np.asarray(num_event_null)
max_event_null = df_grid['max_event']
max_event_null = np.asarray(max_event_null)
max_pressure_null = df_grid['max_pres']
max_pressure_null = np.asarray(max_pressure_null)
total_injected_null = df_grid['total_injected']
total_injected_null = np.asarray(total_injected_null)
df_grid.head()

#%%
# plot the image array
plt.figure(figsize=(4,4)) 
plt.subplots_adjust(left=0.0, bottom=0.0, right=3.0, top=2.5, wspace=0.2, hspace=0.2)

plt.subplot(221)
img = plt.imshow(np.reshape(max_pressure_null,(10, 10)),extent=[xmin,xmax,ymin,ymax])
plt.xlabel('X(ft)')
plt.ylabel('Y(ft)')
plt.title('Maximum Pore Pressure')
cbar = plt.colorbar(img, orientation = 'vertical')
cbar.set_label('Pressure (PSI)', rotation=270, labelpad=20)

plt.subplot(222)
img = plt.imshow(np.reshape(max_event_null,(10,10)),extent=[xmin,xmax,ymin,ymax])
plt.xlabel('X(ft)')
plt.ylabel('Y(ft)')
plt.title('Maximum Event Magnitude')
cbar = plt.colorbar(img, orientation = 'vertical')
cbar.set_label('Magnitude', rotation=270, labelpad=20)

plt.subplot(223)
img = plt.imshow(np.reshape(total_injected_null,(10,10)),extent=[xmin,xmax,ymin,ymax])
plt.xlabel('X(ft)')
plt.ylabel('Y(ft)')
plt.title('Total Injected Volume')
cbar = plt.colorbar(img, orientation = 'vertical')
cbar.set_label('Total Injected Volume', rotation=270, labelpad=20)

plt.subplot(224)
img = plt.imshow(np.reshape(num_event_null,(10,10)),extent=[xmin,xmax,ymin,ymax])
plt.xlabel('X(ft)')
plt.ylabel('Y(ft)')
plt.title('Number of Events')
cbar = plt.colorbar(img, orientation = 'vertical')
cbar.set_label('Number Events', rotation=270, labelpad=20)

#%%
#change the df.grid format; ready for Chi-square test of independence
#Set levels for number of events
df_grid = df_grid.fillna(0)
df_grid['level_events'] = ''
for i in range (len(df_grid['num_event'])):
    
    if df_grid['num_event'][i] >= float(4):
        df_grid['level_events'][i] = 'high'
        
    elif df_grid['num_event'][i] < float(4):
        df_grid['level_events'][i] = 'low'
        print('hi')
        
#%%
#Set levels for max pressure
df_grid['level_pres'] = ''
for i in range (len(df_grid['max_pres'])):
    
    if df_grid['max_pres'][i] < float(60):
        df_grid['level_pres'][i] = 'low'
        
    elif df_grid['max_pres'][i] <= float(100) and  df_grid['max_pres'][i] > float(60):
        df_grid['level_pres'][i] = 'mid'
        
    else:
        df_grid['level_pres'][i] = 'high'
#%%
#Set levels for injection volume
df_grid['level_inj'] = ''
for i in range (len(df_grid['total_injected'])):
    
    if df_grid['total_injected'][i] < float(.8*10**8):
        df_grid['level_inj'][i] = 'low'
        
    elif df_grid['total_injected'][i] >= float(.8*10**8):
        df_grid['level_inj'][i] = 'high'
        
#%%      
#extract a new df from df_grid 

#Compare pressure and events
df_Chi_test = pd.DataFrame({'events':df_grid['level_events'], 'pressure':df_grid['level_pres'], 'injection':df_grid['level_inj']})

#Contingency Table
contingency_table = pd.crosstab(index = df_Chi_test["pressure"], columns = df_Chi_test["events"])

#.apply(lambda r: r/r.sum(), axis=1)

print('contingency_table :-\n',contingency_table)

#Observed Values
Observed_Values = contingency_table.values 
print("Observed Values :-\n",Observed_Values)

#Expected Values
import scipy.stats
b=scipy.stats.chi2_contingency(contingency_table)
Expected_Values = b[3]
print("Expected Values :-\n",Expected_Values)


#Degree of Freedom
no_of_rows=len(contingency_table.iloc[:,0])
no_of_columns=len(contingency_table.iloc[0,:])
df=(no_of_rows-1)*(no_of_columns-1)
print("Degree of Freedom:-",df)

#Significance Level 5%
alpha=0.05

#chi-square statistic - χ2
from scipy.stats import chi2
chi_square=sum([(o-e)**2./e for o,e in zip(Observed_Values,Expected_Values)])
chi_square_statistic=chi_square[0]+chi_square[1]
print("chi-square statistic:-",chi_square_statistic)

#critical_value
critical_value=chi2.ppf(q=1-alpha,df=df)
print('critical_value:',critical_value)

#p-value
p_value=1-chi2.cdf(x=chi_square_statistic,df=df)
print('p-value:',p_value)

print('Significance level: ',alpha)
print('Degree of Freedom: ',df)
print('chi-square statistic:',chi_square_statistic)
print('critical_value:',critical_value)
print('p-value:',p_value)

#compare chi_square_statistic with critical_value and p-value which is the probability of getting chi-square>0.09 (chi_square_statistic)
if chi_square_statistic>=critical_value:
    print("Reject H0,There is a relationship between 2 categorical variables")
else:
    print("Retain H0,There is no relationship between 2 categorical variables")
    
if p_value<=alpha:
    print("Reject H0,There is a relationship between 2 categorical variables")
else:
    print("Retain H0,There is no relationship between 2 categorical variables")


#%%
#Compare injection and events
contingency_table = pd.crosstab(df_Chi_test['injection'],df_Chi_test["events"])
print('contingency_table :-\n',contingency_table)

#Observed Values
Observed_Values = contingency_table.values 
print("Observed Values :-\n",Observed_Values)

#Expected Values
import scipy.stats
b=scipy.stats.chi2_contingency(contingency_table)
Expected_Values = b[3]
print("Expected Values :-\n",Expected_Values)


#Degree of Freedom
no_of_rows=len(contingency_table.iloc[:,0])
no_of_columns=len(contingency_table.iloc[0,:])
df=(no_of_rows-1)*(no_of_columns-1)
print("Degree of Freedom:-",df)

#Significance Level 5%
alpha=0.05

#chi-square statistic - χ2
from scipy.stats import chi2
chi_square=sum([(o-e)**2./e for o,e in zip(Observed_Values,Expected_Values)])
chi_square_statistic=chi_square[0]+chi_square[1]
print("chi-square statistic:-",chi_square_statistic)

#critical_value
critical_value=chi2.ppf(q=1-alpha,df=df)
print('critical_value:',critical_value)

#p-value
p_value=1-chi2.cdf(x=chi_square_statistic,df=df)
print('p-value:',p_value)

print('Significance level: ',alpha)
print('Degree of Freedom: ',df)
print('chi-square statistic:',chi_square_statistic)
print('critical_value:',critical_value)
print('p-value:',p_value)   

#compare chi_square_statistic with critical_value and p-value which is the probability of getting chi-square>0.09 (chi_square_statistic)
if chi_square_statistic>=critical_value:
    print("Reject H0,There is a relationship between 2 categorical variables")
else:
    print("Retain H0,There is no relationship between 2 categorical variables")
    
if p_value<=alpha:
    print("Reject H0,There is a relationship between 2 categorical variables")
else:
    print("Retain H0,There is no relationship between 2 categorical variables")

#%%
#fisher's exact test for injection and events


import scipy.stats as stats

oddsratio, pvalue = stats.fisher_exact([[2, 7],[4,87]])



#%%
#Set levels for max pressure
df_grid['level_pres'] = ''
for i in range (len(df_grid['max_pres'])):
    
    if df_grid['max_pres'][i] <= float(80):
        df_grid['level_pres'][i] = 'low'
        
    elif df_grid['max_pres'][i] > float(80):
        df_grid['level_pres'][i] = 'high'
        
   
#Contingency Table
contingency_table = pd.crosstab(index = df_grid["level_pres"], columns = df_Chi_test["events"])
oddsratio2, pvalue2 = stats.fisher_exact(contingency_table)

#%%






































































