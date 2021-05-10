#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Python3
Created on Mon 1 Mar 2021
for calculating the d/p/pz band center of the selected atoms
Copyright:  Yanyang Qin. XJTU
750881345@qq.com
@author: Yanyang Qin, XJTU
"""

######   import tools   ######
import numpy as np
import datetime
import time
import pandas as pd
from scipy.integrate import simps

print("Please enter the index of the selected atoms :")
select_index = np.array([int(n) for n in input().split()])

######   timing   ######
start = time.time()
print('\n********** Analyzing the selected DOS from Yanyang Qin XJTU **********')
print('\n********** Calculating ... \n')
### current time ###
start_time = datetime.datetime.now()
print("Start time:       " + start_time.strftime('%Y.%m.%d-%H:%M:%S'))   #strftime可以自定义时间的输出格式
np.seterr(divide='ignore',invalid='ignore')

select_atoms_sum = pd.read_csv("DOS"+str(select_index[0])+"_splitted.dat",
                   sep=' ',
                   index_col=0)


for i in range(len(select_index)-1):
    add_data = pd.read_csv("DOS"+str(select_index[i+1])+"_splitted.dat",
                   sep=' ',
                   index_col=0)
    select_atoms_sum = select_atoms_sum + add_data

select_atoms_sum.to_csv('select_atoms_sum_split.dat', sep=' ', float_format='%12.10f')
########## no split ###############
orbital_group = [0,1,2,3,2,3,2,3,4,5,4,5,4,5,4,5,4,5]
select_atoms_sum_nosplit = select_atoms_sum.groupby(orbital_group,axis=1).sum()
select_atoms_sum_nosplit.columns=['s_up', 's_down', 'p_up', 'p_down', 'd_up', 'd_down']
select_atoms_sum_nosplit.to_csv('select_atoms_sum.dat', sep=' ', float_format='%12.10f')

########## get bandcenter & integrate charge #######
integrate_index = select_atoms_sum_nosplit[select_atoms_sum_nosplit.index<0.05].index.tolist() # choose the integrate energy range
integrate_range = len(integrate_index)
integrate_df = select_atoms_sum_nosplit.loc[integrate_index]

x = integrate_index
y_dup = np.array(integrate_df.iloc[:,4].tolist())   # get data as numpy array
y_ddown = np.array(integrate_df.iloc[:,5].tolist())
y_pup = np.array(integrate_df.iloc[:,2].tolist())
y_pdown = np.array(integrate_df.iloc[:,3].tolist())

dbc_up = simps(y_dup*x, x) / simps(y_dup, x)  #integrate to get d-bandcenter
dbc_down = simps(y_ddown*x, x) / simps(y_ddown, x)

pbc_up = simps(y_pup*x, x) / simps(y_pup, x)  #integrate to get p-bandcenter
pbc_down = simps(y_pdown*x, x) / simps(y_pdown, x)

d_charge_integrate_up = simps(y_dup,x)   #integrate to get the quantity of d-charges
d_charge_integrate_down = simps(y_ddown,x)
d_charge_integrate_total = d_charge_integrate_up - d_charge_integrate_down


p_charge_integrate_up = simps(y_pup,x)   #integrate to get the quantity of p-charges
p_charge_integrate_down = simps(y_pdown,x)
p_charge_integrate_total = p_charge_integrate_up - p_charge_integrate_down

########## write output ##########
reults_text = ('Band centers:\n')+('\tdbc_up(eV)\tdbc_down(eV)\n\
    %-12.10f\t%-12.10f\n' %(dbc_up,dbc_down))+\
    ('\tpbc_up(eV)\tpbc_down(eV)\n\
    %-12.10f\t%-12.10f\n' %(pbc_up,pbc_down))+\
    ('d-integrate charge:\n\tspin_up\t\tspin_down\t\tabsolute_total\n\
    %-12.10f\t%-12.10f\t\t%-12.10f\n' %(d_charge_integrate_up,d_charge_integrate_down,d_charge_integrate_total))+\
    ('p-integrate charge:\n\tspin_up\t\tspin_down\t\tabsolute_total\n\
    %-12.10f\t%-12.10f\t\t%-12.10f' % (p_charge_integrate_up, p_charge_integrate_down, p_charge_integrate_total))

########## save output ###########
with open('analysis_results.txt','w') as results_save:
    results_save.write(reults_text)
##########   timing   #############
stop=time.time()
print("running time:     " + str(stop-start) + " seconds")
terminal_time = datetime.datetime.now()
print("Terminal time:    " + terminal_time.strftime('%Y.%m.%d-%H:%M:%S'))   #strftime可以自定义时间的输出格式
print('\n********** Completed !!!\n')
print(reults_text)
print('\n********** Analyzing the selected DOS from Yanyang Qin XJTU **********\n\n')
