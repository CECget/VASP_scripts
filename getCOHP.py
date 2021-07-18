#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat July 17 2021
for separating COHPCAR.lobster to two parts (bonding & antibonding)
Copyright: Yanyang Qin. XJTU
qinyanyang@stu.xjtu.edu.cn
@author: Yanyang Qin, XJTU
"""

######   separate COHPCAR.lobster to two parts (bonding & antibonding)   ######

import pandas as pd
####### reading lobster file #######
df_ini=pd.read_csv('COHPCAR.lobster',
                sep='\s+',
                skiprows=4,
                usecols=[0,1],
                names=['E (eV)','-COHP'],
                dtype='float64',
                )

######## reverse the sign ########
df_ini['-COHP']=df_ini['-COHP']*(-1)

######## create two copy for further calculation #########
df_bond = df_ini.copy()
df_antibond = df_ini.copy()
df_benchmark = pd.DataFrame(0, index=range(df_ini.shape[0]), columns = ['benchmark'])

####### replace the corresponding value with zero for ploting #######

df_bond.loc[df_bond['-COHP'] < 0, '-COHP'] = 0
df_antibond.loc[df_antibond['-COHP'] > 0, '-COHP'] = 0

######## contact the two table #########
df_output=pd.concat((df_bond,df_antibond['-COHP'],df_benchmark), axis=1)

######## save the final table as .dat file with "," as separator #########
df_output.to_csv('cookedCOHP.dat', sep=',', index=False)