#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 2020
for vabrational frequency calculating the correction terms of G
@author: Yanyang Qin, XJTU
"""
# Ref:https://www.bigbrosci.com/2018/11/07/ex69/ 
# email: 750881345@qq.com

import sys
import math
from scipy import constants as con
import linecache as lce

#########################Create the result file with raw frequency data########################################
out=open('tsresult',mode='w')
#natom=2
#str_natom='%d' %natom
#keyword=str_natom+' '+'f'
keyword='2PiTHz'
with open('OUTCAR',mode='r') as getdatafile:
	for line in getdatafile:
		try:
			line_str=line.strip( ).split()
			if '=' in line_str[1]:
				if keyword in line_str[5]:
                                	out.write(line)	
			else:
				if keyword in line_str[6]:
					out.write(line)
		except:
			continue

out.close()
#########################Calculater######################################################
with open('tsresult','rU') as count:
    num=len(count.readlines())
with open('tsresult',mode='a') as rcal:
    rcal.write('\n\n##############################\n%-15s\t%-15s\t%-15s\t%-15s\t%-15s\t\n' %('DOFnumber','S_J/K/mol','TS','Cpi','EmeV'))
Sum_S=float(0)
Sum_TS = float(0)
Sum_meV=float(0)
Sum_Cp=float(0)
for i in range(1,num+1):
    line_str=lce.getline('tsresult',i)
    str=line_str.split()
    if '=' in str[1]:
        wavenumber=str[6]
    EmeV=str[8]
    else:
        wavenumber=str[7]
	EmeV=str[9]
    nu = float(wavenumber) * 100  # convert unit from cm-1 to m-1
    h = con.h  # get Plank Constant
    k = con.k  # get Boltzman Constant
    R = con.R  # get Gas Constant
    c = con.c  # get Lightspeed Constant
    T = 300  # Set Temperature as 300K
    beta = 1 / (k * T)
    def get_pf(nu):
        x = beta * h * c * float(nu)  # Numerator of the first term of the equation
        pf1 = x / (math.exp(x) - 1)  # First term of the equation
        pf2 = math.log(1 - math.exp(-x))  # Second term *
        pf = pf1 - pf2
        entropy = R * pf
        return entropy
    entropy = get_pf(nu)  # unit is J*K-1*mol-1
    TS = entropy * T / 1000 / 96.485  # unit is eV
    nu = float(nu) / 100  # convert unit back
    epsilon_i=h*c*float(nu)*100
    expot_i=epsilon_i*beta
    PCapacity_i=epsilon_i/con.e/(math.exp(expot_i)-1)
    if '=' in str[1]:
	Sum_S=Sum_S
	Sum_TS=Sum_TS
	Sum_meV=Sum_meV
        Sum_Cp=Sum_Cp
    else:
        Sum_S=Sum_S+entropy
        Sum_TS=Sum_TS+TS
	Sum_meV=Sum_meV+float(EmeV)
        Sum_Cp=Sum_Cp+PCapacity_i
    with open('tsresult',mode='a') as rcal:                 #write into the result file
	rcal.write('%-15s\t%-15.4f\t%-15.4f\t%-15.4f\t%-15.4f\t\n' %(str[0],entropy,TS,PCapacity_i,float(EmeV)))
    print(str)
    print('%-15s\t%-15.4f\t%-15.4f\t%-15.4f\t%-15.4f\t' %(str[0],entropy,TS,PCapacity_i,float(EmeV)))
ZPEeV=float(Sum_meV)/2000
Dif=ZPEeV+Sum_Cp-Sum_TS 
with open('tsresult',mode='a') as rcal:                    #write the sumup value       
    rcal.write('##############################\n%-15s\t%-15s\t%-15s\t%-15s\t%-15s\t\n%-15s\t%-15.4f\t%-15.4f\t%-15.4f\t%-15.4f\t\n%-15s\t%-15.4f\n' %('/','Sum_S(J/K/mol)','Sum_TS(eV)','CvT(eV)','ZPE(eV)','Sum/Final',Sum_S,Sum_TS,Sum_Cp,ZPEeV,'TotDiff is',Dif))
