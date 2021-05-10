#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Copyright: 24th, December, Yaqiong Su, Eindhoven
Y-Q.Su@tue.nl; 398168203@qq.com

Modified by Yanyang Qin in 2020
750881345@qq.com

@author: Yaqiong Su, Eindhoven; Yanyang Qin, XJTU
"""

######   split the DOSCAR to each atom   ######
import numpy as np
import datetime
import time
import sys
import linecache as lce
import pandas as pd
from scipy.integrate import simps
import os

np.set_printoptions(suppress=True)
atom_ini, atom_final = sys.argv[1],sys.argv[2]   # pass parameters  atom range

######   timing   ######
start = time.time()
print '********** splitted dos and d-band center from Yaqiong Su Eindhoven **********'
print 'is getting splitted dos and d-band center'
### current time ###
start_time = datetime.datetime.now()
print "Start time:       " + start_time.strftime('%Y.%m.%d-%H:%M:%S')   #strftime可以自定义时间的输出格式

######   reading DOSCAR   ######
#f1 = open('header','wb')
#with open('DOSCAR','rb') as f:
#    i = 0
#    while True:
#        line = f.readline()
#        i+=1
#        if i == 1:
#            atoms_number = line.split()[0]
#        if i == 6:
#            efermi = line.split()[3]
#            epoints = line.split()[2]
#        if i < 7:
#            f1.write(line)
#        if i > 1e5:
#            break
#f1.close()

line1 = lce.getline('DOSCAR', 1)
atoms_number = int(line1.split()[0])
line6 = lce.getline('DOSCAR', 6)
efermi = float(line6.split()[3])
epoints = int(line6.split()[2])
linec = lce.getline('DOSCAR', 6 + epoints + 1 + 1)
column_num = int(len(linec.split()))
print "atoms number:",   atoms_number
print "Fermi energy:",   efermi
print "Energy points:",   epoints
print "column number:",   column_num

f2 = open('total_dos.dat','wb')
f3 = open('atoms_dos0.dat','wb')
f4 = open('atoms_dos.dat','wb')
with open('DOSCAR', 'rb') as f:
    i = 0
    while True:
        i+=1
        line = f.readline()
        if i > 6 and i < 6+epoints+1:
            f2.write(line)
        if i > 6+epoints and i < 1e5:
            f3.write(line)
            if '301   ' in line:
                continue
            f4.write(line)
        if i > 1e5:
            break
f2.close()
f3.close()
f4.close()

Ecorrected = np.array([float(l.split()[0]) for l in open('total_dos.dat','rb')])-efermi
Ecorrected = Ecorrected.reshape(epoints,1)
rows = np.zeros([epoints,column_num-1])
Ef = np.zeros([epoints,1])-efermi
Ef0 = np.hstack((Ef,rows))
ones = np.ones([epoints,column_num])
for i in range(2,column_num,2):
    ones[:,i] = -1   # construct 1 -1 matrix

names = [ 's', 'p_y', 'p_z', 'p_x', 'd_xy', 'd_yz', 'd_z2-r2', 'd_xz', 'd_x2-y2' ]
all_names = []
for n in names:
    all_names.extend( [ '{}_up'.format(n), '{}_down'.format(n) ] )
all_names.insert( 0, 'energy(eV)' )

names0 = [ 's', 'p', 'd' ]
group_names = []
for n in names0:
    group_names.extend( [ '{}_up'.format(n), '{}_down'.format(n) ] )
group_names.insert( 0, 'energy(eV)' )
#print all_names

atoms_dos = np.loadtxt('atoms_dos.dat')
dos = []
for i in range(1,atoms_number+1):
    dos.append('DOS' + str(i))
    dos[i-1] = np.vsplit(atoms_dos,atoms_number)[i-1].astype(float) # vertical split
    dos[i-1] = np.add(dos[i-1],Ef0) # corrected energy by efermi
    dos[i-1] = np.multiply(dos[i-1],ones)
    locals()['DOS' + str(i)] = dos[i-1]
    file = 'DOS' + str(i) + '.dat'
    np.savetxt(file,dos[i-1])   # projected to each atomic orbital orientation
    fi = 'f' + str(i)
    fi = np.loadtxt(file)
    dfi = 'df' + str(i)
    dfi = pd.DataFrame(fi,columns = all_names)
    dfi.set_index('energy(eV)',inplace=True)
    file_splitted = 'DOS' + str(i) + 'splitted' + '.dat'
    dfi.to_csv(file_splitted,sep=' ',float_format='%12.10f') # projected to each atomic orbital orientation
######   sum of orbital dos   ######
    Eni = 'En' + str(i)
    Eni = dos[i-1][:,0].reshape(epoints,1)
    locals()['En' + str(i)] = Eni
    s_upi = 's_up' + str(i)
    s_upi = dos[i-1][:,1].reshape(epoints,1)
    locals()['s_up' + str(i)] = s_upi
    s_downi = 's_down' + str(i)
    s_downi = dos[i-1][:,2].reshape(epoints,1)
    locals()['s_down' + str(i)] = s_downi
    p_upi = 'p_up' + str(i)
    p_upi = (dos[i-1][:,3]+dos[i-1][:,5]+dos[i-1][:,7]).reshape(epoints,1)
    locals()['p_up' + str(i)] = p_upi
    p_downi = 'p_down' + str(i)
    p_downi = (dos[i-1][:,4]+dos[i-1][:,6]+dos[i-1][:,8]).reshape(epoints,1)
    locals()['p_down' + str(i)] = p_downi
    d_upi = 'd_up' + str(i)
    d_upi = (dos[i-1][:,9]+dos[i-1][:,11]+dos[i-1][:,13]+dos[i-1][:,15]+dos[i-1][:,17]).reshape(epoints,1)
    locals()['d_up' + str(i)] = d_upi
    d_downi = 'd_down' + str(i)
    d_downi = (dos[i-1][:,10]+dos[i-1][:,12]+dos[i-1][:,14]+dos[i-1][:,16]+dos[i-1][:,18]).reshape(epoints,1)
    locals()['d_down' + str(i)] = d_downi
    groupi = 'group' + str(i)
    groupi = np.hstack((Eni,s_upi,s_downi,p_upi,p_downi,d_upi,d_downi))
    locals()['group' + str(i)] = groupi
#    print group1
    np.savetxt(file,groupi)   # projected to each atomic orbital
    fi = 'f' + str(i)
    fi = np.loadtxt(file)
    dfi = 'df' + str(i)
    dfi = pd.DataFrame(fi, columns = group_names)
    dfi.set_index('energy(eV)',inplace=True)
    dfi.to_csv(file,sep=' ',float_format='%12.10f') # projected to each atomic orbital
######   sum of dos   ######
    dos_groupi = 'dos_group' + str(i)
    dos_groupi = np.hstack((s_upi,s_downi,p_upi,p_downi,d_upi,d_downi))
    locals()['dos_group' + str(i)] = dos_groupi
    filedos = 'dos' + str(i)
    np.savetxt(filedos,dos_groupi)

dos_sum = np.zeros([epoints,6])
#print dos_sum.shape
for j in range(int(atom_ini),int(atom_final)+1):
    dos_groupj = 'dos' + str(j)
    locals()['dos' + str(j)] = dos_groupj
    foj = 'fo' + str(j)
    foj= np.loadtxt(dos_groupj)
#    print foj.shape
    dos_sum = np.add(dos_sum,foj)
dos_sum = np.hstack((Ecorrected,dos_sum))
ddos = pd.DataFrame(dos_sum, columns = group_names)
ddos.set_index('energy(eV)',inplace=True)
ddos.to_csv('sum_dos.dat',sep=' ',float_format='%12.10f') # pdos sum of selected atoms


###### calculaton of p-band center ######
#file = sys.argv[1] # pass paramter
file = 'sum_dos.dat'  # the selected atom
list_energy = [l.split()[0] for l in open(file,'rb')]
list_energy.pop(0)
energy = np.array([float(i) for i in list_energy])  # extract energy

emin, emax = energy[0], 0.05   # integral energy range
erange = (energy[0],energy[-1])
emask = (energy >= emin) & (energy <= emax) # bool to make a mapping between energy and dos

list_pup = [l.split()[3] for l in open(file,'rb')]
list_pup.pop(0)
p_up = np.array([float(i) for i in list_pup])   # extract p_up
list_pdown = [l.split()[4] for l in open(file,'rb')]
list_pdown.pop(0)
p_down = np.array([float(i) for i in list_pdown])   # extract p_down

x = energy[emask]
y1 = p_up[emask]
y2 = p_down[emask]

pbc_up   = simps(y1*x, x) / simps(y1, x)
pbc_down = simps(y2*x, x) / simps(y2, x)
pbc = []
pbc.append(pbc_up)
pbc.append(pbc_down)

###### calculaton of d-band center ######

## same energy set ##

list_dup = [l.split()[5] for l in open(file,'rb')]
list_dup.pop(0)
d_up = np.array([float(i) for i in list_dup])   # extract d_up
list_ddown = [l.split()[6] for l in open(file,'rb')]
list_ddown.pop(0)
d_down = np.array([float(i) for i in list_ddown])   # extract d_down

y1_d = d_up[emask]
y2_d = d_down[emask]

dbc_up   = simps(y1_d*x, x) / simps(y1_d, x)
dbc_down = simps(y2_d*x, x) / simps(y2_d, x)
dbc = []
dbc.append(dbc_up)
dbc.append(dbc_down)

###### calculation of p_z-band center #######

## same energy set ##

dos_splitted = np.loadtxt('atoms_dos.dat')
dos_s = np.zeros([epoints,column_num])
for j in range(int(atom_ini),int(atom_final)+1):
    dos_s_add = np.vsplit(atoms_dos,atoms_number)[j-1].astype(float)
    dos_s_add = np.add(dos_s_add,Ef0) # corrected energy by efermi
    dos_s_add = np.multiply(dos_s_add,ones)
    dos_s = np.add(dos_s,dos_s_add)    # vertical stack

dos_s[:,0] = energy
df_s = pd.DataFrame(dos_s, columns=all_names)
df_s.set_index('energy(eV)', inplace=True)
df_s.to_csv('select_dos_split.dat', sep=' ', float_format='%-12.10f')  # projected to each atomic orbital orientation

list_pzup = [l.split()[5] for l in open('select_dos_split.dat','rb')]
list_pzup.pop(0)
pz_up = np.array([float(i) for i in list_pzup])   # extract p_up
list_pzdown = [l.split()[6] for l in open('select_dos_split.dat','rb')]
list_pzdown.pop(0)
pz_down = np.array([float(i) for i in list_pzdown])   # extract p_down

y1_pz = pz_up[emask]
y2_pz = pz_down[emask]

pzbc_up   = simps(y1_pz*x, x) / simps(y1_pz, x)
pzbc_down = simps(y2_pz*x, x) / simps(y2_pz, x)
pzbc = []
pzbc.append(pzbc_up)
pzbc.append(pzbc_down)


print 'the slected atoms:', range(int(atom_ini), int(atom_final)+1)
print 'dbc_up(eV), dbc_down(eV)'
print dbc
print 'p-BC_up(eV), p-BC_down(eV)'
print pbc
print 'p_z-BC_up(eV), p_z-BC_down(eV)'
print pzbc

os.system('rm atoms_dos* dos*')
##########   timing   #############
stop=time.time()
print("running time:     " + str(stop-start) + " seconds")
terminal_time = datetime.datetime.now()
print "Terminal time:    " + terminal_time.strftime('%Y.%m.%d-%H:%M:%S')   #strftime可以自定义时间的输出格式
print 'splitted dos and d-band center have been obtained'
print '********** splitted dos and d-band center from Yaqiong Su Eindhoven **********'
