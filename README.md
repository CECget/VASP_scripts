# VASP_scripts
Useful scripts for VASP (Vienna Ab-initio Simulation Package)
# vfcals.py
Function: Calculating the correction terms of free energy G (TS, CvT, ZPE) after vabrational frequency calculation  
Usage: $ vfcals.py  
Check output: $ vi tsresult  
# getdpzbc.py / getdpzbc_p3.py  
Function: Calculating the d/pz band center of certain atoms from DOSCAR, and get the corresponding PDOS for ploting    
Usage: 
1. run single point calculation with high KPOINTS mesh to obtain the DOSCAR   
2. $ mkdir newfile  
   $ cd newfile  
   $ cp ../DOSCAR .   
   $ getdpzbc.py x1 x2 (here use x1, x2 to select the atoms' range you want to calculate,x1 represent the atom index you want to begin with, x2 represent that to end with)
   
Check output: the results will print on screen directly
# magic.sh
Function: Manage all the subfolders of a selected pathway ("path/") in one time  
Usage: $ magic.sh path/
DIY according to your personal needs, in the present script, it will run cifout.sh to output a .cif file in every subfolder.  

# Waiting for update ...
