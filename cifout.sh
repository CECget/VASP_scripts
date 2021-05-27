#!/bin/sh

#####################################################
### Created on Tues Dec 8 2020
### for creating cif files automatically IN ALL CASES
### @author: QYY
#####################################################

filename=`awk '/-J/{print $3}' job_submit`
##if [ -z "`grep "Selective Dynamics" POSCAR`" ];then
##	key_loc=7
##else
##	key_loc=8
##fi

keywords=`awk 'NR==6 {print $1}' POSCAR`

if [ -z "`tail -5 output 2>/dev/null |grep "reached required accuracy" 2>/dev/null`" ];then
	echo "Calculation NOT converged! Continue to get cif file for confirming."	
	if [ ! -f "$filename-cfm.cif" ]; then
		if [ "$keywords" -gt 0 ] 2>/dev/null ; then
			vasp-pos-to-cif POSCAR > $filename-cfm.cif
		else  
     		sed '6d' POSCAR > result 
            vasp-pos-to-cif result > $filename-cfm.cif
		fi
		rm result
		echo "cfm for confirming"
		du -h $filename-cfm.cif
	else
		keywords=`awk 'NR==6 {print $1}' CONTCAR`
		if [ "$keywords" -gt 0 ] 2>/dev/null ; then
			vasp-pos-to-cif CONTCAR > $filename-mid.cif
		else	
			sed '6d' CONTCAR > result
			vasp-pos-to-cif result > $filename-mid.cif
		fi
		rm result
		echo "mid for checking"
		du -h $filename-mid.cif
	fi
	echo "Done!"
	exit 0
else 
	echo "Calculation is CONVERGED! Continue to get cif file as final structure."
	break
fi 
 
keywords=`awk 'NR==6 {print $1}' CONTCAR`

if [ "$keywords" -gt 0 ] 2>/dev/null ; then
	vasp-pos-to-cif CONTCAR > $filename-end.cif
else	
	sed '6d' CONTCAR > result
	vasp-pos-to-cif result > $filename-end.cif
	rm result
fi
echo "end for saving"
du -h $filename-end.cif
echo "Done!"	
