#!/bin/bash

##########for correction terms calculation by QYY

mkdir fqc
cd fqc
cp  ../* .  2>/dev/null
cp CONTCAR POSCAR
rm CONTCAR DOSCAR EIGENVAL IBZKPT OUTCAR OSZICAR PCDAT WAVECAR XDATCAR output slurm-* vasprun.xml *.cif PROCAR 2>/dev/null
cp ~/vaspjobs/yanyang/bin/{freeze,vfcals.py} .



##########改写INCAR用于频率计算

######鲁棒性较差写法
#sed -n '/IBRION/p' INCAR | sed -i 's/2/5/g' INCAR
#sed -n '/POTIM/p' INCAR | sed -i 's/0.10/0.02/g' INCAR
#sed -i "/POTIM/a\\NFREE = 2" INCAR

######改良写法
Index_IBRION=`grep -n IBRION INCAR | awk -F ":" '{print$ 1}'`
sed -i "${Index_IBRION}c IBRION = 5 " INCAR
Index_POTIM=`grep -n POTIM INCAR | awk -F ":" '{print$ 1}'`
sed -i "${Index_POTIM}c POTIM = 0.02 " INCAR
sed -i "${Index_POTIM}a NFREE = 2 " INCAR



##########检查并删除POSCAR坐标锁定关键字
grep "Selective" POSCAR >/dev/null
if [ $? -eq 0 ]; then
	sed -i '/Selective/d' POSCAR
fi



##########冻结所有基底原子
./freeze POSCAR > pos ###冻结标准基于z轴坐标，在上方freeze脚本的目录中设置
mv pos POSCAR

rm fqcCal.sh freeze

echo "Done!"





