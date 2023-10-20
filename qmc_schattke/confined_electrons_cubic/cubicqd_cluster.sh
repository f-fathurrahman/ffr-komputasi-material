#!/bin/sh

#PBS -q serial
#PBS -l mem=500Mb
#PBS -l cput=400:00:00
#PBS -N Jmcc05_sm4
#
JOB="cubicqd_cluster"
WD=/local_scratch/schattke_$JOB
SD=$HOME/table/work/$JOB
HD=$HOME/table

if test ! -e $WD
then
echo "because $WD does not exist, it is generated"
mkdir $WD
fi

if test ! -e $SD
then
echo "because $SD does not exist, it is generated"
mkdir $SD
fi


cd $WD
cp $HD/*.f .

f95 M_variables_cubicqd_cluster.f M_random_cubicqd_cluster.f M_orbital_cubicqd_cluster.f M_determinant_cubicqd_cluster.f M_jastrow_cubicqd_cluster.f M_observables_cubicqd_cluster.f M_output_cubicqd_cluster.f cubicqd_cluster.f -o cubicqd_cluster.x

cp $HD/AKF.DAT .

echo 'Work begins at :' 
date -u 

AB="  "
A="0.0  0.4  0.8  1.2  1.6  2.0  2.4  2.8  3.2"
B="0.0  1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0  9.0  10.0  11.0  12.0  13.0  14.0  15.0  20.0"
CJAS=0.50
C="_"
n=0
for i in $B 
do
  for k in $A
  do
  cat > aaa$i$C$k << EZ!
time ./Lisolid_$JOB.x << [EOF]
2000  64000  200
$i$AB$k
3.0
$CJAS
[EOF]
EZ!
  cat aaa$i$C$k > Li_$i$C$k.job
  chmod 700 Li_$i$C$k.job
  ./Li_$i$C$k.job > Li_out_$JOB$C$i$C$k.log
  done
done

cp $WD/*.f $SD
cp $WD/*.sh $SD
cp $WD/*.DAT $SD
cp $WD/*.log $SD

echo 'Work ends at:' 
date -u 


