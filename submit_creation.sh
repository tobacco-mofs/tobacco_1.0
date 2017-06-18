#!/bin/bash
#PBS -N create_MOFs_BCU
#PBS -l nodes=1:ppn=1
#PBS -r n
#PBS -j oe
#PBS -V
cd $PBS_O_WORKDIR


python $PBS_O_WORKDIR/main_auto.py

