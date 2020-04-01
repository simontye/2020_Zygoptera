#!/bin/bash

#PBS -N test_coi_sh
#PBS -q debug12core
#PBS -j oe
#PBS -m abe
#PBS -M simontye@uark.edu
#PBS -o test_coi_sh.out
#PBS -l nodes=1:ppn=12
#PBS -l walltime=00:00:30:00

### Change directory
cd /storage/simontye/2020_Zygoptera

### Load modules
module load muscle/3.8.31

### Run muscle
muscle -in test_coi_20200328.fasta -out test_coi_20200328_sh.phy -seqtype protein

### Run iqtree
#iqtree -s test_coi_20200328.phy
#iqtree -s odonata_coi_20200328.phy
