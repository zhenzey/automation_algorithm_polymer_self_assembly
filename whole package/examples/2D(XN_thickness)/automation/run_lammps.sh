#!/bin/bash
#SBATCH -J XN_Z	
#SBATCH -o xn_z.out
#SBATCH -n 32
#SBATCH -N 1
#SBATCH -p extended-mem


  mpirun -np 32  /pool/hejin/bin/lmp_minimal_31May17 -in in.asymmetric

