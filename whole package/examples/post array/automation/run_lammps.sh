#!/bin/bash
#SBATCH -J POSTA2B8
#SBATCH -o posta2b8.out
#SBATCH -n 32
#SBATCH -N 1
#SBATCH -p extended-mem


  mpirun -np 32  /pool/hejin/bin/lmp_minimal_31May17 -in in.asymmetric

