#!/bin/bash
#SBATCH -J A2B8singlelayer
#SBATCH -o a2b8singlelayer.out
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p regular-cpu
#SBATCH --reservation=hejin_1

  mpirun -np 8  /pool/hejin/bin/lmp_minimal_31May17 -in in.asymmetric

