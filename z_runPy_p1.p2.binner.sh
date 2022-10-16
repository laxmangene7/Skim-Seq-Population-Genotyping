#!/bin/bash -l
#SBATCH --job-name=Run_Bin_Averager
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=1-00:00:00   # Use the form DD-HH:MM:SS
#SBATCH --mem-per-cpu=15G    # Memory per core, use --mem= for memory per node
#SBATCH --output="%x_%j.o"
#SBATCH --error="%x_%j.e"


#module load Python
#source ~/virtualenvs/Python_Env/bin/activate
#export PYTHONDONTWRITEBYTECODE=1

module load SciPy-bundle/2020.03-foss-2020a-Python-3.8.2

export PYTHONDONTWRITEBYTECODE=1
mpirun python  z_p1.p2.per.bin.py
