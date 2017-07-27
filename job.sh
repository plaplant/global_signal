#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -l h_vmem=16G
#$ -j y
#$ -N global_signal
#$ -o log
#$ -l all
#$ -pe omp 16
#$ -m be
#$ -M plaplant@sas.upenn.edu

# Compiler vairables
source /opt/intel/composerxe-2011.4.191/bin/compilervars.sh intel64
which ifort

# Set environment variables
export KMP_LIBRARY=turnaround
export KMP_SCHEDULE=static,balanced
export KMP_STACKSIZE=128m

# Compile and run job
make clean
make
./calculate_xi.x
