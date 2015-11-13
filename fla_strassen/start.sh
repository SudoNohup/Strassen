tacc_cpufreq 2700000
#module load vtune
#module load intel/14.0.1.106
#export OMP_PROC_BIND=TRUE
#export KMP_AFFINITY=scatter,granularity=fine,verbose
#export KMP_AFFINITY=compact,verbose
export KMP_AFFINITY=compact,verbose
export MKL_NUM_THREADS=16
export OMP_NUM_THREADS=16
