# Job Name
#$ -N calculix_test
# Use current working directory
#$ -cwd

# Parallel Environment request.  Set your number of processors here
#$ -pe mpich 1

# Run job through bash shell
#$ -S /bin/bash

source /etc/profile.d/modules.sh
module load gcc/7.2.0

echo Got $NSLOTS processors.
echo Machines:
echo JOBNAME: $JOBNAME
echo JOBN0: $JOBN0
echo JOB_ID: ${JOB_ID}
cat $TMPDIR/machines

cwd=$(pwd)
echo pwd: $cwd

export OMP_NUM_THREADS=1

SECONDS=0

