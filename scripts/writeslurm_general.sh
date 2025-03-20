#!/bin/bash
#
# Writes a slurm script to be passed into Kennedy
#
read -p 'Username: ' username
read -p 'Jobname:  '  jobname
read -p 'Infile:   '  infile
read -p 'Start command: ' command 
read -p 'Debug queue (y/n)? ' -e -i 'n' yn

echo '#!/bin/bash -l'
echo '#SBATCH --job-name '$jobname
echo '#SBATCH --nodes=1'
case $yn in
   [Yy]* ) echo '#SBATCH --partition small-short';;
   [Nn]* ) echo '#SBATCH --partition large-long';;
   * ) echo '#SBATCH --partition large-long';;
esac
echo '#SBATCH --ntasks=1'
case $yn in 
   [Yy]* ) echo '#SBATCH --cpus-per-task=20';;
   [Nn]* ) echo '#SBATCH --cpus-per-task=50';;
   * ) echo '#SBATCH --cpus-per-task=50';;
esac
case $yn in
   [Yy]* ) echo '#SBATCH --time=20:00:00';;
   [Nn]* ) echo '#SBATCH --time=168:00:00';;
   * ) echo '#SBATCH --time=168:00:00';;
esac
echo '#SBATCH --mail-type=BEGIN,END,FAIL'
echo '#SBATCH --mail-user='$username'@st-andrews.ac.uk'
echo '#SBATCH --mem=80G'

echo 'export OMP_SCHEDULE="dynamic"'
case $yn in
   [Yy]* ) echo 'export OMP_NUM_THREADS=48';;
   [Nn]* ) echo 'export OMP_NUM_THREADS=96';;
   * ) echo 'export OMP_NUM_THREADS=96';;
esac
echo 'export KMP_STACKSIZE=128M'
echo 'ulimit -s unlimited'

echo 'module load openmpi/5.0.5'
echo 'module load gcclibs/11.4.1'
echo 'module load ucx/1.16.0'

echo 'echo $SLURM_NPROCS'
echo 'echo "HOSTNAME = $HOSTNAME"'
echo 'echo "HOSTTYPE = $HOSTTYPE"'
echo 'echo Time is `date`'
echo 'echo Directory is `pwd`'

echo 'echo "starting run"'
echo 'export outfile="'$jobname'.log"'
echo 'echo "writing output to $outfile"'
echo $command' '$infile' >& $outfile'
