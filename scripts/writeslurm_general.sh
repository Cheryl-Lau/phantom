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
echo '#SBATCH -J '$jobname
echo '#SBATCH --nodes=1'
case $yn in
   [Yy]* ) echo '#SBATCH -p debug';;
   [Nn]* ) echo '#SBATCH -p singlenode';;
   * ) echo '#SBATCH -p singlenode';;
esac
echo '#SBATCH --ntasks=1'
echo '#SBATCH --cpus-per-task=32'
case $yn in
   [Yy]* ) echo '#SBATCH --time=2:00:00';;
   [Nn]* ) echo '#SBATCH --time=168:00:00';;
   * ) echo '#SBATCH --time=168:00:00';;
esac
echo '#SBATCH --output='$infile'.qout'
echo '#SBATCH --error='$infile'.err'
echo '#SBATCH --mail-type=BEGIN,END,FAIL'
echo '#SBATCH --mail-user='$username'@st-andrews.ac.uk'
echo '#SBATCH --mem=80G'

echo 'export OMP_SCHEDULE="dynamic"'
echo 'export OMP_NUM_THREADS=32'
echo 'export KMP_STACKSIZE=128M'
echo 'ulimit -s unlimited'

echo 'echo $SLURM_NPROCS'
echo 'echo "HOSTNAME = $HOSTNAME"'
echo 'echo "HOSTTYPE = $HOSTTYPE"'
echo 'echo Time is `date`'
echo 'echo Directory is `pwd`'

echo 'echo "starting run"'
echo 'export outfile='$jobname'.log'
echo 'echo "writing output to $outfile"'
echo $command' '$infile' >& $outfile'
