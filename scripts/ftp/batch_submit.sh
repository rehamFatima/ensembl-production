#!/bin/bash

#SBATCH -t 00:05:00

#SBATCH --mem=20

#SBATCH -o /homes/reham/CARSScripts/SLURM_Os/%j.out --error /homes/reham/CARSScripts/SLURM_Os/%j.err

for alpha in {A..Z}
do
  sbatch -t 7-00:00:00 --mem=3G --wrap="perl paralpha_ftp_copy_new_struct.pl -alpha $alpha" -o /homes/reham/CARSScripts/SLURM_Os/%j.out --error /homes/reham/CARSScripts/SLURM_Os/%j.err
  echo "Submitted chromosome $alpha"
done
