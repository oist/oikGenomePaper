#!/usr/bin/bash

#SBATCH --mail-user="aleksandra.bliznina2@oist.jp"
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH -t 10:00:00
#SBATCH --mem 50G
#SBATCH -c 8

module use /apps/unit/LuscombeU/.modulefiles/
module load augustus/3.3.1_oiko

DIR=/work/LuscombeU/aleksandrabliznina/genome_annotation/augustus/3.3

augustus --species=oikopleura_dioica_okinawa_1 --AUGUSTUS_CONFIG_PATH=$DIR/config --alternatives-from-sampling=true --minexonintronprob=0.08 --minmeanexonintronprob=0.4 --maxtracks=3 --progress=true --softmasking=on --gff3=on i69-4.juicer.fa.masked > augustus.abinitio-3.gff
