#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH --cpus-per-task=150
#SBATCH --mem=300G 
#SBATCH --time=7-02:30:00
#SBATCH --output=job.arachisRM.%J.out 
#SBATCH --error=job.arachisRM.%J.err
#SBATCH -J "arachisRM" 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu 


module load repeatmasker/4.0.7-openmpi3-vruhfro

cd prepareSource

mkdir -p RM.dir
RepeatMasker -pa 75 -s -nolow -lib arachis_organelles.fa -no_is -nocut -noisy -dir RM.dir -html -gff nuclear.transcripts.fa
mv RM.dir/nuclear.transcripts.fa.masked ./nuclear.transcripts.masked.fa

ml py-biopython/1.65-py2-fgq2vqv
python sequenceCleaner.py nuclear.transcripts.masked.fa 75 70
cat clean_nuclear.transcripts.masked.fa arachis.chloroplast.cds.20210610.fa  arachis.mitochondrial.cds.20210610.fa > arachisref.fa

ml  


