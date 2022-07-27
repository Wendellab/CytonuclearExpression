#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH --cpus-per-task=150
#SBATCH --mem=300G 
#SBATCH --time=7-02:30:00
#SBATCH --output=job.arabidopsisRM.%J.out 
#SBATCH --error=job.arabidopsisRM.%J.err
#SBATCH -J "arabidopsisRM" 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu 


module load repeatmasker/4.0.7-openmpi3-vruhfro

cd prepareSource

mkdir RM.dir
RepeatMasker -pa 75 -s -nolow -lib arabidopsis_organelles.fa -no_is -nocut -noisy -dir RM.dir -html -gff nuclear.transcripts.fa
mv RM.dir/nuclear.transcripts.fa.masked ./nuclear.transcripts.masked.fa

ml py-biopython/1.65-py2-fgq2vqv
python sequenceCleaner.py nuclear.transcripts.masked.fa 75 70
cat clean_nuclear.transcripts.masked.fa arabidopsis.chloroplast.cds.20210610.fa arabidopsis.mitochondrial.cds.20210610.fa > arabidopsisref.fa

ml  


