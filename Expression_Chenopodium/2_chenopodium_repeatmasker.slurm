#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH --cpus-per-task=150
#SBATCH --mem=300G 
#SBATCH --time=7-02:30:00
#SBATCH --output=job.cottonRM.%J.out 
#SBATCH --error=job.cottonRM.%J.err
#SBATCH -J "cottonRM" 
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu 


module load repeatmasker/4.0.7-openmpi3-vruhfro

cd prepareSource

#mkdir RM.dir
#RepeatMasker -pa 75 -s -nolow -lib chenopodium_organelles.fa -no_is -nocut -noisy -dir RM.dir -html -gff Cquinoa.cds.primary.fasta
mv RM.dir/Cquinoa.cds.primary.fasta.masked ./Cquinoa.cds.primary.masked.fa

ml py-biopython/1.65-py2-fgq2vqv
python sequenceCleaner.py Cquinoa.cds.primary.masked.fa 75 70
cat clean_Cquinoa.cds.primary.masked.fa chenopodium.chloroplast.cds.20210428.fa chenopodium.mitochondrial.cds.20210428.fa > chenopodiumref.fa

ml  


