cd prepareSource


# chloroplast and mitochondria genomes were downloaded from GenBank 
# NCBI accession numbers are given in the original file names
# gene content in each was modified by Dan Sloan
# file names are now arabidopsis.chloroplast.cds.20210610.fa, arabidopsis.mitochondrial.cds.20210610.fa

# make repeatmasker database

cat arabidopsis.chloroplast.cds.20210610.fa arabidopsis.mitochondrial.cds.20210610.fa NC_000932cp.fasta NC_037304mt.fasta > arabidopsis_organelles.fa

cat Asuecica_Aa.cds.primaryTranscriptOnly.fasta Asuecica_At.cds.primaryTranscriptOnly.fasta > nuclear.transcripts.fa