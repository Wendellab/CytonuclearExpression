cd prepareSource


# chloroplast and mitochondria genomes were downloaded from GenBank 
# NCBI accession numbers are given in the original file names
# gene content in each was modified by Dan Sloan
# file names are now gossypium.chloroplast.cds.20210428, gossypium.mitochondrial.cds.20210428

# make repeatmasker database

cat arachis.chloroplast.cds.20210610.fa  arachis.mitochondrial.cds.20210610.fa  NC_037358cp.fasta > arachis_organelles.fa

mv arahy.Tifrunner.gnm2.ann1.4K0L.cds_primaryTranscript.fna nuclear.transcripts.fa