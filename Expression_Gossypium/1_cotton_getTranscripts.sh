module load samtools

cd prepareSource
curl -O ftp://ftp.bioinfo.wsu.edu/species/Gossypium_hirsutum/UTX-TM1_v2.1/genes/Ghirsutum_527_v2.1.transcript_primaryTranscriptOnly.fa.gz
zcat Ghirsutum_527_v2.1.transcript_primaryTranscriptOnly.fa.gz | sed  -e 's/^.*ID=/>/g' -e 's/[.][1-9][0-9][.]v2.*//g' -e 's/[.][1-9][.]v2.*//g' > nuclear.all.transcripts.fa
grep "Gohir.[AD]"  nuclear.all.transcripts.fa | sed 's/>//g' > Gohir.names
xargs samtools faidx nuclear.all.transcripts.fa < Gohir.names >> nuclear.transcripts.fa
rm Ghirsutum_527_v2.1.transcript_primaryTranscriptOnly.fa.gz nuclear.all.transcripts.fa Gohir.names


# chloroplast and mitochondria genomes were downloaded from GenBank 
# NCBI accession numbers are given in the original file names
# gene content in each was modified by Dan Sloan
# file names are now gossypium.chloroplast.cds.20210428, gossypium.mitochondrial.cds.20210428

# make repeatmasker database

cat gossypium.chloroplast.cds.20210428.fa gossypium.mitochondrial.cds.20210428.fa <(echo) NC_007944cp.fasta JX065074mt.fasta > cotton_organelles.fa
