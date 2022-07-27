ml gmap-gsnap
gmap_build -d rbcs -D . all_rbcS.txt 
gmap -D . -d rbcs -t 20 --cross-species -f 3 --nofails clean_nuclear.transcripts.masked.fa > rbcs.gmap
cut -f9 rbcs.gmap | sed -e '/^#/d' -e 's/ID=//g' -e 's/.path.*//g' | sort | uniq > rbcS.candidates

ml samtools
xargs samtools faidx clean_nuclear.transcripts.masked.fa < rbcS.candidates >> rbcS.candidates.fasta
