mkdir DEanalysis/counts
for a in mapping/*/*.tsv; do name=$(dirname $a | cut -f2 -d '/'); sed -e "s/est_counts/$name/g" -e "s/tpm/$name\_tpm/g" $a > DEanalysis/counts/$(echo $a | tr / _ | sed 's/mapping_//g'); done

cut -f14,20,21 DEanalysis/Chenopodium.mergedQuartets.one-file-to-rule-them-all.txt | grep "Cqui" | sed 's/Cquinoa_//g'  > DEanalysis/gene.pairs
#Targeting	Paternal	Maternal

cut -f1 DEanalysis/gene.pairs | sort | uniq | while read line; do grep $line DEanalysis/gene.pairs | cut -f2,3 | awk '{OFS="\t"} {print $1,$2,$1"_"$2}' | sed 's/\t/\n/g' | sort > DEanalysis/$line.list; done

