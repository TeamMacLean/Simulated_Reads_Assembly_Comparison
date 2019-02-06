#!/bin/bash

# run assembly
source abyss-2.1.2
kmer=$1
coverage=$2
kc=$3
R1=$4
R2=$5

if test -e assembly_completed; then
	echo "Assembly completed"
else
	ABYSS -k $kmer -o contigs.fasta --trim-masked --no-SS --coverage=${coverage}.0 --kc=$kc  $R1 $R2 && touch assembly_completed

fi


# filter the contigs greater or equal to 200 bps

cat contigs.fasta | paste - - | awk '{if($2 >= 200){print $1,$2,$3"\n"$4}}' > contigs.fasta.filtered && mv contigs.fasta.filtered contigs.fasta

#run blast

source blast+-2.3.0

top=1
mismatches=1
gaps=1

if test -e blast_completed; then
	echo "Blast completed"
else
	if test -e contigs.fasta; then
		if test ! -e  ~/workarea/assembly_comparison/data/PGSC_DM_v3_scaffolds.fasta.nhr; then makeblastn -in ~/workarea/assembly_comparison/data/PGSC_DM_v3_scaffolds.fasta  -out ~/workarea/assembly_comparison/data/PGSC_DM_v3_scaffolds.fasta -dbtype nucl; fi

		blastn -db ~/workarea/assembly_comparison/data/PGSC_DM_v3_scaffolds.fasta  -query contigs.fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore' | awk -v mismatches="$mismatches" -v gaps="$gaps" -F "\t" '{if(mismatches <= $5 && gaps <= $6){print}}' | sort -k1,1 -k4nr -k3nr | awk -v top="$top" 'BEGIN{query=""; counter=1}{if(query!=$1){counter=1; print; query=$1;} else if(query==$1 && counter<top){print $0; counter+=1}}' > blast_result.txt && touch blast_completed
	else
		echo "assembly file does not exist. blastn cannot be executed."
	fi
fi

# blast baits to contigs.fasta

if test -e blast_contigs_to_baits_completed; then
	echo "Bait blat completed"
else
	if test -e contigs.fasta; then
		makeblastdb -in contigs.fasta -out contigs.fasta -dbtype nucl
		blastn -db contigs.fasta -query ~/workarea/assembly_comparison/data/renseq-bait-library-sequences.fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore' | awk -v mismatches="$mismatches" -v gaps="$gaps" -F "\t" '{if(mismatches <= $5 && gaps <= $6){print}}' | sort -k1,1 -k4nr -k3nr | awk -v top="$top" 'BEGIN{query=""; counter=1}{if(query!=$1){counter=1; print; query=$1;} else if(query==$1 && counter<top){print $0; counter+=1}}' > bait_blast_result.txt && touch blast_contigs_to_baits_completed && rm contigs.fasta.n*
	fi
fi

#run quast
source quast-5.0.2
source python-3.7.1

if test -e quast_completed; then
	echo "Quast completed"
else
	if test -e contigs.fasta; then
		quast --output-dir ./quast -r ~/workarea/assembly_comparison/data/PGSC_DM_v3_scaffolds.fasta --min-contig 200 --threads 4 --eukaryote --min-alignment 200 contigs.fasta && touch quast_completed
	else
		echo "assembly file does not exist. quast cannot be executed."
	fi
fi

#run transrate

source transrate-1.0.3
source blast+-2.2.28

if test -e transrate_completed; then
	echo "Transrate completed"
else
	if test -e contigs.fasta; then
		transrate --reference ~/workarea/assembly_comparison/data/PGSC_DM_v3_scaffolds.fasta --assembly contigs.fa --threads 4 --output ./transrate && touch transrate_completed
	else
		echo "assembly file does not exist. transrate cannot be executed."
	fi
fi
