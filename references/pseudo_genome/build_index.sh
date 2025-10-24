bedtools getfasta -fo transcripts.fasta -fi pseudo_genome.fasta -bed RNA_index.gff
./generateDecoyTranscriptome.sh --keepDuplicates -m /home/molgenuser/MashMap-3.1.3/build/bin/mashmap -m /home/molgenuser/MashMap-3.1.3/build/bin/mashmap -g pseudo_genome.fasta -t transcripts.fasta -a RNA_index.gff -o decoy
salmon index -t decoy/gentrome.fa -i transcripts_index --decoys decoy/decoys.txt -k 13
