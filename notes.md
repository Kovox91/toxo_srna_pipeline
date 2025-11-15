# LSU/SSU Zuordnung der veränderten Fragmente (wirkt mir wie hauptsächlich SSU?)
- Large Features (LSUFG, LSUDE. LSUA) have no reads because they are to large (fractionOverlapFeature from featureCounts) **General Issue**

- Still need to quantify?
	- percentages
- Most recent assignments!
	- Soldatifari
	

# bei denen mit mehreren subfragmenten scheinen mir nicht immer alle gleich (stark) verändert? Können wir da einen Überblick bekommen bzw ein gff file das man in IGV öffnen kann wäre cool ;)
- Should work with (modified) .gtf file, Try on windows machine.
- Maybe tabular data? Discuss! |> Yes!


# ich fände es interessant die mito reads mal genauer anzuschauen, die nicht auf sRNA features mappen (bei HPR4 waren da ja ziemlich viele oder?). Verändern die sich im ganz allgemeinen und vielleicht auch spezifisch reads die die 3 mRNAs mappen.. laut Northern ist cox3 an Tag 1 schon stark reduziert, wäre interessant ob man mit den Seq Daten zumindest einen Trend bestätigen kann, auch wenn es nicht die präferierte Methode dafür ist
- 5 - 11% of reads are unassigned in the helicase dataset, 3-5% in HPR4.
- Raw reads follow the same downwards trend as all the rRNAs 
- Can featureCount with mRNA features, would also 'solve' the large feature issue -> No Overlap filtering! **General Issue**

# und natürlich: wie verändern sich die antisense RNAs?
- should be relatively simple by featureCounts with strandednes set to 2? 
- Quantification? What would you like to determine? |> Ratios & Normalized counts

# subunits percentages
# Calculate % mito only mapped reads of total lib size

# Pipeline
- Restructure output folder
- Separate .bam files by read direction save as additional files. (IGV only)
	-> bedtools!
	for bam in out/filtered/*_mito_filtered.bam; do
	base=$(basename "${bam%.bam}")   # e.g. sample_mito_filtered
	echo "Processing $bam …"
    	bedtools genomecov
    	         -ibam "$bam"
 	         -strand + -bg
                 > "out/filtered/strand_cov/${base}_plus.bedgraph"
    	bedtools genomecov
    	         -ibam "$bam"
    	         -strand - -bg
    	         > "out/filtered/strand_cov/${base}_minus.bedgraph"; done

- remove -O, add --largestOverlap
- Add second pass, only for long features (filter in R later), --fracOverlap .9



