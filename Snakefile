# Define the deafult conda environment
conda: "environment.yaml"   

# Import necessary modules
import pandas as pd

# Load configuration
configfile: "config.yaml"

# Load sample information from CSV file
samples_df = pd.read_csv(config['samples_csv'])

# Create a list of sample names
sample_names = list(samples_df['sample'])

# Define rule all to specify final output files
rule all:
    input:
        # === Trimming & QC ===
        trimmed_R1 = expand("out/trimmed/{sample}_R1_001_trimmed.fq.gz", sample=sample_names),
        QC_report = expand("out/fastqc/{sample}_R1_001_trimmed_fastqc.html", sample=sample_names),        

        # === Mapping ===
        mito_sam_raw = expand("out/mapped/{sample}_mito_mapped.sam",  sample=sample_names),
        decoy_sam_raw = expand("out/mapped/{sample}_decoy_mapped.sam", sample=sample_names),

        # === Read Filtering ===
        mito_filtered = expand("out/filtered/{sample}_mito_filtered.bam", sample=sample_names),
        decoy_filtered = expand("out/filtered/{sample}_decoy_filtered.bam", sample=sample_names),

        # === Sorting and Indexing ===
        mito_bam = expand("out/filtered/{sample}_mito_filtered.bam", sample=sample_names),
        mito_bai = expand("out/filtered/{sample}_mito_filtered.bam.bai", sample=sample_names),

        decoy_bam = expand("out/filtered/{sample}_decoy_filtered.bam", sample=sample_names),
        decoy_bai = expand("out/filtered/{sample}_decoy_filtered.bam.bai", sample=sample_names),

        # === Counting ===
        counts_mito = expand("out/counts/{sample}_mito_featureCounts.txt", sample=sample_names),
        counts_decoy_unique = expand("out/counts/{sample}_decoy_unique_featureCounts.txt", sample=sample_names),
        counts_decoy_multi = expand("out/counts/{sample}_decoy_multi_featureCounts.txt", sample=sample_names)

# Create the trim and filter rule
rule cutadapt_trim_all:
    input:
        R1="in/{sample}_R1_001.fastq.gz"
    output:
        R1_trimmed="out/trimmed/{sample}_R1_001_trimmed.fq.gz"
    threads: 8
    shell:
        r"""
        cutadapt -j {threads} -q 20 -m 15 --trim-n \
          -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
          --poly-a \
          --report=full \
          -o {output.R1_trimmed} {input.R1} \
          > out/trimmed/{wildcards.sample}_cutadapt_report.txt
        """

rule fastqc:
    input:
        "out/trimmed/{sample}_R1_001_trimmed.fq.gz"
    output:
        html = "out/fastqc/{sample}_R1_001_trimmed_fastqc.html"
    threads: 8
    shell:
        r"""
        fastqc --threads {threads} \
               --outdir out/fastqc \
               {input}
        """


# Map against decoy
rule map_decoy:
    input:  R1_trimmed = "out/trimmed/{sample}_R1_001_trimmed.fq.gz"
    output: decoy_mapped = temp("out/mapped/{sample}_decoy_mapped.sam")
    params:
        index = "references/decoy/ToxoDB-68_TgondiiRH88_Genome", # ebwt basename
    threads: 8
    shell: 
        r"""
        mkdir -p out/mapped
        zcat {input.R1_trimmed} | \
        bowtie -q -q -v 1 -k 4 --best --strata -p {threads} \
            -S {params.index} - > {output.decoy_mapped}
        """



# Map against pseudo genome
rule genome_mapping:
    input:
        R1_trimmed = "out/trimmed/{sample}_R1_001_trimmed.fq.gz"
    output:
        genome_mapped = temp("out/mapped/{sample}_mito_mapped.sam")
    params:
        index = "references/pseudo_genome/pseudo_genome"   # bowtie1 index basename
    threads: 8
    shell:
        r"""
        mkdir -p out/mapped
        zcat {input.R1_trimmed} | \
        bowtie -q -v 1 -a --best --strata -p {threads} \
            -S {params.index} - > {output.genome_mapped}
        """

rule sort_by_name_mito:
    input:  "out/mapped/{sample}_mito_mapped.sam"
    output: "out/mapped/{sample}_mito_mapped_sorted.bam"
    threads: 8
    shell:
        # drop unmapped (-F 4), keep multi-mappers (DO NOT add -F 0x100/0x800)
        "samtools view -@{threads} -b -F 4 {input} | "
        "samtools sort -n -@{threads} -o {output} -"

rule sort_by_name_decoy:
    input:  "out/mapped/{sample}_decoy_mapped.sam"
    output: "out/mapped/{sample}_decoy_mapped_sorted.bam"
    threads: 8
    shell:
        "samtools view -@{threads} -b -F 4 {input} | "
        "samtools sort -n -@{threads} -o {output} -"

rule classify_reads:
    input:
        mito  = "out/mapped/{sample}_mito_mapped_sorted.bam",   # coord-sorted inputs OK
        decoy = "out/mapped/{sample}_decoy_mapped_sorted.bam"
    output:
        # write UNSORTED temp BAMs from the classifier
        mito_uns  = temp("out/filtered/{sample}_mito_filtered.unsorted.bam"),
        decoy_uns = temp("out/filtered/{sample}_decoy_filtered.unsorted.bam")
    threads: 8
    shell:
        r"""
        python scripts/Python/classify_reads.py \
          --mito {input.mito} \
          --decoy {input.decoy} \
          --mito-out {output.mito_uns} \
          --decoy-out {output.decoy_uns}
        """

rule sort_index_mito:
    input:
        "out/filtered/{sample}_mito_filtered.unsorted.bam"
    output:
        bam = "out/filtered/{sample}_mito_filtered.bam",
        bai = "out/filtered/{sample}_mito_filtered.bam.bai"
    threads: 8
    shell:
        r"""
        samtools sort -@ {threads} -o {output.bam} {input}
        samtools index {output.bam}
        """

rule sort_index_decoy:
    input:
        "out/filtered/{sample}_decoy_filtered.unsorted.bam"
    output:
        bam = "out/filtered/{sample}_decoy_filtered.bam",
        bai = "out/filtered/{sample}_decoy_filtered.bam.bai"
    threads: 8
    shell:
        r"""
        samtools sort -@{threads} -o {output.bam} {input}
        samtools index {output.bam}
        """

rule featurecounts_mito:
    input:
        bam = "out/filtered/{sample}_mito_filtered.bam"
    output:
        tsv = "out/counts/{sample}_mito_featureCounts.txt",
        sum = "out/counts/{sample}_mito_featureCounts.txt.summary"
    threads: 8
    params:
        stranded = 1, feature = "exon", attr = "gene_id"
    shell:
        r"""
        featureCounts -T {threads} \
          -a "references/pseudo_genome/RNA_index.gtf" -t {params.feature} -g {params.attr} \
          -s {params.stranded} \
          -M -O --fraction \
          -o {output.tsv} {input.bam}
        """

rule featurecounts_decoy_unique:
    input:
        bam = "out/filtered/{sample}_decoy_filtered.bam"
    output:
        tsv = "out/counts/{sample}_decoy_unique_featureCounts.txt",
        sum = "out/counts/{sample}_decoy_unique_featureCounts.txt.summary"
    threads: 8
    params:
        stranded = 1, feature = "exon", attr = "gene_id"
    shell:
        r"""
        featureCounts -T {threads} \
          -a "references/decoy/ToxoDB-68_TgondiiRH88.gft" -t {params.feature} -g {params.attr} \
          -s {params.stranded} \
          -o {output.tsv} {input.bam}
        """
rule featurecounts_decoy_multi:
    input:
        bam = "out/filtered/{sample}_decoy_filtered.bam"
    output:
        tsv = "out/counts/{sample}_decoy_multi_featureCounts.txt",
        sum = "out/counts/{sample}_decoy_multi_featureCounts.txt.summary"
    threads: 8
    params:
        stranded = 1, feature = "exon", attr = "gene_id"
    shell:
        r"""
        featureCounts -T {threads} \
          -a "references/decoy/ToxoDB-68_TgondiiRH88.gft" -t {params.feature} -g {params.attr} \
          -s {params.stranded} \
          -M --fraction \
          -o {output.tsv} {input.bam}
        """