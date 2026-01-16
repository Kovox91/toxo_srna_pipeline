# Define the deafult conda environment
conda: "environment.yaml"   

# Import necessary modules
import pandas as pd
import os, subprocess, socket, getpass, time

# Load configuration
configfile: "config.yaml"

# Set up Notifications
HOST  = socket.gethostname()
USER  = getpass.getuser()
START = time.time()

# Configure via env vars.
NTFY_TOPIC = os.environ.get("NTFY_TOPIC", "")
NTFY_URL   = os.environ.get("NTFY_URL", "https://ntfy.sh")  # change if you self-host

def _notify(title, msg, tags="snake,white_check_mark"):
    if not NTFY_TOPIC:
        return
    try:
        subprocess.run(
            [
                "curl","-sS",
                "-H", f"Title: {title}",
                "-H", "Priority: high",
                "-H", f"Tags: {tags}",
                "-d", msg,
                f"{NTFY_URL.rstrip('/')}/{NTFY_TOPIC}",
            ],
            check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
    except Exception:
        pass  # never block the workflow on notify

def _fmt_duration(sec):
    h = int(sec // 3600); m = int((sec % 3600) // 60); s = int(sec % 60)
    return f"{h:02d}:{m:02d}:{s:02d}"

onstart:
    _notify(
        "Snakemake ▶️ started",
        f"{USER}@{HOST}\nDir: {os.getcwd()}",
        tags="snake,play_button"
    )

onsuccess:
    _notify(
        "Snakemake ✅ done",
        f"{USER}@{HOST}\nDir: {os.getcwd()}\nDuration: {_fmt_duration(time.time()-START)}",
        tags="snake,white_check_mark"
    )

onerror:
    _notify(
        "Snakemake ❌ failed",
        f"{USER}@{HOST}\nDir: {os.getcwd()}\nAfter: {_fmt_duration(time.time()-START)}\nCheck: .snakemake/log/ and run.log",
        tags="snake,cross_mark"
    )
# ---- end notifications ----


# Load sample information from CSV file
samples_df = pd.read_csv(config['samples_csv'])

# Create a list of sample names
sample_names = list(samples_df['sample'])

# Define rule all to specify final output files
rule all:
    input:
        # === Sorting and Indexing ===
        mito_bam = expand("out/filtered/{sample}_mito_filtered.bam", sample=sample_names),
        mito_bai = expand("out/filtered/{sample}_mito_filtered.bam.bai", sample=sample_names),

        decoy_bam = expand("out/filtered/{sample}_decoy_filtered.bam", sample=sample_names),
        decoy_bai = expand("out/filtered/{sample}_decoy_filtered.bam.bai", sample=sample_names),

        # === Summarising .bam Files strand specific ===
        mito_bed_plus = expand("out/filtered/strand_cov/{sample}_mito_filtered_plus.bedgraph", sample=sample_names),
        mito_bed_minus = expand("out/filtered/strand_cov/{sample}_mito_filtered_minus.bedgraph", sample=sample_names),

        # === Start and End Positions ===
        mito_read_start_end = expand("out/filtered/strand_cov/{sample}_mito_filtered_read_start_end.bed", sample=sample_names),

        # === Counting ===
        counts_mito_short_sense = expand("out/counts/mito/shortFeatures/sense/{sample}_mito_featureCounts.txt", sample=sample_names),
        counts_mito_short_antisense = expand("out/counts/mito/shortFeatures/antisense/{sample}_mito_featureCounts.txt", sample=sample_names),
        counts_mito_long = expand("out/counts/mito/longFeatures/{sample}_mito_featureCounts.txt", sample=sample_names),

        counts_decoy = expand("out/counts/decoy/{sample}_decoy_featureCounts.txt", sample=sample_names),

        # === PolyA Analysis ===
        polyAsummary = expand("out/polyA/{sample}_polyA_lengths.tsv", sample=sample_names)


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
        bowtie -q -v 1 -k 4 --best --strata -p {threads} \
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
        mito  = "out/mapped/{sample}_mito_mapped_sorted.bam",   
        decoy = "out/mapped/{sample}_decoy_mapped_sorted.bam"
    output:
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

# Summarise Filtered mito .bam files by strand
rule strand_cov:
    input:
        bam = "out/filtered/{sample}_mito_filtered.bam"
    output:
        plus  = "out/filtered/strand_cov/{sample}_mito_filtered_plus.bedgraph",
        minus = "out/filtered/strand_cov/{sample}_mito_filtered_minus.bedgraph"
    shell:
        """
        bedtools genomecov -ibam {input.bam} -strand + -bg > {output.plus}
        bedtools genomecov -ibam {input.bam} -strand - -bg > {output.minus}
        """

rule start_end_positions:
    input:
        bam = "out/filtered/{sample}_mito_filtered.bam"
    output:
        reads_start_end = "out/filtered/strand_cov/{sample}_mito_filtered_read_start_end.bed"
    shell:
        r"""
        bedtools bamtobed -split -i {input.bam} \
            | awk 'BEGIN{OFS="\t"}
            {
            if ($6 == "+") {
                s = $2 + 1; e = $3
            } else {
                s = $3;     e = $2 + 1
            }
            start[$1 FS s FS $6]++
            end[$1 FS e FS $6]++
            }
            END {
            for (k in start) {
                split(k, a, FS)
                print a[1], a[2], start[k], a[3], "start"
            }
            for (k in end) {
                split(k, a, FS)
                print a[1], a[2], end[k], a[3], "end"
            }
            }' > {output.reads_start_end}
        """

rule featurecounts_mito_short_antisense:
    input:
        bam = "out/filtered/{sample}_mito_filtered.bam"
    output:
        tsv = "out/counts/mito/shortFeatures/antisense/{sample}_mito_featureCounts.txt",
        sum = "out/counts/mito/shortFeatures/antisense/{sample}_mito_featureCounts.txt.summary"
        rep = "out/counts/mito/shortFeatures/sense/{sample}_mito_filtered.bam.featureCounts"
    threads: 8
    params:
        stranded = 2, feature = "rRNA", attr = "rRNA_id"
    shell:
        r"""
        featureCounts -T {threads} \
          -a "references/pseudo_genome/RNA_index_rebuild.gtf" -t {params.feature} -g {params.attr} \
          -s {params.stranded} \
          --fracOverlapFeature 0.75 \
          --largestOverlap \
          -o {output.tsv} {input.bam}
        """

rule featurecounts_mito_long:
    input:
        bam = "out/filtered/{sample}_mito_filtered.bam"
    output:
        tsv = "out/counts/mito/longFeatures/{sample}_mito_featureCounts.txt",
        sum = "out/counts/mito/longFeatures/{sample}_mito_featureCounts.txt.summary"
    threads: 8
    params:
        stranded = 1, feature = "gene", attr = "gene_id"
    shell:
        r"""
        featureCounts -T {threads} \
          -a "references/pseudo_genome/RNA_index_long_only.gtf" -t {params.feature} -g {params.attr} \
          -s {params.stranded} \
          --fracOverlap 0.8 \
          -o {output.tsv} {input.bam}
        """

rule featurecounts_decoy:
    input:
        bam = "out/filtered/{sample}_decoy_filtered.bam"
    output:
        tsv = "out/counts/decoy/{sample}_decoy_featureCounts.txt",
        sum = "out/counts/decoy/{sample}_decoy_featureCounts.txt.summary"
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

rule featurecounts_mito_short_sense:
    input:
        bam = "out/filtered/{sample}_mito_filtered.bam"
    output:
        tsv = "out/counts/mito/shortFeatures/sense/{sample}_mito_featureCounts.txt",
        sum = "out/counts/mito/shortFeatures/sense/{sample}_mito_featureCounts.txt.summary"
    threads: 8
    params:
        stranded = 1, feature = "rRNA", attr = "rRNA_id"
    shell:
        r"""
        featureCounts -T {threads} \
          -a "references/pseudo_genome/RNA_index_rebuild.gtf" -t {params.feature} -g {params.attr} \
          -s {params.stranded} \
          --fracOverlapFeature 0.8 \
          --largestOverlap \
          -R CORE \
          -o {output.tsv} {input.bam}
        """

rule filter_assigned_ids:
    input:
        report="out/counts/mito/shortFeatures/sense/{sample}_mito_filtered.bam.featureCounts"
    output:
        filtered_ids=temp("out/polyA/{sample}.filtered_ids.txt")
    shell:
        r"""
        grep -v '^#' {input.report} \
        | awk '$2 == "Assigned"' \
        > {output.filtered_ids}

        rm {input.report}
        """

rule trim_adapters_only:
    input:
        R1="in/{sample}_R1_001.fastq.gz"
    output:
        R1_trimmed=temp("out/polyA/{sample}_R1_001_trimmed.fq.gz")
    threads: 8
    shell:
        r"""
        cutadapt -j {threads} -q 20 -m 15 --trim-n \
          -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
          -o {output.R1_trimmed} {input.R1}
        """

rule extract_assigned_reads:
    input:
        R1_trimmed="out/polyA/{sample}_R1_001_trimmed.fq.gz",
        filtered_ids="out/polyA/{sample}.filtered_ids.txt"
    output:
        assigned_fastq=temp("out/polyA/{sample}_assigned.fastq")
    shell:
        r"""
        # extract only the read IDs (1st column)
        awk '{{print $1}}' {input.filtered_ids} > out/polyA/{wildcards.sample}.idlist.txt

        # subset FASTQ to those IDs
        seqtk subseq {input.R1_trimmed} out/polyA/{wildcards.sample}.idlist.txt > {output.assigned_fastq}

        rm out/polyA/{wildcards.sample}.idlist.txt
        """

rule polyA_lengths:
    input:
        assigned_fastq = "out/polyA/{sample}_assigned.fastq",
        id_feature     = "out/polyA/{sample}.filtered_ids.txt"
    output:
        stats = temp("out/polyA/{sample}_polyA_lengths_with_features.tsv")
    shell:
        r"""
        awk 'NR==FNR {{
                 # First pass: read FASTQ and compute per-read length & polyA tail
                 if (FNR%4==1) {{
                     id=$1
                     sub(/^@/, "", id)
                     cur=id
                 }} else if (FNR%4==2) {{
                     seq=$1
                     L=length(seq)
                     tail=0
                     for (i=L; i>0; i--) {{
                         c=substr(seq, i, 1)
                         if (c=="A" || c=="a") tail++
                         else break
                     }}
                     len[cur]=L
                     poly[cur]=tail
                 }}
                 next
             }}
             # Second pass: go through filtered_ids and append L and tail
             {{
                 id=$1
                 L = (id in len ? len[id] : "NA")
                 tail = (id in poly ? poly[id] : "NA")
                 print $0 "\t" L "\t" tail
             }}' {input.assigned_fastq} {input.id_feature} > {output.stats}
        """

rule summarise_polyA:
    input:
        stats = "out/polyA/{sample}_polyA_lengths_with_features.tsv"
    output:
        summary = "out/polyA/{sample}_polyA_lengths.tsv"
    shell:
        r"""
        cut -f4,6 {input.stats} \
        | sort \
        | uniq -c \
        | awk '{{print $2"\t"$3"\t"$1}}' \
        > {output.summary}
        """