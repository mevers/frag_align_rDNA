# Fragment the rDNA sequence and generate reads of length `len_frag`
# Two different methods are implemented:
#   1. `method = "uniform_sample"` to sample uniformly `n_reads`
#   2. `method = "sliding_window"` to use a sliding window approach with step
#      size `step`; `n_reads` is ignored.
# Here we set `len_frag` = `step` to enforce non-overlapping windows. This can
# be changed by setting a different `step`.
rule create_rDNA_frags:
    input:
        lambda wildcards: expand(
            "00_ref_sequences/{{organism}}/rDNA_GenBank/{id}.fa.gz",
            id = config["data"][wildcards.organism]["GenBank"]["accession_no"])
    output:
        "02_rDNA_frags/{organism}/rDNA_frags_len{len_frag}_step{step}.fa.gz"
    conda:
        "../envs/fragment_seq.yaml"
    params:
        method = "sliding_window",
        len_frag = lambda wildcards: wildcards.len_frag,
        step = lambda wildcards: wildcards.step
    script:
        "../scripts/fragment_seq.py"


# Align rDNA fragments to reference genome
# We use default parameters to be consistent with typical workflows
# involving `bowtie2`; explicitly: `-N 0`; `--end-to-end`
# We report all alignments: `--all`
# `bowtie2` results are piped directly into samtools for filtering out unmapped
# reads with `-F 4`; even though it is possible (see reference below) to pipe
# `bowtie2` and/or `samtools view` output directly into `samtools sort`, this
# seems to be very slow, so sorting will be done in a separate rule.
# Ref:
# - https://www.biostars.org/p/374173/,
# - https://biology.stackexchange.com/questions/59493/how-to-convert-bwa-mem-output-to-bam-format-without-saving-sam-file
rule align_rDNA_frags:
    input:
        idx = expand("01_bowtie2_ref/{{organism}}/{{name}}.{idx}.bt2",
            idx = range(1, 5)),
        frags = "02_rDNA_frags/{organism}/"
            "rDNA_frags_len{len_frag}_step{step}.fa.gz"
    output:
        "03_alignment/{organism}/{name}/"
        "unsorted_rDNA_frags_len{len_frag}_step{step}.bam"
    conda:
        "../envs/bowtie2.yaml"
    log:
        "logs/bowtie2_{organism}_{name}_rDNA_frags_len{len_frag}_step{step}.log"
    params:
        index_basename = "01_bowtie2_ref/{organism}/{name}"
    threads: 4
    shell:
        """
            bowtie2 \
            -x {params.index_basename} \
            -U {input.frags} \
            -f \
            -N 0 \
            --end-to-end \
            --all \
            --threads {threads} \
            2> {log} \
            | samtools view -bSh -@ {threads} -F4 - > {output}
        """


# Filter and sort BAM file using `samtools sort`
rule sort_BAM_alignment:
    input:
        "03_alignment/{organism}/{name}/"
        "unsorted_rDNA_frags_len{len_frag}_step{step}.bam"
    output:
        "03_alignment/{organism}/{name}/"
        "rDNA_frags_len{len_frag}_step{step}.bam"
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        """
            samtools sort -@ {threads} -o {output} {input}
        """
