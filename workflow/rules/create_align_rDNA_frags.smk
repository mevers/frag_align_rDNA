# Fragment the rDNA sequence and generate `n_reads` of length `len_frag`
rule create_rDNA_frags:
    input:
        "00_ref_sequences/rDNA_GenBank/BK000964.3.fa.gz"
    output:
        "02_rDNA_frags/rDNA_frags.fa.gz"
    params:
        n_reads = 1e6,
        len_frag = 100
    script:
        "../scripts/fragment_seq.py"


# Align rDNA fragments to reference genome
# We use default parameters to be consistent with typical workflows
# involving bowtie2; explicitly: `-N 0`; `--end-to-end`
# We report all alignments: `--all`
rule align_rDNA_frags:
    input:
        idx = expand("01_bowtie2_ref/{{id}}/{{name}}.{idx}.bt2",
            idx = range(1, 5)),
        frags = "02_rDNA_frags/rDNA_frags.fa.gz"
    output:
        "03_alignment/{id}/rDNA_frags_{name}.bam"
    conda:
        "../envs/bowtie2.yaml"
    log:
        "logs/bowtie2_{id}_{name}.log"
    params:
        index_basename = "01_bowtie2_ref/{id}/{name}"
    shell:
        """
            bowtie2 \
            -x {params.index_basename} \
            -U {input.frags} \
            -f \
            -N 0 \
            --end-to-end \
            --all \
            2> {log} \
            | samtools view -bSu -F4 - | samtools sort - > {output}
        """
