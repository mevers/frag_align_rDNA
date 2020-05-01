rule concatenate_genome_ref:
    input:
        genome = lambda wildcards: expand(
            "00_ref_sequences/{{organism}}/chr_Ensembl/{units}",
            units = config["data"][wildcards.organism]["Ensembl"]["units"])
    output:
        "01_bowtie2_ref/{organism}/{name}.fa.gz"
    shell:
        """
            cat {input.genome} > {output}
        """


# This creates a small index. Note that if the reference file contains more
# than 4 billion nucleotides, `bowtie2-build` will automatically create a large
# index which will result in different output files, causing this rule to fail.
# This can be fixed by either including a `make_genome_ref_large_index` rule,
# or by forcing `bowtie2-build` to always generate a large index with
# `--large-index`. The current workflow uses reference sequences < 4 billion
# nucleotides, so this is not an issue.
rule make_genome_ref_index:
    input:
        "01_bowtie2_ref/{organism}/{name}.fa.gz"
    output:
        expand("01_bowtie2_ref/{{organism}}/{{name}}.{idx}.bt2",
            idx = range(1, 5)),
        expand("01_bowtie2_ref/{{organism}}/{{name}}.rev.{idx}.bt2",
            idx = range(1, 3))
    conda:
        "../envs/bowtie2.yaml"
    log:
        "logs/bowtie2-build_{organism}_{name}.log"
    params:
        index_basename = "01_bowtie2_ref/{organism}/{name}"
    threads: 4
    shell:
        """
            bowtie2-build --threads {threads} \
            {input} {params.index_basename} > {log}
        """
