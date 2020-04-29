rule concatenate_genome_ref:
    input:
        genome = lambda wildcards: expand(
            "00_ref_sequences/chr_Ensembl/{{id}}/{units}",
            units = config["genome_reference"][wildcards.id]["units"])
    output:
        "01_bowtie2_ref/{id}/{name}.fa.gz"
    shell:
        """
            cat {input.genome} > {output}
        """


rule make_genome_ref_index:
    input:
        "01_bowtie2_ref/{id}/{name}.fa.gz"
    output:
        expand("01_bowtie2_ref/{{id}}/{{name}}.{idx}.bt2",
            idx = range(1, 5)),
        expand("01_bowtie2_ref/{{id}}/{{name}}.rev.{idx}.bt2",
            idx = range(1, 3))
    conda:
        "../envs/bowtie2.yaml"
    log:
        "logs/bowtie2-build_{id}_{name}.log"
    params:
        index_basename = "01_bowtie2_ref/{id}/{name}"
    shell:
        """
            bowtie2-build {input} {params.index_basename} > {log}
        """
