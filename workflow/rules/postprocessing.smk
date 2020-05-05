# Index BAM file
rule index_bam_file:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
            samtools index {input}
        """


# Create BigWig coverage track
# Note that `bamCoverage` needs a BAM file along with an index, so the index
# is a critical input file
rule create_bw_coverage_track:
    input:
        expand("03_alignment/{{organism}}/{{name}}/"
        "rDNA_frags_len{{len_frag}}_step{{step}}.{ext}",
        ext = ["bam", "bam.bai"])
    output:
        "04_output_data/{organism}/{name}/"
        "rDNA_frags_len{len_frag}_step{step}.bw"
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/bamCoverage_{organism}_{name}_len{len_frag}_step{step}.log"
    shell:
        """
            bamCoverage -b {input[0]} -o {output} 2>{log} 1>&2
        """


# Create BED file from BAM alignment
# This uses bedtools' bamtobed
rule create_bed_from_bam:
    input:
        "03_alignment/{organism}/{name}/"
        "rDNA_frags_len{len_frag}_step{step}.bam"
    output:
        "04_output_data/{organism}/{name}/"
        "rDNA_frags_len{len_frag}_step{step}.bed"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
            bedtools bamtobed -i {input} -ed -bed12 > {output}
        """
