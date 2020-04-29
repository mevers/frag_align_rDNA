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
rule create_bw_coverage_track:
    input:
        "03_alignment/{id}/rDNA_frags_{name}.bam"
    output:
        "04_coverage_tracks/{id}/rDNA_frags_{name}.bw"
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/bamCoverage_{id}_{name}.log"
    shell:
        """
            bamCoverage -b {input} -o {output} 2> {log} 1>&2
        """
