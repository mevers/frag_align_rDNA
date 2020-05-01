rule generate_circos_plot:
    input:
        bam = "03_alignment/{organism}/{name}/"
        "rDNA_frags_len{len_frag}_step{step}.bam",
        bai = "03_alignment/{organism}/{name}/"
        "rDNA_frags_len{len_frag}_step{step}.bam.bai",
        rDNA = lambda wildcards: expand(
            "00_ref_sequences/{{organism}}/rDNA_GenBank/{acc}.fa.gz",
            acc = config["data"][wildcards.organism]["GenBank"]["accession_no"]),
        cyto = "05_circos_plots/{organism}/{name}/cytoband.txt.gz"
    output:
        pdf = "05_circos_plots/{organism}/{name}/"
        "circos_rDNA_frags_len{len_frag}_step{step}.pdf",
        png = "05_circos_plots/{organism}/{name}/"
        "circos_rDNA_frags_len{len_frag}_step{step}.png"
    conda:
        "../envs/generate_circos_plot.yaml"
    script:
        "../scripts/generate_circos_plot.R"
