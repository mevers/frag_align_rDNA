from os.path import join
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
FTP = FTPRemoteProvider()
NCBI = NCBIRemoteProvider(email = "someone@example.com")
HTTP = HTTPRemoteProvider()


# This is very VERY slow, see:
# https://bitbucket.org/snakemake/snakemake/issues/1275/ftp-remote-inputs-seem-to-slow-down-dag
# It would be faster to use wget or curl, but this adds a dependency
rule download_genome_ref_units:
    input:
        lambda wildcards:
            FTP.remote(join(
                config["data"][wildcards.organism]["Ensembl"]["base_url"],
                "{units}"))
    output:
        "00_ref_sequences/{organism}/chr_Ensembl/{units}"
    shell:
        """
            mv {input} {output}
        """


# bgzip without intermediate file thanks to:
# https://superuser.com/a/325505
rule download_GenBank_rDNA:
    input:
        lambda wildcards:
            NCBI.remote(
                config["data"][wildcards.organism]["GenBank"]["accession_no"] +
                ".fasta",
                db = "nuccore")
    output:
        "00_ref_sequences/{organism}/rDNA_GenBank/{accession_no}.fa.gz"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
            cat {input} | bgzip > {output}
        """


# Download cytoband data from UCSC
rule download_UCSC_cytoband_data:
    input:
        lambda wildcards:
            HTTP.remote(join(
                config["data"][wildcards.organism]["UCSC"]["base_url"],
                config["data"][wildcards.organism]["UCSC"]["cytoband"]),
            keep_local = True)
    output:
        "05_circos_plots/{organism}/{name}/cytoband.txt.gz"
    shell:
        """
            mv {input} {output}
        """
