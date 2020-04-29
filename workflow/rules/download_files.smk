from os.path import join
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
FTP = FTPRemoteProvider()
NCBI = NCBIRemoteProvider(email = "someone@example.com")


# This is very VERY slow, see:
# https://bitbucket.org/snakemake/snakemake/issues/1275/ftp-remote-inputs-seem-to-slow-down-dag
# It would be faster to use wget or curl, but this adds a dependency
rule download_genome_ref_units:
    input:
        lambda wildcards:
            FTP.remote(join(
                config["genome_reference"][wildcards.id]["base_url"],
                "{units}"))
    output:
        "00_ref_sequences/chr_Ensembl/{id}/{units}"
    shell:
        """
            mv {input} {output}
        """


# gzip without intermediate file thanks to:
# https://superuser.com/a/325505
rule download_GenBank_rDNA:
    input:
        NCBI.remote("BK000964.3.fasta", db = "nuccore")
    output:
        "00_ref_sequences/rDNA_GenBank/BK000964.3.fa.gz"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
            cat {input} | bgzip > {output}
        """
