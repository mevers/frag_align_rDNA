# frag_align_rDNA workflow

This is the repository for the `snakemake` workflow `frag_align_rDNA`. In a
nutshell, the workflow creates artificial reads from the canonical mouse rDNA
sequence [BK000964](https://www.ncbi.nlm.nih.gov/nuccore/BK000964) and aligns
these reads to a hard-masked and unmasked mouse reference genome to help locate
rDNA-like regions in the mouse reference genome.


## Summary

`snakemake` rule DAG

![rulegraph](workflow_rulegraph.png)

The workflow consists of the following steps:

- Download reference genome files from Ensembl: chr1-19, X, Y, MT; do this for
the hard-masked and unmasked sequence files.
- Concatenate individual genome reference files and create a single reference
file for `bowtie2` for both the hard-masked and unmasked reference
- Index the single reference files with `bowtie2-build`
- Download BK000964 rDNA reference sequence from GenBank
- Fragment rDNA sequence by uniformly sampling 1e6 subsequences of length
100 bp, and store resulting sequences in a gzip'ed FASTA file; this is done
using a custom Python script [`workflow/scripts/fragment_seq.py`](workflow/scripts/fragment_seq.py)
- Align rDNA subsequences to the reference genome with `bowtie2`, using default
parameters except for switch `--all` to report *all* alignments.
- Index resulting BAM file with `samtools index`
- Create BigWig coverage track from BAM file using `deeptools` `bamCoverage`


## Details

Mouse rDNA tandem repeats are located on chromosomes (11,) 12, 15, 16, 18 and 19.

References:

- [Gibbons et al., *Concerted copy number variation balances ribosomal DNA dosage in human and mouse genomes*, PNAS 112, 2485 (2015)](https://www.pnas.org/content/112/8/2485)
- [Henderson et al., *The chromosomal location of ribosomal DNA in the mouse*, Chromosoma 49, 155 (1974)](https://link-springer-com.virtual.anu.edu.au/article/10.1007/BF00348887)
- [Xu et al., *Ribosomal DNA copy number loss and sequence variation in cancer*, PLOS Genetics 1006771 (2017)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006771)


## Workflow deployment and requirements

The `snakemake` workflow should be entirely self-contained and reproducible.

Provided `snakemake` is installed, the workflow can be deployed to a new system
by following these steps

```
# Clone repository
git clone https://github.com/mevers/frag_align_rDNA
cd frag_align_rDNA

# Execute workflow
./run_workflow.sh
# Or explicitly:
# snakemake --use-conda
```

External dependencies include `bowtie2`, `samtools` and `deeptools`, and will
be met automatically through `conda` environments. There is no need to manually
install `bowtie2` etc (although it doesn't affect the workflow if these tools
are already installed).


## Contact

Maurits Evers (maurits.evers@gmail.com)

Please raise any questions, concerns, bugs as an Issue on the GitHub project
site.
