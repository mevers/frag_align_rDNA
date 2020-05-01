#!/usr/bin/env python

# Author: Maurits Evers (maurits.evers@gmail.com)


import gzip
import numpy as np
from Bio import SeqIO, bgzf
from Bio.SeqRecord import SeqRecord


# Snakemake parameters
infile = snakemake.input[0]
outfile = snakemake.output[0]
method = snakemake.params.method
if method == "uniform_sample":
    n_reads = int(snakemake.params.n_reads)
len_frag = int(snakemake.params.len_frag)
step = int(snakemake.params.step)


# Create FASTA record from substring
def get_frag(seq, pos, len_frag, i):
    sseq = seq[pos:(pos + len_frag)]
    label = "read_%i_pos%i-%i" % (i, pos, pos + len_frag)
    return (SeqRecord(sseq, id = label, description = ""))


# Main
with gzip.open(infile, "rt") as fh:
    record = SeqIO.read(fh, "fasta")
    seq = record.seq
    len = len(seq)

    if method == "uniform_sample":
        #np.random.seed(2020)
        idx = np.random.choice(len - len_frag, n_reads)
    elif method == "sliding_window":
        n_reads = int(((len - len_frag) / step) + 1)
        idx = np.arange(0, n_reads * step, step)

    subseq = [get_frag(seq, pos, len_frag, i)
        for (i, pos) in enumerate(idx, start = 1)]


# Write to bgzip'ed FASTA
with bgzf.BgzfWriter(outfile, "wb") as fh:
    SeqIO.write(subseq, fh, "fasta")
