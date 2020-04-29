#!/usr/bin/env python

# Author: Maurits Evers (maurits.evers@gmail.com)


import gzip
import numpy as np
from Bio import SeqIO, bgzf
from Bio.SeqRecord import SeqRecord


# Snakemake parameters
infile = snakemake.input[0]
outfile = snakemake.output[0]
n_reads = int(snakemake.params.n_reads)
len_frag = int(snakemake.params.len_frag)


def get_frag(seq, pos, len_frag, i):
    sseq = seq[pos:(pos + len_frag)]
    label = "read_%i_pos%i-%i" % (i, pos, pos + len_frag)
    return (SeqRecord(sseq, id = label, description = ""))


with gzip.open(infile, "rt") as fh:
    record = SeqIO.read(fh, "fasta")
    seq = record.seq
    len = len(seq)
    np.random.seed(2020)
    idx = np.random.choice(len - len_frag, n_reads)
    subseq = [get_frag(seq, pos, len_frag, i)
        for (i, pos) in enumerate(idx, start = 1)]


with bgzf.BgzfWriter(outfile, "wb") as fh:
    SeqIO.write(subseq, fh, "fasta")
