#!/bin/bash

bedtools intersect \
	-a ../04_output_data/human/GRCh38/rDNA_frags_len100_step100.bed \
	-b ./human/GRCh38/RepeatMasker_GRCh38.bed.gz \
	-wo \
	-f 0.5 > human/GRCh38/overlap_rDNA_frags_len100_step100.tsv

bedtools intersect \
	-a ../04_output_data/human/GRCh38/rDNA_frags_len500_step500.bed \
	-b ./human/GRCh38/RepeatMasker_GRCh38.bed.gz \
	-wo \
	-f 0.5 > human/GRCh38/overlap_rDNA_frags_len500_step500.tsv

bedtools intersect \
	-a ../04_output_data/mouse/GRCm38/rDNA_frags_len100_step100.bed \
	-b ./mouse/GRCm38/RepeatMasker_GRCm38.bed.gz \
	-wo \
	-f 0.5 > mouse/GRCm38/overlap_rDNA_frags_len100_step100.tsv

bedtools intersect \
	-a ../04_output_data/mouse/GRCm38/rDNA_frags_len500_step500.bed \
	-b ./mouse/GRCm38/RepeatMasker_GRCm38.bed.gz \
	-wo \
	-f 0.5 > mouse/GRCm38/overlap_rDNA_frags_len500_step500.tsv
