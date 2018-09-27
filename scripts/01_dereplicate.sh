#!/usr/bin/env bash

cat data/reads/no_spacer/*.fa \
    > output/reads/all_reads.fasta

(
cd output/reads
singularity exec \
    ../../img/mothur_1.40.5 \
    mothur "#unique.seqs(fasta=all_reads.fasta)"
 )

