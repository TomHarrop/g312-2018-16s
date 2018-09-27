#!/usr/bin/env bash

cat data/reads/with_spacer/*fa |
singularity exec img/bbmap_38.00 \
    filterbyname.sh \
    in=stdin.fasta \
    names=output/annotated/aligned_otus.txt \
    out=output/annotated/aligned_otus.fasta \
    include=t

cp output/annotated/aligned_otus.fasta \
    output/taxonomy/aligned_otus.fasta

(
cd output/taxonomy
singularity exec \
    ../../img/mothur_1.40.5 \
    mothur "#classify.seqs(fasta=aligned_otus.fasta, template=../../data/silva/silva.seed_v132.align, taxonomy=../../data/silva/silva.seed_v132.tax,  processors=20)"
 )
