#!/usr/bin/env bash

cat data/reads/with_spacer/*fa |
singularity exec img/bbmap_38.00 \
    filterbyname.sh \
    in=stdin.fasta \
    names=output/filter/kept_otus.txt \
    out=output/filter/kept_otus.fasta \
    include=t

cp output/filter/kept_otus.fasta \
    output/annotated/kept_otus.fasta

(
cd output/annotated
singularity exec \
    ../../img/mothur_1.40.5 \
    mothur "#align.seqs(candidate=kept_otus.fasta, template=../../data/silva/silva.seed_v132.align, processors=20, flip=f)"
 )
