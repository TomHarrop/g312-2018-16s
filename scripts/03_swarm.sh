#!/usr/bin/env bash

singularity exec img/swarm_2.2.2 \
    swarm \
    -t 20 \
    -d 6 \
    -o output/swarm/all_reads.swarm \
    output/swarm/all_reads.fasta
