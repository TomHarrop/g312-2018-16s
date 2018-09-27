#!/usr/bin/env python3

from Bio import SeqIO
import csv
import re

mothur_names = 'output/reads/all_reads.names'
swarm_results = 'output/swarm/all_reads.swarm'
dereplicated_fasta = 'output/reads/all_reads.unique.fasta'
output_names = 'output/swarm/cluster.names'
output_fasta = 'output/swarm/cluster.fasta'
long_table_file = 'output/swarm/long_table.tab'

# generate a dict of names from the mothur namefile
with open(mothur_names, 'rt') as mothur:
    split_strings = [str(x).rstrip('\n').split('\t') for x in mothur]
    mothur_name_dict = {x[0]: x[1] for x in split_strings}

# read the swarm results
with open(swarm_results, 'rt') as swarm:
    readlists = [str(x).rstrip('\n').split(' ') for x in swarm]

# discard singletons
kept_swarms = {x[0]: x for x in readlists
               if len(x) > 1}

# remove the trailing number of reads from the header and the read list
renamed_swarms = {}
for key in kept_swarms:
    new_key = re.sub('_\d+$', '', key)
    renamed_swarms[new_key] = [re.sub('_\d+$', '', x)
                               for x in kept_swarms[key]]

# retrieve the full list of reads from the mothur names dict
names_output = []
for key in renamed_swarms:
    names_list = ''
    for read_name in renamed_swarms[key]:
        if mothur_name_dict[read_name]:
            names_list += mothur_name_dict[read_name] + ','
    names_output.append([key, names_list.rstrip(',')])

# read over the renamed_swarms a second time to generate a long-form data.frame
# for loading into R
with open(long_table_file, 'wt') as long_table:
    long_table.write('otu_id\tread_name\n')
    for key in renamed_swarms:
        read_list = []
        for read_name in renamed_swarms[key]:
            if mothur_name_dict[read_name]:
                mothur_reads = mothur_name_dict[read_name].split(',')
                for mothur_read in mothur_reads:
                    read_list.append(mothur_read)
        for read in read_list:
            long_table.write('{0}\t{1}\n'.format(key, read))


# write as names file
with open(output_names, 'w', newline='') as mothur:
    csvwriter = csv.writer(mothur, delimiter=' ')
    csvwriter.writerows(names_output)

# filter fasta by mothur_output
fasta_input = SeqIO.parse(dereplicated_fasta, 'fasta')
records = [x for x in fasta_input]
output_records = [x for x in records if x.id in renamed_swarms.keys()]
SeqIO.write(output_records, output_fasta, 'fasta')
