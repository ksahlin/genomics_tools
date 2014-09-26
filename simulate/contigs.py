import random
from genomics_tools.simulate.misc_functions import reverse_complement

def generate_contigs(genome,min_c_len,max_c_len,min_distance,max_distance):
    position = 0
    index = 0
    while True:
        contig_len = random.randrange(min_c_len,max_c_len)
        position += random.randrange(min_distance,max_distance)
        if position + contig_len > len(genome):
            break

        rev_comp = random.randrange(0,2)
        if rev_comp:
            yield '>c{0},pos:{1}-{2},rc:1\n{3}\n'.format(index,
                position,position+contig_len, reverse_complement( genome[position:position+contig_len]))
        else:
            yield '>c{0},pos:{1}-{2},rc:0\n{3}\n'.format(index,
                position,position+contig_len, genome[position:position+contig_len])
        index += 1
        position += contig_len

def generate_contigs_two_sizes(genome,small_size,large_size,min_distance,max_distance,distr_weight_large_ctgs):
    position = 0
    index = 0
    while True:
        r = random.uniform(0,1)
        if r < distr_weight_large_ctgs:
            contig_len = large_size
        else:
            contig_len = small_size

        position += random.randrange(min_distance,max_distance)
        if position + contig_len > len(genome):
            break

        rev_comp = random.randrange(0,2)
        if rev_comp:
            yield '>c{0},pos:{1}-{2},rc:1\n{3}\n'.format(index,
                position,position+contig_len, reverse_complement( genome[position:position+contig_len]))
        else:
            yield '>c{0},pos:{1}-{2},rc:0\n{3}\n'.format(index,
                position,position+contig_len, genome[position:position+contig_len])
        index += 1
        position += contig_len

