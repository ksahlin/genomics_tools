
def to_AGP(AGP_file,acc,positions,seqs, gaps, directions):
	"""
		Takes a write to file object and four indexable structures (e.g. lists, tuples).
		positions, seqs and directions are of length n and gaps are of length n-1.

		AGP_file		- A file object with write acess
		acc				- A string thet is the accession of the contig/scaffold/exon/or other structure
		positions		- indexable datastructure of tuples (start_pos, end_pos) for each of the sequences in the structure
		sequences 		- indexable datastructure of strings that are the sequences
		gaps 			- indexable datastructure of gap lengths between two strings
		directions		- indexable datastructure of booleans denoting the direction of sequences, '+' for forward and '-' for reverse.

	"""
	assert len(positions) == len(seqs) == len(directions) == len(gaps)+1
	component_count = 0
	for i in range(len(seqs)-1):
	    sign = '+' if directions[i] else '-'           
	    if i > 0 and gaps[i-1] > 0:
	        component_count += 1
	        print >> AGP_file, acc + '\t' + str(positions[i-1][1] + 1) + '\t' + str(positions[i][0]-1) + '\t' + str(component_count) + '\t' + 'N\t' + str(gaps[i]) + '\tfragment\tyes\t'
	    component_count += 1
	    print >> AGP_file, acc + '\t' + str(positions[i][0]) + '\t' + str(positions[i][1]) + '\t' + str(component_count) + '\t' + 'W\t' + 'ctg{0}'+str(i) + '\t1\t' + str(positions[i][1] - positions[i][0] + 1) + '\t' + sign

def to_GFF(GFF_file, seq_name,source,feature,start_coord,end_coord, score, strand, frame, attribute):
	info_tuple =(seq_name, source, feature, start_coord, end_coord, score, strand, frame, attribute)
	GFF_file.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(*info_tuple))
