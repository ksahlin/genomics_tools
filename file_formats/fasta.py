
def get_sequences_iter(fasta_file):

    raise NotImplementedError



def fasta_iter(fasta_file):
    """
        Reads a fasta file into memory.

        Arguments:
        fasta_file - A python file object. The file should be in 
        fasta format.

        Returns:
            an iterator over accession, sequence.

    """  

    k = 0
    temp = []
    accession = ''
    for line in fasta_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            k += 1
        elif line[0] == '>':
            temp = ''.join(temp)
            yield accession, temp
            temp = []
            accession = line[1:].strip().split()[0]
        else:
            temp.append(line.strip())
    
    temp = ''.join(temp)
    yield accession, temp