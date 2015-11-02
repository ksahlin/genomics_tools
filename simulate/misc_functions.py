import random
##
# General functions

def reverse_complement(string):
    """
        Reverse complements a DNA string
        Arguments:
        string  - A DNA string

        Returns:
            A python string object that represents the
            reverse complement of the input (DNA) string.

    """
    rev_nuc={'A':'T','C':'G','G':'C','T':'A','N':'N','X':'X', 
    'a':'t', 't':'a','g':'c','c':'g'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def read_in_fasta_file(fasta_file):
    """
        Reads a fasta file into memory.

        Arguments:
        fasta_file - A python file object. The file should be in 
        fasta format.

        Returns:
            A python dictionary with accessions as keys and
            sequences as values.

    """  
    seqs = {}
    k = 0
    temp = []
    accession = ''
    for line in fasta_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            seqs[accession] = ''
            k += 1
        elif line[0] == '>':
            temp = ''.join(temp)
            seqs[accession] = temp
            temp = []
            accession = line[1:].strip().split()[0]
        else:
            temp.append(line.strip())
    
    
    temp = ''.join(temp)
    seqs[accession] = temp
    return(temp)


def insertion():
    return(''.join( random.choice('AGCT') for i in range( random.randint(1,10) ) ))

def deletion():
    return( abs(int(random.gauss(4,2))))

def mutation():
    return(random.choice('AGCT'))