# -*- coding: utf-8 -*-



import random
#import numpy
#import sys
#import check_GC_content
import argparse
from collections import Counter

# Markov Chain: Transition matrix
transition = {
    'Sub' : {'Sub':0.99, 'TDup':0.001, 'IDup':0.001, 'Ins':0.004,'Inv':0.004},
    'TDup' : {'Sub':0.99, 'TDup':0.006, 'IDup':0.001, 'Ins':0.001,'Inv':0.001},
    'IDup' : {'Sub':0.99, 'TDup':0.001, 'IDup':0.003, 'Ins':0.005,'Inv':0.001},
    'Ins' : {'Sub':0.99, 'TDup':0.002, 'IDup':0.002, 'Ins':0.003,'Inv':0.003},
    'Inv' : {'Sub':0.99, 'TDup':0.001, 'IDup':0.002, 'Ins':0.003,'Inv':0.004}}

# transition = {
#     'Sub' : {'Sub':0.95, 'TDup':0.00, 'IDup':0.03, 'Ins':0.01,'Inv':0.01},
#     'TDup' : {'Sub':0.95, 'TDup':0.00, 'IDup':0.03, 'Ins':0.01,'Inv':0.01},
#     'IDup' : {'Sub':0.95, 'TDup':0.00, 'IDup':0.03, 'Ins':0.01,'Inv':0.01},
#     'Ins' : {'Sub':0.95, 'TDup':0.00, 'IDup':0.03, 'Ins':0.01,'Inv':0.01},
#     'Inv' : {'Sub':0.95, 'TDup':0.00, 'IDup':0.03, 'Ins':0.01,'Inv':0.01}}

transition = {
    'Sub' : {'Sub':0.995, 'TDup':0.00, 'IDup':0.003, 'Ins':0.01,'Inv':0.001},
    'TDup' : {'Sub':0.995, 'TDup':0.00, 'IDup':0.003, 'Ins':0.01,'Inv':0.001},
    'IDup' : {'Sub':0.995, 'TDup':0.00, 'IDup':0.003, 'Ins':0.01,'Inv':0.001},
    'Ins' : {'Sub':0.995, 'TDup':0.00, 'IDup':0.003, 'Ins':0.01,'Inv':0.001},
    'Inv' : {'Sub':0.995, 'TDup':0.00, 'IDup':0.003, 'Ins':0.01,'Inv':0.001}}



def Sub(genome):

# Substitutes one nt to another
# Arguments: Genome

    #print('sub')
    gen = len(genome)
    number = random.randrange(0,gen)
    substitutes = ['A','C','G','T']
    replace = random.choice(substitutes)
    genome[number] = replace
    return(genome)
    





def TDup(genome):

# Makes tandem duplications; repeats of a sequence next to it
# Arguments: Genome

    #print('T')
    gen = len(genome)
    copy = []
    copy_length = abs(int(random.normalvariate(2, gen/20)))
    position = random.randrange(0,gen-copy_length)
    #u = position

    stop = position + copy_length
    copy = [ genome[u] for u in range(position,stop) ]
    # for u in range (position,stop):
    #     copy.append(genome[u])
    #     u += 1
    genome[stop:stop] = copy
    return(genome)






def IDup(genome):
    
# Makes an interspread duplication; repeats of a sequence on another place in the genome
# Arguments: Genome

    #print('I')
    gen = len(genome)
    copy = []
    copy_length = abs(int(random.normalvariate(2, gen/20)))
    position = random.randrange(0,gen-copy_length)
    #u = position
    stop = position + copy_length
    copy = [ genome[u] for u in range(position,stop) ]
    # for u in range (position,stop):
    #     copy.append(genome[u])
    #     u += 1
    inter = random.randrange(0,gen)
    genome[inter:inter] = copy
    return(genome)






def Ins(genome):

# Inserts a sequence
# Arguments: Genome

    #print('Ins')
    gen = len(genome)
    #copy = []
    copy_length = abs(int(random.expovariate(1/10.0)))
    #gen2 = len(copy)
    nucleobases = ['A','C','G','T']
    copy = [ random.choice(nucleobases) for i in range(copy_length) ]
    # while gen2 <= copy_length-1:
    #     copy.append(random.choice(nucleobases))
    #     gen2 = len(copy)
    inter = random.randrange(0,gen)
    genome[inter:inter] = copy
    return(genome)
    





def Inv(genome):

# Inserts a inverted sequence in the genome
# Arguments: Genome

    #print('Inv')
    gen = len(genome)
    copy = []
    copy_length = abs(int(random.normalvariate(2, gen/20)))
    position = random.randrange(0,gen-copy_length)
    #u = position
    stop = position + copy_length
    # print genome[position:stop]
    # for u in range (position,stop):
    #     copy.insert(0,genome[u])
    #     u += 1
    # print copy
    copy = [genome[nucl] for nucl in range(stop -1, position -1, -1)]

    #print inv
    #assert inv == copy
    inter = random.randrange(0,gen)
    #print genome
    genome[inter:inter] = copy
    return genome




def probability_func(state):

# Markov chain - change current state depending on what the last current state was. Uses the transition matrix
# Argument: dictionary with states and probabilities

    count = 0
    choice = random.random()
    for state_key, prob in state.iteritems():
        count += prob
        if choice <= count:
            current_state = state_key
            return current_state



def new_state(genome, finish_length):

# Sends the genome to the different functions above depending on the current state
# Arguments: The genome
    
    current_state = 'Sub'
    genome_length = len(genome)
    genome = list(genome)
    counter = 0
    while genome_length < finish_length:
        counter += 1
        #if counter % 100 == 0:
        #    print genome_length
        current_state = probability_func(transition[current_state])
        if current_state == 'Sub':
            genome = Sub(genome)
            genome_length = len(genome)
        elif current_state == 'TDup':
            genome = TDup(genome)
            genome_length = len(genome)
        elif current_state == 'IDup':
            genome = IDup(genome)
            genome_length = len(genome)
        elif current_state == 'Ins':
            genome = Ins(genome)
            genome_length = len(genome)
        elif current_state == 'Inv':
            genome = Inv(genome)
            genome_length = len(genome)
        else:
            print('something is wrong')
    #print('hej3')
    return genome


def count_kmers(genome,k):
    kmer_counter = Counter()
    for i in range(len(genome)-k+1):
        kmer_counter[genome[i:i+k]] += 1
        #if i % 10000 == 0:
        #    print i
    
    return Counter(map(lambda x: x[1],kmer_counter.most_common()))


def generate_genome(genome, length):
    #genome = [random.choice(['A','C','G','T']) for i in range(200)]
    genome = new_state(genome, length)
    genome_str = "".join(genome)
    return genome_str


class Genome(object):
    """docstring for Genome"""
    def __init__(self,probs, length,accession):
        super(Genome, self).__init__()
        self.probs = probs
        self.length = length
        self.accession = accession

        self.genome() # generate a strand

    def generate_bp(self, bp):
        """
        Generates base pairs according to the given
        distribution.
        
        @param probs A list of probabilities for each base pair.
        @param bp A list of base pairs.
        """
        rand = random.random( )
        total = 0.0
        for i in range( len( self.probs ) ):
            total += self.probs[ i ]
            if rand <= total:
                return bp[ i ]

    def genome(self):
        """
        Generates the genome a returns it.
            
        @param probs Base pair probabilities.
        @param length Length of the genome.
        
        @return The sequence of the generated genome.
        
        """
        bp = "ACGT"
        self.sequence = ''.join( self.generate_bp( bp ) for i in range( self.length ) )

    def diploid_copy(self,insertion_rate, deletion_rate, mutation_rate,accession):
        self.copy = ''
        self.diploid_accession = accession
        for i in range(len(self.sequence)):
            if random.uniform(0,1) < mutation_rate:
                self.copy = ''.join([i for i in [self.sequence[0:i],mutation(),self.sequence[i+1:]]])
                continue
            if random.uniform(0,1) < deletion_rate:
                self.copy = ''.join([i for i in [self.sequence[0:i-deletion()],self.sequence[i:]]])
                continue
            if random.uniform(0,1) < insertion_rate:
                self.copy = ''.join([i for i in [self.sequence[0:i],insertion(),self.sequence[i:]]])
                continue

    def genome_fasta_format(self):
        return '>{0}\n{1}'.format(self.accession,self.sequence)

    def diploid_copy_fasta(self):
        return '>{0}\n{1}'.format(self.diploid_accession,self.copy)
        



def main(args):

# Makes a random genome with 500 nt and then increases the length of the genome using the functions above.
# Arguments: length of the random genome, length of the finished genome)

    nucleobases = ['A','C','G','T']
    i = 0
    genome = ''
    while i < args.start_length:
        rand = random.choice(nucleobases)
        genome = genome + rand
        i += 1
    length_genome = len(genome)
    # print genome
    print ('length of genome:', length_genome)
    genome = new_state(genome, args.finish_length)


    genome_str = "".join(genome)
    #print('hej1')
    outfile = open(args.outfile,'w')
    print >>outfile, genome_str
    count_kmers(genome_str,2000)
    # print genome
    #print ('length of genome', len(genome))
    #GC_content = check_GC_content.check_GC_content(genome)
    #print GC_content

    #genome2_str = Genome2.diploid(genome, 0.0001, 0.0001, 0.0001)    
    #print('hej2')
    #outfile2 = open(args.outfile+'copy2.fa','w')
    #print >>outfile2, genome2_str




if __name__ == '__main__':

# Take care of input


    parser = argparse.ArgumentParser(description = "Simulate a genome of desired length")
    parser.add_argument('start_length', type=int, help='The length of the first genome that is random. ')
    parser.add_argument('finish_length', type=int, help='The approximate size of the finished genome. ')
    parser.add_argument('outfile', type=str, help='outfile prefix. ')

    args = parser.parse_args()
    main(args)



