import sys
import os
import random
import operator

from genomics_tools.simulate.misc_functions import reverse_complement, read_in_fasta_file
        


class PairedEndRead(object):
    """docstring for PairedEndRead"""
    def __init__(self,distribution = 'normal',mean=None,sigma=None,read_length= None, min_size=None, max_size=None):
        super(PairedEndRead, self).__init__()
        self.distribution = distribution
        self.mean = mean
        self.sigma = sigma
        self.read_length = read_length
        self.min_size = min_size
        self.max_size = max_size

    def generate(self,reference_accession, reference_sequence, read_index):
        if self.distribution == 'normal':
            self.fragment_length = int(random.gauss(self.mean,self.sigma))
        elif self.distribution == 'uniform':
            self.fragment_length = int(random.uniform(self.min_size,self.max_size))

        if self.fragment_length >= len(reference_sequence): 
            raise Exception("To short reference sequence length for \
                simulated read. \nRead fragment: {0}\nTranscript \
                length:{1}".format(self.fragment_length,len(reference_sequence)))
        
        self.start_pos = random.randrange(len(reference_sequence) - self.fragment_length)
        self.read1 = reference_sequence[self.start_pos : self.start_pos + self.read_length]
        self.read2 = reverse_complement(reference_sequence[self.start_pos + self.fragment_length - self.read_length : self.start_pos+self.fragment_length])
        self.reference_accession = reference_accession
        self.read_index = read_index

    def fastq_format(self):
        r1= '@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\
        \n@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{3}\n'.format(self.read_index,
            self.reference_accession,self.read1,'J'*self.read_length)
        r2= '@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\
        \n@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{3}\n'.format(self.read_index,
            self.reference_accession,self.read2,'J'*self.read_length)        
        yield r1
        yield r2

    def fasta_format(self):
        r1= '>HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\n'.format(self.read_index,
            self.reference_accession,self.read1)
        r2= '>HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{2}\n'.format(self.read_index,
            self.reference_accession,self.read2)        
        yield r1
        yield r2

class MatePairRead(object):
    """docstring for MatePairRead"""
    def __init__(self, distribution = 'normal',mean=None,sigma=None,read_length= None,min_size=None, max_size=None):
        super(MatePairRead, self).__init__()
        self.distribution = distribution
        self.mean = mean
        self.sigma = sigma
        self.read_length = read_length
        self.min_size = min_size
        self.max_size = max_size

    def generate(self,reference_accession, reference_sequence, read_index):
        if self.distribution == 'normal':
            self.fragment_length = int(random.gauss(self.mean,self.sigma))
        elif self.distribution == 'uniform':
            self.fragment_length = int(random.uniform(self.min_size,self.max_size))

        if self.fragment_length >= len(reference_sequence): 
            raise Exception("To short reference sequence length for \
                simulated read. \nRead fragment: {0}\nTranscript \
                length:{1}".format(self.fragment_length,len(reference_sequence)))
        self.start_pos = random.randrange(len(reference_sequence) - self.fragment_length)
        self.read1 = reverse_complement(reference_sequence[self.start_pos : self.start_pos + self.read_length])
        self.read2 = reference_sequence[self.start_pos + self.fragment_length - self.read_length : self.start_pos+self.fragment_length]
        self.reference_accession = reference_accession
        self.read_index = read_index

    def fastq_format(self):
        r1= '@HWUSI-EAS100R:6:73:941:{0},(RevCompOriented)genome:{1}\n{2}\
        \n@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{3}\n'.format(self.read_index,
            self.reference_accession,self.read1,'J'*self.read_length)
        r2= '@HWUSI-EAS100R:6:73:941:{0},(RevCompOriented)genome:{1}\n{2}\
        \n@HWUSI-EAS100R:6:73:941:{0},genome:{1}\n{3}\n'.format(self.read_index,
            self.reference_accession,self.read2,'J'*self.read_length)        
        yield r1
        yield r2

    def fasta_format(self):
        r1= '>HWUSI-EAS100R:6:73:941:{0},(RevCompOriented)genome:{1}\n{2}\n'.format(self.read_index,
            self.reference_accession,self.read1)
        r2= '>HWUSI-EAS100R:6:73:941:{0},(RevCompOriented)genome:{1}\n{2}\n'.format(self.read_index,
            self.reference_accession,self.read2)        
        yield r1
        yield r2




class Transcriptome(object):
    """docstring for Transcriptome"""
    def __init__(self):
        super(Transcriptome, self).__init__()

    def generate_splice_variants(self, nr_transcripts,genome,
        max_intron_size,min_intron_size,max_exon_size, min_exon_size):
        self.transcripts = []
        for i in range(nr_transcripts):
            self.transcripts.append(Transcript(genome,max_intron_size,min_intron_size,max_exon_size, min_exon_size))

        return self.transcripts   

    def fasta_format(self):
        for transcript in self.transcripts:
            yield str(transcript)


class Transcript(object):
    """
    Class for an Transcript object
    
    Attributes:

    genome              - A python string object that is the sequence of
                        the genome.
    max_intron_size     - Maximum size of intron
    min_intron_size     - Minimun size of intron
    max_exon_size       - Maximum size of exon
    min_exon_size       - Minimun size of exon    
    gene_distance       - Distance between genes

    """
    def __init__(self, genome_strand,max_intron_size,min_intron_size,max_exon_size, min_exon_size):
        super(Transcript, self).__init__()
        self.genome_strand = genome_strand
        self.max_intron_size = max_intron_size
        self.min_intron_size = min_intron_size
        self.max_exon_size = max_exon_size
        self.min_exon_size = min_exon_size

        self.get_sequence()
        


    def get_sequence(self):
        """
        Generates an Transcript from a genome
        Returns:
        An exon
        """

        nr_exons = random.randrange(1,5)
        self.intron_length = [random.randrange(self.min_intron_size,self.max_intron_size) for i in range(nr_exons)]
        self.exon_lengths = [random.randrange(self.min_exon_size,self.max_exon_size) for i in range(nr_exons)]
        self.reverse_complement = random.randrange(2)

        if sum(self.intron_length) + sum(self.exon_lengths) >= len(self.genome_strand):
            self.get_sequence()

        if self.reverse_complement:
            self.start_position = random.randrange(sum(self.intron_length) + sum(self.exon_lengths),len(self.genome_strand))
        else:
            self.start_position = random.randrange(0,len(self.genome_strand)-sum(self.intron_length)-sum(self.exon_lengths))


        position = self.start_position
        self.positions = []
        self.sequence = ''
        for e_len,i_len in zip(self.exon_lengths,self.intron_length):
            if self.reverse_complement:
                self.sequence += reverse_complement( self.genome_strand[position-e_len:position])
                self.positions.append((position-e_len,position))
                position -= e_len + i_len

            else:

                self.sequence += self.genome_strand[position:position+e_len]
                self.positions.append((position,position+e_len))
                position += e_len + i_len

        self.accession = '>spliced_variant{0},rc={1}'.format(self.positions,self.reverse_complement)


    def __str__(self):
        return '{0}\n{1}\n'.format(self.accession, self.sequence)
        

class DNAseq(object):
    """docstring for DNAseq"""
    def __init__(self,read_length,coverage, mean=None,stddev=None, min_size=None, max_size = None, contamination_rate = 0,contamine_mean= None,contamine_std_dev= None, distribution='normal'):
        super(DNAseq, self).__init__()
        self.distribution = distribution
        
        if self.distribution == 'normal':
            self.mean = mean
            self.stddev = stddev
        elif self.distribution == 'uniform':
            self.min_size = min_size
            self.max_size = max_size

        self.read_length = read_length
        self.contamination_rate = contamination_rate
        self.contamine_mean = contamine_mean
        self.contamine_std_dev = contamine_std_dev
        self.coverage = coverage
    
    def simulate_pe_reads(self,genome):
        """
        Arguments:

        """
        genome_length = len(genome.sequence)
        number_of_reads=(genome_length*self.coverage)/(2*self.read_length)     #Specifiels the number of simulated read pairs (related to insertion size length of genome and coverage
    
        self.reads = []
        
        
        for i in range(number_of_reads):
            if self.distribution == 'normal':
                read_pair = PairedEndRead(distribution = self.distribution, mean=self.mean,sigma=self.stddev,read_length=self.read_length)
            elif self.distribution == 'uniform':
                read_pair = PairedEndRead(distribution = self.distribution, min_size=self.min_size,max_size=self.max_size,read_length=self.read_length)
            
            read_pair.generate(genome.accession.replace(" ",""), genome.sequence, i)        
            self.reads.append(read_pair)

    def simulate_mp_reads(self,genome):
        """
        Arguments:

        """
        genome_length = len(genome.sequence)
        number_of_reads=(genome_length*self.coverage)/(2*self.read_length)     #Specifiels the number of simulated read pairs (related to insertion size length of genome and coverage
    
        self.reads = []
        
        
        for i in range(number_of_reads):
            r = random.uniform(0,1)
            if r < self.contamination_rate:
                # contamination distribution is allways normally distribbuted in this implementation
                read_pair = PairedEndRead(mean=self.contamine_mean, sigma=self.contamine_std_dev, read_length=self.read_length)
                read_pair.generate(genome.accession.replace(" ",""), genome.sequence, i)

            else:
                if self.distribution == 'normal':
                    read_pair = MatePairRead(distribution = self.distribution, mean=self.mean,sigma=self.stddev,read_length=self.read_length)
                elif self.distribution == 'uniform':
                    read_pair = MatePairRead(distribution = self.distribution, min_size=self.min_size, max_size=self.max_size,read_length=self.read_length)
                
                read_pair.generate(genome.accession.replace(" ",""), genome.sequence, i)
            
            self.reads.append(read_pair)


    def fastq_format(self):
        for pe_read in self.reads:
            for mate in pe_read.fastq_format():
                yield mate

    def fasta_format(self):
        for pe_read in self.reads:
            for mate in pe_read.fasta_format():
                yield mate



class RNAseq(object):
    """docstring for RNAseq"""
    def __init__(self,lib_mean,lib_std_dev,lib_read_length):
        super(RNAseq, self).__init__()

        self.lib_mean = lib_mean
        self.lib_std_dev = lib_std_dev
        self.lib_read_length = lib_read_length

    def simulate_pe_reads(self,transcriptome,coverage):
        """
        Arguments:

        transcriptome   - A Transcriptome object
        """
        transcriptome_length = reduce(lambda x,y: x+y, [len(transcript.sequence) for transcript in transcriptome.transcripts]) # total transcriptome length
        number_of_reads=(transcriptome_length*coverage)/(2*self.lib_read_length)     #Specifiels the number of simulated read pairs (related to insertion size length of genome and coverage
        
        self.reads = []
        
        
        for i in range(number_of_reads):
            transcript_index = int(random.randrange(len(transcriptome.transcripts)))
            transcript = transcriptome.transcripts[transcript_index]
            pe = PairedEndRead()
            pe.generate(transcript.accession[1:].replace(" ",""), transcript.sequence, i, self.lib_mean,self.lib_std_dev,self.lib_read_length)
            self.reads.append(pe)

    def fastq_format(self):
        for pe_read in self.reads:
            for mate in pe_read.fastq_format():
                yield mate

    def fasta_format(self):
        for pe_read in self.reads:
            for mate in pe_read.fasta_format():
                yield mate




if __name__=='__main__':
    pass