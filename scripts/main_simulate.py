import sys,os

import argparse
from genomics_tools.file_formats import bam
from genomics_tools.simulate import genome,contigs,reads
from mapping import align

def simulate_instance(args):
    if args.burnin >= args.genomelen:
        print 'Burn-in longer than genome length, exiting...'
        sys.exit()
    print 'Started simulating'
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
    genome_path = os.path.join(args.output_path, 'genome.fa')
    contig_path = os.path.join(args.output_path, 'ctgs.fa')
    read1_path = os.path.join(args.output_path, 'reads1.fa')
    read2_path = os.path.join(args.output_path, 'reads2.fa')
    bam_path = os.path.join(args.output_path, 'mapped')

    #genome
    g = genome.Genome([0.25]*4,args.genomelen,'genome1')
    g.genome()
    print >> open(genome_path,'w'), g.genome_fasta_format()

    #contigs/scaffolds
    if args.scaffolds:
    	scafs = open(contig_path,'w')
        scafs.write('>scf0\n{0}\n'.format(g.sequence[0:args.burnin]))
    	scaffold = ''
    	for i,x in enumerate(range(args.burnin, args.genomelen,(args.contiglen + args.gaplen))):
    		scaffold += g.sequence[x:x+args.contiglen]+ 'N'* (args.gaplen-args.error)
        scafs.write('>scf{0}\n{1}\n'.format(i+1,scaffold))   	
    else:
    	ctgs = open(contig_path,'w')
        ctgs.write('>ctg0\n{0}\n'.format(g.sequence[0:args.burnin]))
    	for i,x in enumerate(range(args.burnin,args.genomelen,(args.contiglen + args.gaplen))):
        	ctgs.write('>ctg{0}\n{1}\n'.format(i+1,g.sequence[x:x+args.contiglen]))

    #reads
    lib = reads.DNAseq(args.read_length ,args.coverage, args.mean,args.sd)
    lib.simulate_pe_reads(g)
    reads1 = open(read1_path,'w')
    reads2 = open(read2_path,'w')
    i=0
    for read in lib.fasta_format():
        if i%2==0:
            reads1.write(read)
        else:
            reads2.write(read)
        i+=1

    print 'Started mapping'
    #mapping
    align.map_paired_reads(read1_path, read2_path, contig_path, bam_path, args)

def main(args):
	simulate_instance(args)

if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Simulate paired end reads with a given mean and standard deviation on mean insert size. A given coverage should also be provided.")
    parser.add_argument('genomelen', type=int, help='Length of genome. ')
    parser.add_argument('contiglen', type=int, help='Length of contigs. ')
    parser.add_argument('gaplen', type=int, help='Length of gap. ')

    parser.add_argument('mean', type=int, help='Mean insert. ')
    parser.add_argument('sd', type=int, help='Stddev insert.')
    parser.add_argument('coverage', type=int, help='Coverage. ')
    parser.add_argument('read_length', type=int, help='Length of read fragments within a pair. ')
    parser.add_argument('output_path', type=str, help='path to folder output. ')
    parser.add_argument( '-sort', dest='sort', action='store_true', default=False, help='Coordinate sort the reads in the bam file' )
    parser.add_argument( '-sam', dest='sam', action='store_true', default=False, help='Output a samfile (default is bam)' )
    parser.add_argument( '-scafs', dest='scaffolds', action='store_true', default=False, help='scaffolds are simulated instead of contigs' )
    parser.add_argument( '-error', dest='error', type=int, default=False, help='gap distance error' )
    parser.add_argument( '-burnin', dest='burnin', type=int, default=False, help='Estimation window' )

   
    #parser.add_argument('coverage', type=int, help='Coverage for read library. ')
    #parser.add_argument('outpath', type=str, help='Path to output location. ')
    #parser.add_argument('experiments', type=int, help='Number of experiment for each line in sim_in.txt file to run. ')
    #parser.add_argument('genomelen', type=int, help='Length of the reference sequence. ')
    #parser.add_argument('c_len', type=int, help='Contig length. ')

    args = parser.parse_args()
    main(args)
