import sys,os,subprocess

import argparse
from genomics_tools.file_formats import bam
from genomics_tools.simulate import genome,contigs,reads
from mapping import align
from genomics_tools.file_formats.various_annotations import to_AGP,to_GFF

def simulate_instance(args):
    print 'Started simulating'
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)
    genome_path = os.path.join(args.output_path, 'genome.fa')
    contig_path = os.path.join(args.output_path, 'ctgs.fa')
    read1_path = os.path.join(args.output_path, 'reads1.fa')
    read2_path = os.path.join(args.output_path, 'reads2.fa')
    bam_path = os.path.join(args.output_path, 'mapped')
    gff_path = os.path.join(args.output_path, 'true_error_pos.gff')
    gff_file = open(gff_path,'w')

    #genome
    genomelen = args.burnin + ( (args.contiglen+args.gaplen)*(args.nrgaps + 1 ) + args.contiglen ) * (len(args.errorsize) + 1)
    print genomelen
    g = genome.Genome([0.25]*4,genomelen,'genome1')
    g.genome()
    print >> open(genome_path,'w'), g.genome_fasta_format()

    #contigs/scaffolds
    if args.scaffolds:
    	scafs = open(contig_path,'w')
        scafs.write('>scf_burnin{0}\n{1}\n'.format(args.gaplen,g.sequence[0:args.burnin]))
    	scaffold = ''
        pos = args.burnin

        for error in args.errorsize:
            scaffold_coord = 0
            for i,x in enumerate(range(pos, pos + (args.nrgaps + 1)*(args.contiglen + args.gaplen ), args.contiglen + args.gaplen)):
                #print 'pos:', x
                if (args.gaplen + error) > 0:
                    if i < args.nrgaps:
                        scaffold += g.sequence[x:x+args.contiglen]+ 'N'* (args.gaplen + error) 
                        scaffold_coord = len(scaffold)
                        error_start = scaffold_coord - (args.gaplen + error) 
                        error_stop = scaffold_coord  # error is anywhere in the introduced gap (either contraction or expansion)
                    else:
                        scaffold += g.sequence[x:x+args.contiglen]
                else:
                    #scaffold += g.sequence[i*(args.gaplen + error) + x : x + args.contiglen + (i+1)*(args.gaplen + error)] 
                    scaffold += g.sequence[x : x + args.contiglen + (args.gaplen + error)] 

                    scaffold_coord = len(scaffold)
                    error_start = scaffold_coord
                    error_stop = scaffold_coord+1 # error is at a specific position where a contraction has occured

                if error < 0 and i < args.nrgaps:
                    to_GFF(gff_file, 'scf_gap{1}_errorsize_minus{2}'.format(i+1, args.gaplen, abs(error)), 'TRUTH','FCD', error_start, error_stop, 1, '+', '.', 'Note=Error:Contraction {0}bp'.format(abs(error)))
                elif error > 0 and i < args.nrgaps:
                    to_GFF(gff_file, 'scf_gap{1}_errorsize{2}'.format(i+1, args.gaplen, abs(error)), 'TRUTH','FCD', error_start, error_stop, 1, '+', '.', 'Note=Error:Expansion {0}bp'.format(abs(error)))
                else:
                    pass

            if error <0:
                scafs.write('>scf_gap{1}_errorsize_minus{2}\n{3}\n'.format(i+1, args.gaplen, abs(error), scaffold)) 
            else:
                scafs.write('>scf_gap{1}_errorsize{2}\n{3}\n'.format(i+1, args.gaplen, error, scaffold))   
	
            scaffold = ''
            pos = x + 2*args.contiglen  
        # dummy sequences to prevent bwa tor remove any of our scaffolds
        for i in range(10):
            dummy = genome.Genome([0.25]*4,1000,'z_dummy{0}'.format(i+1))
            dummy.genome()
            scafs.write('>z_dummy{0}\n{1}\n'.format(i+1, dummy.sequence)) 
            
    else:
    	ctgs = open(contig_path,'w')
        ctgs.write('>ctg0\n{0}\n'.format(g.sequence[0:args.burnin]))
    	for i,x in enumerate(range(args.burnin,genomelen,(args.contiglen + args.gaplen))):
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
    #align.bwa_mem(read1_path, read2_path, contig_path, bam_path, args)

def main(args):
    successful_experiments = 0
    while successful_experiments < 1: 
        try:
            simulate_instance(args)
        except subprocess.CalledProcessError:
            continue

        successful_experiments += 1
	

if __name__ == '__main__':
    ##
    # Take care of input
    parser = argparse.ArgumentParser(description="Simulate paired end reads with a given mean and standard deviation on mean insert size. A given coverage should also be provided.")
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
    parser.add_argument( '-errors', dest='errorsize', type=int, nargs='+', default=False, help='gap distance error' )
    parser.add_argument( '-burnin', dest='burnin', type=int, default=False, help='Estimation window' )
    parser.add_argument( '-nrgaps', dest='nrgaps', type=int, default=False, help='Number of gaps' )
    parser.add_argument('--threads', type=str, dest='threads', default='8', required=False, help='Number of threads for bwa mem.')


   
    #parser.add_argument('coverage', type=int, help='Coverage for read library. ')
    #parser.add_argument('outpath', type=str, help='Path to output location. ')
    #parser.add_argument('experiments', type=int, help='Number of experiment for each line in sim_in.txt file to run. ')
    #parser.add_argument('genomelen', type=int, help='Length of the reference sequence. ')
    #parser.add_argument('c_len', type=int, help='Contig length. ')

    args = parser.parse_args()
    main(args)
