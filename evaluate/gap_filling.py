import argparse
import os

from genomics_tools.file_formats import quast


def parse_GapFiller(infile_):
	pass		
def parse_GapCloser(infile_):
	pass

def parse_Salmela(infile_):
	pass

def run_quast

def main(args):
	if not os.path.exists(args.outpath):
		os.mkdir(args.outpath)
	quast.quast_parse_stdout(open(args.scaffolds,'r'))



if __name__ == '__main__':
    ##
    # Take care of input
	parser = argparse.ArgumentParser(description="Evaluate gaps filled by the Salmela-Ford algorithm.")
	parser.add_argument('scaffolds', type=str, help='A quast contig_reports.stdout.')
	parser.add_argument('outpath', type=str, help='Folder for output.')
	# parser.add_argument('sigma', type=float, help='Stddev of library')
	# parser.add_argument('cov', type=float, help='Mean coverage of library')
	# parser.add_argument('r', type=float, help='read length  of library')
	# parser.add_argument('s', type=float, help='Maximum allowed softclipped bases')

	args = parser.parse_args()
	main(args)