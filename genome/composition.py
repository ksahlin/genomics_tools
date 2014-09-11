import os
import argparse

from matplotlib import pyplot

from genomics_tools.file_formats import fasta
from collections import defaultdict,Counter

def batch_gen(data, batch_size):
    for i in range(0, len(data), batch_size):
    	if i + batch_size <= len(data):
            yield data[i:i+batch_size]
        else:
        	raise StopIteration

class KmerInfo(object):
	"""docstring for KmerInfo
	 distance_histogram - A dictionary with distance (int) as key
	 and number of occurences of that distance as value. Distance
	 is the distance between closest neighbors of tha same kmer """
	def __init__(self):
		super(KmerInfo, self).__init__()

		self.distance_histogram = Counter() 
		self.frequency = 1
		self.current_distance = 0

	def update(self):
		"""
			When the kmer is seen on the genome
		"""
		self.frequency += 1
		self.distance_histogram[self.current_distance] += 1
		self.current_distance = 0


		

class DistanceContainer(object):
	"""docstring for DistanceContainer"""
	def __init__(self):
		super(DistanceContainer, self).__init__()
		self.kmers = defaultdict(KmerInfo)

	def parse_genome(self,genome, k):
		for kmer in batch_gen(genome, k):
			# update distance for all seen kmers 
			for seen_kmer in self.kmers:
				self.kmers[seen_kmer].current_distance += 1
				

			# special treatment for the current kmer	
			if kmer not in self.kmers:
				self.kmers[kmer] = KmerInfo()
			else:
				self.kmers[kmer].update()

	def plot(self):
		tot_kmers = sum([self.kmers[kmer].frequency for kmer in self.kmers] ) 
		for kmer in self.kmers:
			ratio = self.kmers[kmer].frequency / float(tot_kmers)
			# print kmer, ratio, self.kmers[kmer].distance_histogram
			distance, frequency = zip(*sorted(self.kmers[kmer].distance_histogram.items()))
			#distance = self.kmers[kmer].distance_histogram.keys()
			#frequency = self.kmers[kmer].distance_histogram.values()
			pyplot.plot(distance, frequency, label=kmer)
			pyplot.savefig(os.path.join(args.folder_path,kmer+'-'+str(ratio)+'.pdf'))
			e_x = sum([(i+1)*frequency[i] for i in range(len(frequency))]) / float(sum(frequency))
			print 'mean:', e_x




		
def main(args):
	if not os.path.exists(args.folder_path):
		os.mkdir(args.folder_path)
	for acc, seq in fasta.fasta_iter(open(args.genome,'r')):
		accession,genome = acc, seq
		break

	#print seq
	dist = DistanceContainer()
	dist.parse_genome(seq,args.k)
	dist.plot()



if __name__ == '__main__':
    ##
    # Take care of input
	parser = argparse.ArgumentParser(description="Plots the distance between identical closest k-mers in a genome.")
	parser.add_argument('k', type=int, help='k-mer size ')
	parser.add_argument('genome', type=str, help='Genome in fasta format. ')
	parser.add_argument('folder_path', type=str, help='Location of output')
	# parser.add_argument('sigma', type=float, help='Stddev of library')
	# parser.add_argument('cov', type=float, help='Mean coverage of library')
	# parser.add_argument('r', type=float, help='read length  of library')
	# parser.add_argument('s', type=float, help='Maximum allowed softclipped bases')

	args = parser.parse_args()
	main(args)
