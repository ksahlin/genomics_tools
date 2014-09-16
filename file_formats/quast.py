import re

class ScaffoldAlignment(object):
	"""docstring for ScaffoldAlignemt"""
	def __init__(self, scaffold_text):
		super(ScaffoldAlignment, self).__init__()
		self.scaffold_text = scaffold_text
	def get_info(self):
		pass

	def get_data(self):
		stringA = "CONTIG:"
		stringB = "("
		self.accession = re.search(re.escape(stringA)+"(.*?)"+re.escape(stringB),self.scaffold_text).group(1)
		self.length = int(re.search("\([0-9]+",self.scaffold_text).group(0)[1:])
		#pattern = "CONTIG:(.+)"
		#print self.scaffold_text[:20]
		#print re.search(pattern,self.scaffold_text).group(0)		
	def get_extensive_missassembly(self):
		for line1,line2,line3 in zip(self.scaffold_text.splitlines()[:-2],self.scaffold_text.splitlines()[1:-1],self.scaffold_text.splitlines()[2:]):
			print line2.lstrip()[:4]
			if line2.lstrip()[:4] == 'Exte':
				print line2
			if line1.lstrip()[:4] == 'Real' and line3.lstrip()[:4] == 'Real' and line2.lstrip()[:4] == 'Exte':
				return line2# print line2

		

class ScaffoldContainer(object):
	"""docstring for ScaffoldContainer"""
	def __init__(self):
		super(ScaffoldContainer, self).__init__()
		self.scaffolds = []
	def add_scaffold(self, scaffold_alignment):
		self.scaffolds.append(scaffold_alignment)


def quast_parse_stdout(infile):
	container = ScaffoldContainer()
	scaffold_string = ''
	for line in infile:
		if line[:7] == 'CONTIG:':
			scaffold_string = line 
			if scaffold_string[:7] == 'CONTIG:':
				scaf = ScaffoldAlignment(scaffold_string)
				scaf.get_info()
				container.add_scaffold(scaf)

		else:
			scaffold_string += line

	scaf = ScaffoldAlignment(scaffold_string)
	scaf.get_info()
	container.add_scaffold(scaf)

	print container
	for scaf in container.scaffolds:
		if scaf.get_extensive_missassembly():
			print scaf
			scaf.get_data()
			print scaf.accession, scaf.length
			print scaf.get_extensive_missassembly()
	#scaffold_alignment = "CONTIG: (.+) CONTIG:"
	#print re.findall(scaffold_alignment,text,re.MULTILINE | re.DOTALL)
	#for scaffold_text in re.findall(scaffold_alignment,text):
	#	print scaffold_text
