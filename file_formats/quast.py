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
		self.accession = re.search(re.escape(stringA)+"(.*?)"+re.escape(stringB),self.scaffold_text).group(1).strip()
		self.length = int(re.search("\([0-9]+",self.scaffold_text).group(0)[1:])
		#pattern = "CONTIG:(.+)"
		#print self.scaffold_text[:20]
		#print re.search(pattern,self.scaffold_text).group(0)		
	def get_extensive_missassembly(self):
		self.breakpoints = []
		for line1,line2,line3 in zip(self.scaffold_text.splitlines()[:-2],self.scaffold_text.splitlines()[1:-1],self.scaffold_text.splitlines()[2:]):
			#print repr(line2.lstrip()[:4])

			#res = re.findall('misasse',self.scaffold_text)
			#if res:
			#	print res.group(0)
			#if line2.lstrip(" \t\n\r\x0b\x0c")[0:4] == 'Exte':
			#	print repr(line2)
			if line1.lstrip(' \t\n\r')[:4] == 'Real' and line3.lstrip()[:4] == 'Real' and line2.lstrip()[:4] == 'Exte':
				#print line1
				#print line2
				#print line3
				#print 'BP =', max( map(lambda x: int(x), line1.split('|')[1].strip().split(' ')))
				self.breakpoints.append(max( map(lambda x: int(x), line1.split('|')[1].strip().split(' '))))

				#return line2# print line2
		#return 1

		

class ScaffoldContainer(object):
	"""docstring for ScaffoldContainer"""
	def __init__(self):
		super(ScaffoldContainer, self).__init__()
		self.scaffolds = []
	def add_scaffold(self, scaffold_alignment):
		self.scaffolds.append(scaffold_alignment)


def quast_parse_stdout(infile):
	container = ScaffoldContainer()
	k= 0
	scaffold_string = ''
	for line in infile:
		if line[:7] == 'CONTIG:' and k==0:
			k = 1
			scaffold_string = line
		elif line[:7] == 'CONTIG:':
			if scaffold_string[:7] == 'CONTIG:' and k==1:
				scaf = ScaffoldAlignment(scaffold_string)
				scaf.get_info()
				container.add_scaffold(scaf)
				scaffold_string =line
		else:
			scaffold_string += line




	scaf = ScaffoldAlignment(scaffold_string)
	scaf.get_info()
	container.add_scaffold(scaf)

	print len(container.scaffolds)
	for scaf in container.scaffolds:
		scaf.get_data()
		scaf.get_extensive_missassembly()
		if scaf.breakpoints:
			#print scaf
			print scaf.accession
			for bp in scaf.breakpoints:
				print bp

			#print scaf.accession, scaf.length
			#print repr(scaf.scaffold_text)
		#	i+=1
		#if i>340:
		#	break

			#print scaf.get_extensive_missassembly()
	#scaffold_alignment = "CONTIG: (.+) CONTIG:"
	#print re.findall(scaffold_alignment,text,re.MULTILINE | re.DOTALL)
	#for scaffold_text in re.findall(scaffold_alignment,text):
	#	print scaffold_text
