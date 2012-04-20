#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Various ways of making matches or building trees from seqs.

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

## CONSTANTS & DEFINES ###

### IMPLEMENTATION ###

def db_seq_to_seqrec (rec):
	"""
	Convert a database record to a BioPython SeqRecord.
	
	There's some mighty odd magic in the gene_region field of the sequences.
	I don't even pretend to know what it is or does, I'm just copying Tanya's
	lead.
	"""
	
	#sensible_name = rec['accession_number'] or rec['isolate_name'] or rec['id'] \
	#	or '(unlabelled)'
	sensible_name = "%s (%s)" % (rec['isolate_name'], rec['accession_number'])
	seq = rec['nucleotide']
	
	# gene_region seems to encode the actual region to be extracted
	gene_region = rec['gene_region']
	
	start = end = None
	for region in gene_region.split('|'):
		array_region = region.split(";")
		if (2 < len(array_region)>2):
			thisgeneregion = array_region[2]
			if (thisgeneregion == selectedgeneregion):
				start = int(array_region[0]) - 1
				end = int(array_region[1]) - 1
				break
	subseq = sequence[start:end]
						 
	return SeqRecord (id=sensible_name, name=sensible_name, seq=Seq(subseq))
	
	
	
	
def build_tree (method, selected_genes, target_seq):
	## Preconditions & preparation:
	
	bioseqs = [db_seq_to_seqrec (g) for g in selected_genes]
	if target_seq:
		bioseqs.append (SeqRecord (id='QUERY', name='QUERY', seq=Seq (bioseq))
	
	# align and build tree
	results = []
	
	if (method == 'nj'):
		mafft_app = mafft.MafftCline()
		mafft_app.run_fftns (bioseqs)
		bp_alignment = mafft_app.extract_results()
		cli = qjoincline.QjoinCline()
		cli.run (bp_alignment)
		phyl = cli.extract_results()
		results.append (('tree', phyl))
		
	elif (method == 'fasta'):
		fcline = fastacline.FastaCline()
		fcline.run (bioseqs[-1], bioseqs[:-1])
		f_results = fcline.extract_results()
		results.append (('table', f_results))
	else:
		raise StandardError, "unknown matching method '%s'" % method
	
	
### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	_doctest()


### END ####################################################################

