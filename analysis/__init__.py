#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Various ways of making matches or building trees from seqs.

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

import re
import string

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import fastacline
import mafft
import qjoincline


## CONSTANTS & DEFINES ###

SEQNAME_RE = re.compile (r'[%s]' % re.escape (string.punctuation))
SPACE_RE = re.compile (r'\s+')


### IMPLEMENTATION ###

def clean_seqname (n):
	n = SEQNAME_RE.sub (' ', n)
	n = SPACE_RE.sub ('_', n.strip())
	return n


def db_seq_to_seqrec (rec, sel_regions):
	"""
	Convert a database record to a BioPython SeqRecord.
	
	There's some mighty odd magic in the gene_region field of the sequences.
	I don't even pretend to know what it is or does, I'm just copying Tanya's
	lead.
	"""
	
	#sensible_name = rec['accession_number'] or rec['isolate_name'] or rec['id'] \
	#	or '(unlabelled)'
	sensible_name = "%s (%s)" % (rec['isolate_name'], rec['accession_number'])
	seq = rec['nucleotide_sequence']
	
	# gene_region seems to encode the actual region to be extracted
	gene_region = rec['gene_region']
	
	start = end = None
	for region in gene_region.split('|'):
		array_region = region.split(";")
		if (2 < len(array_region)>2):
			curr_region = array_region[2]
			if (curr_region in sel_regions):
				start = int (array_region[0]) - 1
				end = int (array_region[1]) - 1
				break
	subseq = rec['nucleotide_sequence'][start:end]
						 
	return SeqRecord (id=sensible_name, name=sensible_name, seq=Seq(subseq))
	
	
	
	
def match_seqs (method, selected_genes, sel_regions, target_seq, exepaths={}):
	## Preconditions & preparation:
	## Main:
	# convert everything to Seqrecs
	bioseqs = [db_seq_to_seqrec (g, sel_regions) for g in selected_genes]
	if target_seq:
		bioseqs. append (SeqRecord (id='QUERY', name='QUERY', seq=Seq (target_seq)))
	
	# align and build tree
	results = []
	
	if (method == 'nj'):
		# clean up names for later display
		for b in bioseqs:
			b.id = clean_seqname (b.id)
			b.name = clean_seqname (b.name)
			
		# now run
		mafft_app = mafft.MafftCline (exepath=exepaths['mafft'])
		mafft_app.run_fftns (bioseqs)
		bp_alignment = mafft_app.extract_results()
		cli = qjoincline.QjoinCline (exepath=exepaths['qjoin'])
		cli.run (bp_alignment)
		phyl = cli.extract_results()
		results.append (('newicktree', phyl))
		
	elif (method == 'fasta'):
		fcline = fastacline.FastaCline (exepath=exepaths['fasta'])
		# because the target seq is last, nominate it as target
		fcline.run (bioseqs[-1], bioseqs[:-1])
		f_results = fcline.extract_results()
		results.append (('fastamatchs', f_results))
		
	else:
		raise StandardError, "unknown matching method '%s'" % method
	
	## Postconditions & return:
	return results
	
	
### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	_doctest()


### END ####################################################################

