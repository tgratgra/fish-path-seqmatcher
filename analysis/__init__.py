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
import fasttreecline
import config


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
	
	
	
	
def match_seqs (method, selected_genes, sel_regions, target_seq, outgroup, exepaths={}):
	"""
	Return similarities between a set of sequences, as a table or tree.

	:Parameters:
		method
			a string giving how to match seqs
		selected_genes
			database ids for the sequences involved
		sel_regions
			required for the db lookup
		target_seq
			the raw seq data for a query sequence
		outgroup
			the id of the sequence to be set as outgroup in a tree
		exe_paths
			hash of paths to the various programs used

	:Returns:
		Either a Fasta table or an ETE tree

	This should really have been two seperate functions, one for tables and
	one for trees. Now it's a frozen mistake.
	"""
	## Preconditions & preparation:

	## Main:
	# gather ref seqs as seqrecs, record outgroup
	bioseqs = []
	out_seqrec = None
	for g in selected_genes:
		sr = db_seq_to_seqrec (g, sel_regions)
		bioseqs.append (sr)
		if outgroup and (g['id'] == int(outgroup)):
			out_seqrec = sr
	# if target seq, add it too
	if target_seq:
		bioseqs.append (SeqRecord (id='QUERY', name='QUERY', seq=Seq (target_seq)))
	
	results = []
	
	if (method == 'nj'):
		# align and build tree
		# clean up names for later display
		for b in bioseqs:
			b.id = clean_seqname (b.id)
			b.name = clean_seqname (b.name)
			
		# align seqs
		mafft_app = mafft.MafftCline (exepath=exepaths['mafft'])
		mafft_app.run_fftns (bioseqs)
		bp_alignment = mafft_app.extract_results()

		# build tree
		# TODO: number of bootstraps should be config val
		cli = qjoincline.QjoinCline (exepath=exepaths['qjoin'])
		cli.run (bp_alignment, num_bootstraps=config.NUM_BOOTSTRAPS)
		phyl = cli.extract_results()

		#cli = fasttreecline.FasttreeCline (exepath=exepaths['fasttree'])
		#cli.run (bp_alignment)
		#phyl = cli.extract_results()

		# now reroot tree if required
		if outgroup:
			assert out_seqrec, "should have outgroup seqrec not None"
			outgroup_name = out_seqrec.name
			possible_outgroup_names = [outgroup_name, "'%s'" % outgroup_name]
			for n in phyl.traverse():
				if n.name in possible_outgroup_names:
					og_node = n
					break
			phyl.set_outgroup (n)

		results.append (('tree', phyl))
		
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

