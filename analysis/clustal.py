#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SHORT DESCRIPTION.

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

from Bio import SeqIO, Clustalw

from relais.dev import scratchfile

__all__ = [
	'align_with_clustal',
]


### CONSTANTS & DEFINES ###

### IMPLEMENTATION ###	

def align_with_clustal (bp_seqs):
	"""
	Align a set of Biopython sequences with Clustal.
	
	:Parameters:
		bp_seqs
			An iterable of Biopython SeqRecords.
			
	:Returns:
		A Biopython multiple alignment.
		
	This presents a friendly wrapper around the Biopython clustal machinery,
	obscuring it's oddities (e.g. needing a physical input file instead of 
	just a group of sequences). It runs the clustal task in a temporary
	directory and returns the result.
	
	Note: this function cannot accept Biopython Seqs, because they have no id
	or name.
	
	"""
	tmpdir = scratchfile.make_scratch_dir()
	infile = scratchfile.make_scratch_file ('in.fasta', tmpdir)
	infile_hndl = open (infile, "w")
	SeqIO.write (bp_seqs, infile_hndl, "fasta")
	infile_hndl.close()
	cmdline = Clustalw.MultipleAlignCL (infile)
	cmdline.set_output ('out.aln')
	return Clustalw.do_alignment (cmdline)


	
### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	main()


### END ######################################################################
