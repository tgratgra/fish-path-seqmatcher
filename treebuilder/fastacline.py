#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Comparing sequences via the FASTA program.


"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

import textwrap, os

from Bio import SeqRecord

from relais.dev.common import *
from relais.dev import clineapp, scratchfile

import fastareader


## CONSTANTS & DEFINES ###

QUERYFILE_NAME = 'query.fasta'
REFFILE_NAME = 'ref.fasta'
OUTPUTFILE_NAME = 'fasta.out'

SEQID = 0


### IMPLEMENTATION ###

def get_new_seqid ():
	"""
	Generate a new sequence id if there are collisions.
	"""
	global SEQID
	SEQID += 1
	return 'refseq_%d' % SEQID


def count_ambigs (seq):
	"""
	Count the number of ambiguties in a DNA sequence.
	
	There's probably a better way to do this, but this will do with the moment.
	Note that we assume the seq is lower case.
	"""
	MSG (seq)
	count = 0
	for base in seq:
		if base not in 'acgt':
			count = count + 1
	return count
	

class FastaCline (clineapp.ClineApp):
	"""
	A class for calling a FASTA similarity comparison on biosequences.
	
	def __init__ (self, exepath='/Users/tgra/Sites/iah2/fasta-35.4.7/bin/fasta35'):
	def __init__ (self, exepath='/home/tgray/apps/fasta-35.4.11/bin/fasta35'):
	"""
	
	def __init__ (self, exepath='/var/www/vhosts/symantix.eu/tmp/fasta-35.4.11/bin/fasta35'):
		
		clineapp.ClineApp.__init__ (self, exepath, use_workdir=True,
			remove_workdir=False, check_requirements=False)
		self._orig_ids_seen = []
		self._new_ids_to_orig_seq= {}
		self._new_ids_to_ambigs= {}
						
	def run (self, target_seq, ref_seqs, maxmatch=20, wordlen=6):
		"""
		Run a FASTA comparison, looking for one sequence in a group of another.
		
		This function accepts Biopython SeqRecords.
		
		:Params:
			target_seq
				The sequence to be matched
			ref_seqs
				The group of sequences to be searched.
			maxmatch
				The maximum number of matches (comparisons) to be returned.
			wordlen
				The wordlength (or ktup) to be used in the comparision.
			
		:Returns:
			to be decided
		
		"""
		# NOTE: SK used '<queryfile> <reffile> <worklen> ktup -O <outfile>
		# -d <maxreturn> -Q' as the commandline but this doesn't work, least not
		# under fasta35. '-Q -O <outfile>' have to go first.
		# NOTE: '-Q' is quiet (no user interaction) and ktup number comes first.
		
		## Preconditions:
		# list what ids have been seen
		self._orig_ids_seen = []
		# map the new ids to original sequences
		self._new_ids_to_orig_seq = {}
		# maps new ids to the number of ambigs
		self._new_ids_to_ambigs = {}
		## Main:
		# record seqs
		self._target_seq = target_seq
		self._ref_seqs = ref_seqs
		cmdfmt = '-Q  -O %(out)s %(qry)s %(ref)s %(wlen)d ktup -d %(max)d'
		self.call_cmdline (cmdfmt % {
			'qry': QUERYFILE_NAME,
			'ref': REFFILE_NAME,
			'wlen': wordlen,
			'out': OUTPUTFILE_NAME,
			'max': maxmatch,
		})
				
	def extract_results (self):
		"""
		Obtain the output produced by FASTA.
		
		We call this as a seperate function, so the caller has a chance to 
		check the status and error output first.
		"""
		# NOTE: '_new_ids_to_orig_seq' - map the new ids to original sequences
		# NOTE: '_new_ids_to_ambigs' - maps new ids to the number of ambigs
		
		## Preconditions:
		# make sure that cline has actually run & output exists
		assert (self._curr_cline)
		output_path = os.path.join (self._curr_workdir, OUTPUTFILE_NAME)
		assert (os.path.exists (output_path)), \
			"can't find outfile %s" % output_path
		## Main:
		# extract the data
		reader = fastareader.FastaReader(output_path)
		fasta_results = reader.read()
		for item in fasta_results:
			fastaid = item.seq_id
			item.res_ambigs = self._new_ids_to_ambigs[fastaid]
			item.seq = self._new_ids_to_orig_seq[fastaid]
		## Postconditions:
		return fasta_results
			
				
	def setup_workdir (self):
		"""
		Perpare the necessary input files for fasta.
		
		This creates a temporary working area, and writes files for the
		query and reference sequences.
		
		"""
		## Preconditions & preparations:
		self._orig_ids_seen = []
		self._new_ids_to_orig_seq= {}
		self._new_ids_to_ambigs= {}
		
		# create workdir and filepaths
		clineapp.ClineApp.setup_workdir (self)
		self._queryfile_path, self._reffile_path = scratchfile.make_scratch_files (
			[QUERYFILE_NAME, REFFILE_NAME], self._curr_workdir)
		# write query workfile
		queryfile_hndl = open (self._queryfile_path, 'w')
		self.write_seq_to_file (self._target_seq, queryfile_hndl)
		queryfile_hndl.close()
		# write ref workfile
		reffile_hndl = open (self._reffile_path, 'w')
		for seq in self._ref_seqs:
			self.write_seq_to_file (seq, reffile_hndl)
		reffile_hndl.close()		

	def write_seq_to_file (self, seq, hndl):
		"""
		Write the passed sequences to the file handle in FASTA format.
		
		Note that we record duplicate ids here and chnage them so we can
		translate them back later, as well as recording the number of
		ambiguities in the source sequence.
		"""
		# NOTE: '_orig_ids_seen' - list what ids have been seen
		# NOTE: '_new_ids_to_orig_seq' - map the new ids to original sequences
		# NOTE: '_new_ids_to_ambigs' - maps new ids to the number of ambigs


		
		## Main:
		# extract id and actual sequence
		assert (isinstance (seq, SeqRecord.SeqRecord))
		seqid = seq.id or seq.name
		seqstr = str (seq.seq.tostring())

		# record data, check for duplicate ids
		seqstr = seqstr.lower()
		while (seqid in self._orig_ids_seen):
			seqid = get_new_seqid()
		self._orig_ids_seen.append (seqid)
		self._new_ids_to_orig_seq[seqid] = seq
		self._new_ids_to_ambigs[seqid] = count_ambigs (seqstr)
		# write it
		hndl.write ('>' + seqid + '\n')
		hndl.write (textwrap.fill (seqstr) + '\n')
		
	def extract_diagnostics (self):
		diag = {}
		filenames = [
			QUERYFILE_NAME,
			REFFILE_NAME,
			OUTPUTFILE_NAME,
		]
		for item in filenames:
			fpath = os.path.join (self._curr_workdir, item)
			diag[item] = utils.file_to_string (fpath)
		return diag


### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	_doctest()


### END ######################################################################
