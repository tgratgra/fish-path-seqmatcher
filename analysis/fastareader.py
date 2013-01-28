#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Reading result file output by Fasta.

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

import re

#from relais.dev.io.readers import singlereader


## CONSTANTS & DEFINES ###

COMPARED_RE = re.compile ("[0-9]+\s+nt\s+overlap")
IDENTITY_RE = re.compile ("[0-9]+.[0-9]*.\s+identity")
ZSCORE_RE = re.compile (" Z\-score: (\d+.\d*) ")
BITS_RE = re.compile (" bits: (\d+.\d*) ")
NAME_RE = re.compile ("\(+[0-9]+\s+nt\)+$")
NUM_RE = re.compile("(?P<int>\d+)\.*(\d*)")


### IMPLEMENTATION ###

class FastaMatch (object):
	"""
	A set of results returned from a Fasta search.
	
	This is just a handy container, with little extended functionality. 'frac'
	(fraction) members are floats from 0.0 to 1.0 and indicate proportions.
	'res' (residue) members are integers, range from 0 up and indicate a number
	of residues.
	
	"""
	def __init__ (self):
		self.seq_id = None
		self.frac_identity = None
		self.res_compared = None
		self.res_ambigs= None
		self.zscore = None
		self.bits = None
		self.seq = None
		
	def _get_frac_difference (self):
		return (1.0 - self.frac_identity)

	def _get_res_matched (self):
		return int (round (self.res_compared * self.frac_identity))
	
	frac_difference = property (_get_frac_difference)
	res_matched = property (_get_res_matched)


class FastaReader (object):
	def __init__ (self, hndl, mode='rb'):
		if (isinstance (hndl, basestring)):
			hndl = open (hndl, mode)
			hndl_opened = True
		else:
			hndl_opened = False
		self.hndl = hndl
		
	def read (self):
		"""
		Return the parsed date from the Fasta output file.
		"""
		# read and split records
		# as '>>' appears in the header, we split on the line start
		buf = self.hndl.read()
		seqmatchs = buf.split ('\n>>')[1:]
		# extract info from each record
		raw_res = []
		for rec in seqmatchs:
			raw_res.append (self.extract_match (rec))
		return raw_res
			
	def extract_match (self, rec):
		"""
		Retrieve sequence id, overlap and identity from a segment of Fasta file.
		
		(no. of nt compared)
		and the percentage identity. The latter two will be obtained by matching
		with regex, which will yield the raw subsection of relevant text.
		"""
		new_fmatch = FastaMatch()
		
		reclines = rec.split ('\n')
		new_fmatch.seq_id = self.extract_id (reclines[0])
		for item in reclines[1:]:
			compare_match = COMPARED_RE.search (item)
			if compare_match:
				new_fmatch.res_compared = int (
					self.extract_identity (compare_match.group()))
			identity_match = IDENTITY_RE.search (item)
			if identity_match:
				new_fmatch.frac_identity = float (
					self.extract_identity (identity_match.group())) / 100.0
			zscore_match = ZSCORE_RE.search (item)
			if zscore_match:
				new_fmatch.zscore = float (zscore_match.group (1))
			bits_match = IDENTITY_RE.search (item)
			if bits_match:
				new_fmatch.bits = float (
					self.extract_identity (bits_match.group()))
		
		## Return:
		return new_fmatch
		
	def extract_id (self, rawname):
		"""
		Strip out the sequence name.
		
		Removes the sequence length and trailing spaces that the fasta program
		appends to the sequence name.
		"""
		return NAME_RE.sub ('', rawname).rstrip()
	
	def extract_identity (self, line):
		"""
		Return the raw identity (or compared) lines as integers.
		"""
		matchnum = NUM_RE.search (line)
		if matchnum:
			return matchnum.group()
		return None


### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod()


### MAIN ###

if __name__ == '__main__':
	_doctest()


### END ######################################################################
		
