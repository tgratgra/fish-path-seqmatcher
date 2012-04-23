#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Phylogenetic reconstruction via the qjoin commandline program.
"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

import os

from Bio import AlignIO

import clineapp
import scratchfile


__all__ = [
	'QJoinCline',
]


## CONSTANTS & DEFINES ###

INALIGN_NAME = 'in.sth'
OUTTREE_NAME = 'out.txt'


### IMPLEMENTATION ###

class QjoinCline (clineapp.ClineApp):
	"""
	A class for calling qjoin for phylogenetic reconstruction.
	
	"""
	def __init__ (self, exepath='/usr/local/bin/qjoin'):
		clineapp.ClineApp.__init__ (self, exepath, use_workdir=True,
			remove_workdir=False, check_requirements=False)
			
	def setup_workdir (self):
		"""
		Prepare the necessary input files for Qjoin.

		This creates a temporary working area, and writes input alignment
		files.

		"""
		## Preconditions & preparations:
		## Main:
		# create workdir and filepaths
		clineapp.ClineApp.setup_workdir (self)
		self._inalign_path = scratchfile.make_scratch_file (INALIGN_NAME,
			self._curr_workdir)
		# write infile workfile
		#MSG (self._curr_workdir, self._inalign_path)
		infile_hndl = open (self._inalign_path, 'w')
		AlignIO.write ([self._in_align], infile_hndl, 'stockholm')

	def run (self, align, num_bootstraps=0):
		"""
		Run Qjoin to reconstruct a tree from an alignment.

		:Params:
			align : Biopython alignment
				The sequence alignment to be built into a tree.
			num_bootstraps : integer or False
				The number of bootstraps to run. If False (the default), no
				bootstraps are done.

		"""
		## Preconditions & preparation:
		self.set_input_alignment (align)
		## Main:
		cmdopts = [
			INALIGN_NAME,
			OUTTREE_NAME,
		]
		#if (num_bootstraps):
		#	cmdopts.append ('--bootstrap=%s' % num_bootstraps) 
		
		self.call_cmdline (*cmdopts)
		
	def extract_results (self):
		"""
		Obtain the output produced by Qjoin.

		We call this as a seperate function, so the caller has a chance to 
		check the status and error output first.
		
		:Returns:
			A ptree object.
		
		"""
		## Preconditions:
		# make sure that cline has actually run & output exists
		assert (self._curr_cline)
		output_path = os.path.join (self._curr_workdir, OUTTREE_NAME)
		assert (os.path.exists (output_path)), \
			"can't find outfile %s" % output_path
		## Main:
		# extract the data
		output_hndl = open (output_path, 'rU')
		tree_str = output_hndl.read()
		output_hndl.close()

		tmp_tree = tree_str.split('\n\n')[0]
		dist_tree = tmp_tree.split(':', 1)[1]
		## Postconditions:
		return dist_tree
			
			
	def set_input_alignment (self, align):
		"""
		Record the input alignment as a series of SeqRecords.

		This is necessary since the IO functions expect a series of SeqRecords.
		This should be called from the ``run`` command for each class to 
		convert the input data. It stores the result on the member
		``_in_align``.

		:Params:
			align
				Either a BioPython alignment or a sequence of BioPython SeqRecords.

		"""
		self._in_align = align
			
			

### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	_doctest()


### END ######################################################################
