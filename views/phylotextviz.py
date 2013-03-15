#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Visualise a phylogeny as simple text (a Newick tree).

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

from htmltags import *
import baseqryviz


__all__ = [
	'PhyloTextViz',
]


### CONSTANTS & DEFINES ###

### IMPLEMENTATION ###	
		
class PhyloTextViz (baseqryviz.BaseQryViz):
	id = __module__ + '.PhyloTextViz'.lower()
	label='Text representation of phylogeny'
	resources = []
	
	def render (self):
		# NOTE: data is an ete2 tree so we do things differently
		return tag_with_contents ('textarea',
			self.data.write(),
			class_="revi_treestring"
		)
	


### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	_doctest()


### END ######################################################################
