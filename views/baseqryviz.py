#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Base class for query result visualisations.

This is all far too ornate for what's needed in FishPath, but is a bowdlerized
version of the original & more complex code that it has been ripped from.
"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

from htmltags import *


__all__ = [
	'BaseQryViz',
]


### CONSTANTS & DEFINES ###

### IMPLEMENTATION ###	

class BaseQryViz (object):
	"""
	Base class for query result visualisations.
	
	"""
	# all these should be overriden by derived classes
	id = __module__ + '.BaseQryViz'.lower()
	label='Base query visualisation'
	description ="This is a description of BaseQryViz."
	input_type = 'This is the sort of results we expect'
	resources = []
	
	def __init__ (self, data, context=None):
		self.data = data
		self.context = context
		self.files = {}
		self.images = {}

	def render (self):
		return u'<div>blah</blah>'
	

### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	_doctest()


### END ######################################################################
