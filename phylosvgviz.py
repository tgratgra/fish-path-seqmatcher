#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Visualise a phylogeny  in SVG using jsPhyloSVG

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

import re

from htmltags import *
import baseqryviz


__all__ = [
	'PhyloSvgViz',
]


### CONSTANTS & DEFINES ###

### IMPLEMENTATION ###	
		
class PhyloSvgViz (baseqryviz.BaseQryViz):
	id = __module__ + '.PhyloSvgViz'.lower()
	label='Visualise a phylogeny  in SVG'
	resources = []
	
	def render (self):
		tree_str = self.data
		
		# need to clean up data
		tree_str = re.sub ("'", '', tree_str)
		tree_str = re.sub (r'\s+', '', tree_str)
		
		# return appropriate 
		return """
<div id="svgCanvas"> </div>
<script type="text/javascript">
	var dataObject = { newick: '%s' };
	phylocanvas = new Smits.PhyloCanvas(
		dataObject,
		'svgCanvas', 
		500, 500
	);
</script>
	""" % tree_str



### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	_doctest()


### END ######################################################################
