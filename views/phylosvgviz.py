#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Visualise a phylogeny in SVG using jsPhyloSVG

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

import re
from copy import deepcopy
from math import log10, floor

from htmltags import *
import baseqryviz
from utils import treeutils, jsutils

import config

__all__ = [
	'PhyloSvgViz', 
]


### CONSTANTS & DEFINES ###

### IMPLEMENTATION ###	
		
class PhyloSvgViz (baseqryviz.BaseQryViz):
	id = __module__ + '.PhyloSvgViz'.lower()
	label='Visualise a phylogeny  in SVG'
	resources = []
	
	def copyAndCleanTree (self):
		"""
		Make the tree suitable for use by jsPhyloSVG
		"""
		# TODO: Need to do several things here:
		# - NoNames
		# - copy support scores to internal branch names

		## Main:
		# Copy the tree so as not to damage original
		ete_tree = deepcopy (self.data)

		# set root branch to zero, make change later
		ete_tree.dist = 0.0

		# find max / min branchlength for diagnostic purposes
		# doesn't use negative or zero branch lengths
		# Also clean names
		max_bl = None
		min_bl = None
		for n in ete_tree.traverse ("postorder"):
			if (0.0 < n.dist):
				if (max_bl is None) or (max_bl < n.dist):
					max_bl = n.dist
				if (min_bl is None) or (n.dist < min_bl):
					min_bl = n.dist
			clean_name = n.name.strip()
			if (clean_name[0] == "'") and (clean_name[-1] == "'"):
				clean_name = clean_name[1:-1]
			n.name = clean_name

		# set all branches to be at least 1/100 of the largest or 1/10 the
		# smallest, whichever is larger
		default_bl = max (max_bl / 100, min_bl/10)
		for n in ete_tree.traverse ("postorder"):
			if (n.dist <= 0.0):
				n.dist = default_bl

		# get support values on tree by setting supprt as name
		for n in ete_tree.traverse ("postorder"):
			# if an internal node
			if (not n.is_leaf()):
				n.name = config.SUPPORT_FMT % n.support	

		# very hacky - calc appropriate scale bar size and stick on root
		magn = int (floor (log10 (max_bl)))
		scale_size = 10**magn
		ete_tree.scale_size = scale_size

		## Postcondtions & return:int ( floor ( log10 (x)))
		return ete_tree


	def render (self):
		ete_tree = self.copyAndCleanTree()

		tree_xml = treeutils.tree_to_phyloxml (ete_tree)
		js_xml = jsutils.make_js_str (tree_xml)
		
		debug_str = """1:
			<textarea rows="20" cols="60">%s</textarea>
			2:
			<textarea rows="20" cols="60">%s</textarea>
		""" % (tree_xml, js_xml)

		# return appropriate 
		return """
<div id="svgCanvas"> </div>
<script type="text/javascript">
	Smits.PhyloCanvas.Render.Style.bootstrap = {
        "font-family":  'Verdana',
        "font-style":  'italic',
        "font-size":    12,
        "text-anchor":  'start'
    };
	Smits.PhyloCanvas.Render.Parameters.Rectangular.showScaleBar = %(SCALE)s;
	var dataObject = { phyloxml: '%(TREE_STR)s' };
	phylocanvas = new Smits.PhyloCanvas(
		dataObject,
		'svgCanvas', 
		%(HT)s, %(WT)s,
		'%(SHAPE)s'
	);
</script>
	""" % {
		'SCALE'     : ete_tree.scale_size,
		'TREE_STR'  : js_xml,
		'SHAPE'     : config.TREE_SHAPE,
		'HT'        : config.TREE_HT,
		'WT'        : config.TREE_WT,
	}



### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	_doctest()


### END ######################################################################
