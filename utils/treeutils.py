#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
various utility functions for working with trees, usuall ete2

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

import ete2


__all__ = [
	'tree_to_phyloxml'
]


### CONSTANTS & DEFINES ###

### IMPLEMENTATION ###	
		
def tree_to_phyloxml (ete_tree):
	"""
	Convert an Ete2 tree to PhyloXML.

	:Parameters:
		ete_tree
			An Ete2 format tree

	:Returns:
		PhyloXML markup as text

	While Ete2 reads and writes Newick and PhyloXML, it cannot convert between
	the two. Hence the need for this function.

	Note that some prettyprinting takes place and that this tree is actually
	just an XML fragment and not strictly legal. It's expected it will be cleaned
	by the client function (e.g. jsPhyloSVG).

	"""
	# XXX: should have used etree and XML? 

	from cStringIO import StringIO
	buffer = StringIO()

	def visit_node (node, buf, indent=0):
		buf.write ("   " * indent)
		buf.write ("<phy:clade>\n")
		buf.write ("   " * (indent+1))
		buf.write ("<phy:name>%s</phy:name>\n" % node.name)
		buf.write ("   " * (indent+1))
		buf.write ("<phy:branch_length>%s</phy:branch_length>\n" % node.dist)
		buf.write ("   " * (indent+1))
		buf.write ("<phy:confidence type='branch_support'>%s</phy:confidence>\n" % node.support)

		for c in node.get_children():
			visit_node (c, buf, indent=indent+1)
	
		buf.write ("   " * indent)
		buf.write ("</phy:clade>\n")
		
	buffer.write ("<phy:Phyloxml xmlns:phy='http://www.phyloxml.org/1.10/phyloxml.xsd'>\n")
        buffer.write ("<phy:phylogeny>\n")
        buffer.write ("<phy:name>test_tree</phy:name>\n")

	visit_node (ete_tree.get_tree_root(), buffer)

	buffer.write ("</phy:phylogeny>\n")
	buffer.write ("</phy:Phyloxml>\n")

	return buffer.getvalue()


### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	_doctest()


### END ######################################################################
