#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
various utility functions for working with trees, usuall ete2

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

# A hairy problem: ete2 has a number of depedencies that it doesn't
# necessarily depend upon: MySQLdb, Numpy, PyQt, etc. If you don't use
# the associated functionality, you won't need these dependencies. But
# ete2 tries to import them anyway and write warning messages to stdout
# and stderr for each one, before importing ete2 successfully. Stray output
# could create problems in a cgi envionment, so these have to be
# suppressed. Hence the convuluted import below.

import os, sys

class SuppressAllOutput (object):
	def __enter__(self):
		sys.stderr.flush()
		self.old_stderr = sys.stderr
		sys.stderr = open('/dev/null', 'a+', 0)
		sys.stdout.flush()
		self.old_stdout = sys.stdout
		sys.stdout = open('/dev/null', 'a+', 0)

	def __exit__(self, exc_type, exc_value, traceback):
		sys.stderr.flush()
		sys.stderr = self.old_stderr
		sys.stdout.flush()
		sys.stdout = self.old_stdout

with SuppressAllOutput():
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
