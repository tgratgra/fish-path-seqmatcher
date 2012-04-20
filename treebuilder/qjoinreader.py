#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Read trees from Qjoin output.

Qjoin output looks something like this::

	"Inferred tree with branch lengths:\n(('RSA71-11':0.0110195,
	'RSA71-08':0.000518949):0.0604745, (('SWZ00-AB':0.00604864,
	'SAR00-AF':0.00199426):0.0537074, 'RSA75-05':0.0666492):0.0074563,
	((((('TEST':0.577465, ('UGA99-03':0.115258,
	'UGA70-21':0.0712217):0.0517766):0.0232916, ('ZIM03-02':0.072245,
	'ZAM00-04':0.0601751):0.052884):0.00460002, ('ZAM88-02':0.0920128,
	((('ZAM93-04':0.0157011, 'ZAW96-10':0.0427784):0.0214736,
	('TAN96-01':0.0403111, 'TAN71155':0.050598):0.00579538):0.0328792,
	'TAN99-43':0.0555978):0.0113797):0.0315955):0.0340985, ('SAR81-09':0.00234728,
	'RSA81-09':0.00882591):0.0836709):0.00645376, (('SAR00-AD':0.00531535,
	'SAR00-AB':0.00518333):9.58443e-05,
	'SAR00-AC':0.0051535):0.0670438):0.016692);\n\nInferred tree with bootstrap
	values:\n(('RSA71-11':0.0110195, 'RSA71-08':0.000518949):0.0604745,
	(('SWZ00-AB':0.00604864, 'SAR00-AF':0.00199426):0.0537074,
	'RSA75-05':0.0666492):0.0074563, ((((('TEST':0.577465, ('UGA99-03':0.115258,
	'UGA70-21':0.0712217):0.0517766):0.0232916, ('ZIM03-02':0.072245,
	'ZAM00-04':0.0601751):0.052884):0.00460002, ('ZAM88-02':0.0920128,
	((('ZAM93-04':0.0157011, 'ZAW96-10':0.0427784):0.0214736,
	('TAN96-01':0.0403111, 'TAN71155':0.050598):0.00579538):0.0328792,
	'TAN99-43':0.0555978):0.0113797):0.0315955):0.0340985, ('SAR81-09':0.00234728,
	'RSA81-09':0.00882591):0.0836709):0.00645376, (('SAR00-AD':0.00531535,
	'SAR00-AB':0.00518333):9.58443e-05,
	'SAR00-AC':0.0051535):0.0670438):0.016692);\n"

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

from phylotree.io.newick import newickreader
from relais.dev.io.readers import singlereader
from relais.dev.common import *
from relais.dev.io.utils import readable_from_string


## CONSTANTS & DEFINES ###

### IMPLEMENTATION ###

class QjoinReader (singlereader.SingleReader):
	"""
	Read trees from Qjoin output.
	
	Unlike most readers, this doesn't read a physical file but the
	semi-structured output from a Qjoin run. There are a few peculiarities.
	First, qjoin always outputs two trees, being one with distances and one
	with support values in the place of distances. If no bootstraps are
	requested, this second tree still appears, but with the distances in the
	place of the support values.
	
	"""
	def __init__ (self, src):
		singlereader.SingleReader.__init__ (self, src, mode='rb')
		
	def read (self):
		# grab the actual text
		buffer = singlereader.SingleReader.read (self)
		# split out the text for the distance & support trees
		dist_buff, supp_buff = buffer.split('\n\n')
		dist_buff = readable_from_string (dist_buff.split (':', 1)[1].strip())
		supp_buff = readable_from_string (supp_buff.split (':', 1)[1].strip())
		#MSG ("dist_buff", dist_buff.getvalue())
		#MSG ("supp_buff", supp_buff.getvalue())
		# parse the trees
		rdr = newickreader.NewickReader()
		dist_tree = rdr.read (dist_buff)
		supp_tree = rdr.read (supp_buff)
		# locate the 'same' node on both trees
		#MSG (dist_tree.root)
		for n in dist_tree.iter_nodes():
			if (n.get ('title')):
				tmp_root_dist = n
				root_title = tmp_root_dist['title']
				break
		for n in supp_tree.iter_nodes():
			if (n.get ('title') == root_title):
				tmp_root_supp = n
				break
		# pair 'same' branches and copy support scores to distance tree
		branches = zip (
			[n for n in dist_tree.iter_branches_postorder (tmp_root_dist)],
			[n for n in supp_tree.iter_branches_postorder (tmp_root_supp)],
		)
		for br_pair in branches:
			dist_br, supp_br = br_pair
			dist_br['support'] = supp_br.get ('distance')
		## Postconditions & return:
		return dist_tree
		


### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod()


### MAIN ###

if __name__ == '__main__':
	_doctest()


### END ######################################################################
