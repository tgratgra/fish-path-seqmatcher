#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Renders Fasta results as a table.

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

from htmltags import *
import baseqryviz


__all__ = [
	'FastaMatchViz',
]


### CONSTANTS & DEFINES ###

### IMPLEMENTATION ###	
		
class FastaMatchViz (baseqryviz.BaseQryViz):
	id = __module__ + '.FastaMatchViz'.lower()
	label='Table of Fasta results'
	resources = []
	
	def render (self):
		if (len (self.data) == 0):
			return tag_with_contents ('p', 'FASTA found no close matches.')
			
		rows = []
		
		header_names = [
			'Seq. ID',
			'Nt. overlap',
			'Identity',
			'Nt. matched',
			# 'Nt. Ambig. or gap',
			'Z score',
			'Bits',
		]
		header_cells = [tag_with_contents ('th', x, class_="nosort column") for x in header_names]
		rows.append (header_cells)
		
		for record in self.data:
			curr_row = []
			for val in [
				record.seq_id,
				'%d' % record.res_compared,
				'%.2f' % record.frac_identity,
				'%d' % record.res_matched,
				#'%d' % record.res_ambigs,
				'%.2f' % record.zscore,
				'%.2f' % record.bits,
			]:
				curr_row.append (tag_with_contents ('td', val))
			rows.append (curr_row)
			
		row_text = [tag_with_contents ('tr', ''.join (row)) for row in rows]
		thead_text = tag_with_contents ('thead', row_text[0])
		tbody_text = tag_with_contents ('tbody', '\n'.join (row_text[1:]))
		table_text = tag_with_contents ('table', thead_text + tbody_text,
			class_="listing center", id_="fasta_res")
		
		leader_txt = tag_with_contents ('p', "The following are the best %s matches with the test sequence:" % len (self.data))
		
		return tag_with_contents ('div', leader_txt + table_text)
	


### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	_doctest()


### END ######################################################################
