#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Visualise a phylogeny as a image.

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

import StringIO

from phylotree.simpledraw import *

from relais.dev.pilutils import crop_to_content

from relais.dev.htmltag import *
import baseqryviz


__all__ = [
	'PhyloImageViz',
]


### CONSTANTS & DEFINES ###

### IMPLEMENTATION ###	
		
class PhyloImageViz (baseqryviz.BaseQryViz):
	id = __module__ + '.PhyloImageViz'.lower()
	label='Image representation of phylogeny'
	resources = []
	
	def render (self):
		tree = self.data
		
		coords = plot_radial_coords (tree)
		fit_coords (coords)
		canvas_size = int (float (tree.count_tip_nodes())**0.5) * 150 + 50
		drwr = TreeDrawer (canvas_size)
		drwr.draw_tree (tree, coords)
		
		img = drwr.save_to_pilimage ()
		bg_clr = img.getpixel ((0, 0))
		new_img = crop_to_content (img, bg_clr, 15)
		buff = StringIO.StringIO()
		new_img.save (buff, format='PNG')
		self.images['phylo_match.png'] = buff.getvalue()
		
		results_text = start_tag ('div')
		
		results_text += closed_tag ('img',
			src='phylo_match.png',
			class_='result_img align_center',
			style="max-width: 90%; margin-left: auto; margin-right: auto; display: block;",
		)
		
		results_text += tag_with_contents ('a',
			'Click to enlarge',
			href='phylo_match.png',
			target='_blank',
			style="text-align: center; display: block;",
		)
		
		results_text += stop_tag ('div')
		return results_text
	


### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	_doctest()


### END ######################################################################
