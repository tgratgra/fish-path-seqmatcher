#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Layout for and generation of the webapp forms.

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

import types
import exceptions

import templates
import htmltags
import config

from fastamatchviz import *
from phylotextviz import *
from phylosvgviz import *


### CONSTANTS & DEFINES ###

### IMPLEMENTATION ###	

def form_field (label, inputs, helptext=''):
	""""
	Consistent formatting for for fields.
	
	:Parameters:
		label
			name of the field
		content
			HTML for the field controls or inputs
		helptext
			descriptive text for using the field
			
	:Returns:
		HTML for the field, to be included in a form
		
	"""
	return templates.FORM_FIELD % {
		'LABEL': label,
		'INPUTS': inputs,
		'HELPTEXT': helptext,
	}


def radio_input (label, name, vals, helptext='', env={}, default=None):
	"""
	:Parameters:
		vals
			a list of value - label pairs
		env
			the environ or request dict
	"""
	## Preconditions & preparation:
	checked_vals = env.get (name, [vals[0][0]])

	## Main:
	radio_vals = []
	for val, choice_label in vals:
		if val in checked_vals:
			checked = 'checked'
		else:
			checked = ''
		radio_vals.append ([val, choice_label, checked])
	
	radio_body = htmltags.tag_with_contents ('div',
		_class = 'radio',
		contents = "<br />\n".join ([
			"<input type='radio' name='%s' value='%s' %s />%s" % (name, val, check, choice_label) for 
			val, choice_label, check in radio_vals
		])
	)

	## Postconditions & return:
	return form_field (label, radio_body, helptext)


def checkbox_input (label, name, vals, helptext='', env={}):
	"""
	:Parameters:
		name
			the value being captured
		vals
			a list of value - label pairs
		env
			the environ or request dict
	"""
	checked_vals = env.get (name, [])
	radio_vals = []
	for val, n in vals:
		if val in checked_vals:
			checked = 'checked'
		else:
			checked = ''
		radio_vals.append ([val, n, checked])
	
	radio_tmpl = "<input type='checkbox' name='%s' value='%s' %s />%s"
	
	check_body = "<div class='checkbox'>%s</div>" % (
		"<br />\n".join ([radio_tmpl % (name, val, check, choice_label) for 
			val, choice_label, check in radio_vals])
	)
	
	return form_field (label, check_body, helptext)


def hidden_input (name, val):
	"""
	Produce the HTML for a hidden field
	
	:Parameters:
		name
			name of the field or value
		val
			value of the field, can be a list or tuple of multiple values
			
	:Returns:
		HTML for hidden input. If multiple values are passed, multiple fields are
		returned, concatenated.
			
	Used for passing form values between stages.
	"""
	if type (val) not in (types.ListType, types.TupleType):
		val = list(val)
	return '\n'.join (["<input type='hidden' name='%s' value='%s' />" %
		(name, v) for v in val])



def seq_table_input (gene_choices, env={}):
	if gene_choices:
		
		def format_vals (v):
			if v is None:
				return ''
			else:
				return "%s" % v
			

		col_headers = [
			'accession_number',
			'isolate_name',
			'sequencing_method',
			'author'
		]
		
		theader = htmltags.tag_with_contents ('thead',
			htmltags.tag_with_contents ('tr',
				' '.join ([
					htmltags.tag_with_contents ('td', 
						h.replace('_', ' ').title()
					) for h in [''] + col_headers
				])
			)
		)
		
		checked_vals = dict (*[(v, 'checked') for v in env.get ('refseqs', [])])
		tbody = htmltags.tag_with_contents ('tbody',
			'\n'.join ([
				htmltags.tag_with_contents ('tr',
					' '.join (
						[
							"<td><input type='checkbox' name='refseqs' value='%s' %s /></td>" % \
								(g['id'], checked_vals.get(g['id'], ''))
						] + 
						[htmltags.tag_with_contents ('td', format_vals (g[h])) for h in col_headers]
					)
				) for g in gene_choices
			])
		)
		
		tab = htmltags.tag_with_contents ('table',
			theader + '\n' + tbody,
			_class="refseq_table"
		)
		helptext = """
			<p class='helptext'>Select at least 3 reference sequences
			to be matched. Select 
			<a href="javascript:SetAllCheckBoxes('dvifish', 'refseqs', true)">all</a> or
			<a href=\"javascript:SetAllCheckBoxes('dvifish', 'refseqs', false)">none</a>.
			</p>
		"""
		return tab + helptext

	

### FORMS / PAGES
# The three pages to show

def select_region_form (d, methods, region_choices):

	seq_inputs = """
		<textarea name="seq" cols="60" rows="5" wrap="physical">%s</textarea>
	""" % d.get ('seq', '')
	
	fields = '\n'.join ([
		form_field ('Sequence to be matched', seq_inputs, """Enter a
			molecular sequence of DNA, stop codon not included. If no sequence is
			entered, matching will only be done within the alignment of sequences
			selected later selected below."""),
		radio_input ('Method', 'match_by', methods, env=d,
			helptext="How to match sequences.", default=None),
		checkbox_input ('Region', 'regions', region_choices, env=d,
			helptext="""
				The genomic region to be matched against. Select 
				<a href="javascript:SetAllCheckBoxes('dvifish', 'regions', true)">all</a> or
				<a href=\"javascript:SetAllCheckBoxes('dvifish', 'regions', false)">none</a>.
			"""
		),
	])

	controls = """
		<input type="hidden" name="_form_submitted" value="True" />
		<input type="reset" value="Reset" />
		<input type="submit" name="submit" value="%s" />	
	""" % config.SUBMIT_SELECT_GENES
	
	## Return:
	return templates.FORM_BODY % {
		'TITLE': 'Step 1: Enter sequence, select method &amp; region',
		'DESC': '',
		'FIELDS': fields,
		'CONTROLS': controls,
	}
	
	
def select_genes_form (d, gene_choices):

	hidden_inputs = '\n'.join ([
		hidden_input ('seq', [d.get ('seq', '')]),
		hidden_input ('regions', d['regions']),
		hidden_input ('match_by', [d['match_by']]),
		])
	
	#fields = '\n'.join ([
	#	checkbox_input ('Select genes', 'genes', gene_choices, env=d,
	#		helptext="The genomic region to be matched against."),
	#])

	fields = seq_table_input (gene_choices)
	
	controls = """
		<input type="hidden" name="_form_submitted" value="True">
		<input type="reset" value="Reset">
		<input type="submit" name="submit" value="%s">
		<input type="submit" name="submit" value="%s">
	""" % (config.SUBMIT_SELECT_REGIONS , config.SUBMIT_MATCH_GENES)
	
	## Return:
	return templates.FORM_BODY %  {
		'TITLE': 'Step 2: Select genes &amp; compare',
		'DESC': '',
		'FIELDS': hidden_inputs + '\n\n' + fields,
		'CONTROLS': controls,
	}


def show_results_form (d, results):
	
	results_viz = ''
	
	# so we can handle multiple results
	for r in results:
		print "R", r
		if r[0] in ['tree', 'newicktree']:
			results_viz += PhyloTextViz (r[1]).render()
			results_viz += PhyloSvgViz (r[1]).render()
		elif r[0] in ['fastamatchs']:
			results_viz += FastaMatchViz (r[1]).render()
		else:
			raise ValueError, "unknown result type '%s'" % results[0]
		
	hidden_inputs = '\n'.join ([
		hidden_input ('seq', [d.get ('seq', '')]),
		hidden_input ('regions', d['regions']),
		hidden_input ('refseqs', d['refseqs']),
		hidden_input ('match_by', [d['match_by']]),
		])
	
	controls = """
		<input type="hidden" name="_form_submitted" value="True">
		<input type="submit" name="submit" value="%s">
		<input type="submit" name="submit" value="%s">
	""" % (config.SUBMIT_SELECT_REGIONS, config.SUBMIT_SELECT_GENES)
	
	## Return:
	return templates.FORM_BODY %  {
		'TITLE': 'Step 3: Show results',
		'DESC': '',
		'FIELDS': '\n\n'.join ([results_viz, hidden_inputs]),
		'CONTROLS': controls,
	}
	
	


### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	main()


### END ######################################################################

