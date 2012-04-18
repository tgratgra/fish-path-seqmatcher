#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Layout for and generation of the webapp forms.

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

import templates
import htmltags
import config


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


def render_hidden (name, val):
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
	if type (val) not in (types.ListType, typesdTupleType):
		val = List(val)
	return '\n'.join (["<input type='hidden' name='%s' value='SEQ' />" %
		(name, v) for v in vals])


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
			helptext="The genomic region to be matched against."),
	])

	controls = """
		<input type="hidden" name="_form_submitted" value="True" />
		<input type="reset" value="Reset" />
		<input type="submit" name="submit" value="%s" />	
	""" % config.SUBMIT_SELECT_REGIONS
	
	## Return:
	return templates.FORM_BODY %  {
		'TITLE': 'Step 1: Enter sequence, select method &amp; region',
		'DESC': '',
		'FIELDS': fields,
		'CONTROLS': controls,
	}
	
	
def enter_and_select_genes_form (d, gene_choices):
	controls = """
		<input type="hidden" name="_form_submitted" value="True">
		<input type="reset" value="Reset">
		<input type="submit" name="select_genes" value="Select genes">	
	"""
	hidden_inputs = '\n'.join ([render_hidden (n, d.get(n, '')) for n in
		['seq', 'region', 'match_by']])
	
	fields = '\n'.join ([
		render_check_input ('Select genes', 'genes', gene_choices, env=d,
			helptext="The genomic region to be matched against."),
	])

	controls = """
		<input type="hidden" name="_form_submitted" value="True">
		<input type="reset" value="Reset">
		<input type="submit" name="submit" value="%s">
		<input type="submit" name="submit" value="%s">
	""" % (config.SUBMIT_RESELECT_REGIONS , config.SUBMIT_MATCH_GENES)
	
	## Return:
	return templates.FORM_BODY %  {
		'TITLE': 'Step 2: Select genes &amp; compare',
		'DESC': '',
		'FIELDS': hidden_inputs + fields,
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

