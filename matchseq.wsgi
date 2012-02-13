#!/usr/bin/env python

### IMPORTS

from wsgiref.simple_server import make_server
from cgi import parse_qs, escape


### CONSTANTS & DEFINES

TEST_ADDR = "158.119.147.40"
TEST_PORT = 9123

CSS_URL = '/matchseq.css'
SITE_TITLE = "FishPathogens.Eu"
PAGE_TITLE = "Sequence matcher"

TEST_REGIONS = [
	('1', 'Region 1'),
	('2', 'Region 2'),
	('3', 'Region 3'),
]

SCRIPT_NAME = __file__

ALL_METHODS = [
	('fasta', 'FASTA'),
	('nj', 'neighbour-joining')
]

PAGE_TMPL = """
<html>
	<head>
		<title>%(SITE_TITLE)s</title>
		<link rel="stylesheet" type="text/css" href="%(CSS_URL)s" />
	<head>
	<body>
	
		<h1>Sequence Matching</h1>
	
		<p class="description">
			This webservice accepts an input sequence and attempts to match it against
			selection of database sequences, using either FASTA matching or by
			building a phylogeny using Neighbour-joining.
		</p>
	
		%(MESSAGES)s
	
		%(RESULTS)s
	
		<hr />
	
	   <form method="get" action="%(SCRIPT_NAME)s">
			%(FORM)s
		</form>
	
		<hr />
		THE FOOTER WILL APPEAR HERE
		
	</body>
</html>
""" 
	



### IMPLEMENTATION ###

def render_radio_input (name, vals, env):
	"""
	:Parameters:
		vals
			a list of value - label pairs
		env
			the environ or request dict
	"""
	checked_vals = env.get (name, [])
	radio_vals = []
	for val, label in vals:
		if val in checked_vals:
			checked = 'checked'
		else:
			checked = ''
		radio_vals.append ([val, label, checked])
	
	radio_tmpl = "<input type='radio' name='%s' value='%s' %s />%s"
	
	return "<div class='radio'>%s</div>" % (
		"<br />".join ([radio_tmpl % (name, val, check, label) for 
			val, label, check in radio_vals])
	)
	

def render_checkbox_input (name, vals, env):
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
	for val, label in vals:
		if val in checked_vals:
			checked = 'checked'
		else:
			checked = ''
		radio_vals.append ([val, label, checked])
	
	radio_tmpl = "<input type='checkbox' name='%s' value='%s' %s />%s"
	
	return "<div class='checkbox'>%s</div>" % (
		"<br />".join ([radio_tmpl % (name, val, check, label) for 
			val, label, check in radio_vals])
	)
	
	
def render_form_select_region (d):
	return """
		<H2>Step 1: Enter sequence &amp; select region</H2>
		<p class="description">%(DESC)s</p>
		
		<div class="form_field">
			<label>Match sequence</label>
			<textarea name="seq" cols="60" rows="5" wrap="physical">%(SEQ)s</textarea>
			<p class="helptext">
				Enter a molecular sequence of DNA, stop codon not included.
				If no sequence is entered, matching will only be done within
				the alignment will be sequences selected below
			</p>
		</div>
		
		<div class="form_field">
			<label>Method</label>
				%(METHOD_RADIOS)s
			<p class="helptext">
				The method to use to match sequences by.
			</p>
		</div>
		
		<div class="form_field">
			<label>Select region</label>
			<input name="region" type="radio" value="ABC">ABC<br />
			<input name="region" type="radio" value="DEF">ABC<br />
			<input name="region" type="radio" value="ABC">ABC<br />
			<p class="helptext">
				Choose a genomic region for the sequence to be matched against.
			</p>
		</div>
		
		<div class="form_controls">
			<input type="hidden" name="_form_submitted" value="True">
			<input type="reset" value="Reset">
			<input type="submit" name="select_genes" value="Select genes ->">
		</div>
	""" %  {
		'DESC': '',
		'SEQ': d.get (seq, ''),
		'METHOD_RADIOS': render_radio_input ('match_by', ALL_METHODS, d),
		'REGION_CHECKS': render_check_input ('match_by', ALL_METHODS, d),
	}
	
	
def render_form_select_genes (d):
	return """
		<H2>Step 1: Select genes</H2>
		<p class="description"></p>
		
		<input type="hidden" name="seq" value="SEQ" />
		<input type="hidden" name="region" value="RGN" />
		<input type="hidden" name="match_by" value="MTHD" />
		
		<div class="form_field">
			<label>Select genes</label>
			<input name="region" type="checkbox" value="ABC">ABC<br />
			<input name="region" type="checkbox" value="DEF">ABC<br />
			<input name="region" type="checkbox" value="ABC">ABC<br />
				Enter a molecular sequence of DNA. Use only 'ACGT'.
			</p>
		</div>
		
		<div class="form_controls">
			<input type="hidden" name="_form_submitted" value="True">
			<input type="reset" value="Reset">
			<input type="submit" name="select_region" value="<- Select gene region">
			<input type="submit" name="find_match" value="Find match ->">
		</div>
	"""


def application(environ, start_response):

	# grab args as a dict containing lists as values
	args = parse_qs(environ['QUERY_STRING'])
	print args
	
	# In this idiom you must issue a list containing a default value.
	#age = d.get('age', [''])[0] # Returns the first age value.
	#hobbies = d.get('hobbies', []) # Returns a list of hobbies.
	
	# Always escape user input to avoid script injection
	#age = escape(age)
	#hobbies = [escape(hobby) for hobby in hobbies]
	
	status = '200 OK'

	messages = []
	results = []
	
	if ((args.get ("_form_submitted", False) == False) or 
			(args.get ("select_region", False))):
		# if we are new to the form
		form_body = render_form_select_region (args)
	elif args.get ("select_genes", False):
		# have selected region, now want to select genes
		# TODO: check args
		if (True):
			form_body = render_form_select_genes (args)
		else:
			# go back to the page
			# TODO: append error message
			form_body = render_form_select_region (args)
	elif args.get ("find_match", False):
		# have selected genes, now want to do the match
		# TODO: check args
		if (True):
			# TODO: render results
			form_body = render_form_select_genes (args)
		else:
			# go back to the page
			# TODO: append error message
			form_body = render_form_select_genes (args)
	else:
		# if we are new to the form
		form_body = render_form_select_region (args)
	
	response_body = PAGE_TMPL % {
		'SITE_TITLE': 'FishPathogens.eu',
		'CSS_URL': CSS_URL,
		'MESSAGES': '\n'.join (['<p class="%s">%s</p>' % (r[0], r[1]) for r in messages]),
		'RESULTS': '',
		'SCRIPT_NAME': __file__,
		'FORM': form_body,
	}
	response_headers = [
		('Content-Type', 'text/html'),
		('Content-Length', str (len (response_body))),
	]
	start_response(status, response_headers)
	
	return [response_body]


### TEST & DEBUG

# make a server & run the application to test it
if __name__ == '__main__':
	httpd = make_server (TEST_ADDR, TEST_PORT, application)
	# for Windows kill in Task Manager, in Linux use Ctrl-C 
	httpd.serve_forever()


### END ###


