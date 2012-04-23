#!/usr/bin/env python
"""
A simple webservice for matching sequences to a pre-existing library of seqs.

It presents a two-step form. In the first, the user (may* enter a sequence, how
sequences are to be compared and what genomic region is to be examined. In the
second step, genes from the genomic region are selected. These are then
compared as requested (fasta or neighbour joining) and the results compared.

Input should be checked and cleaned. Users should be able to step back and
forward in the steps of the program as desired, with values being populated
appropriately (e.g. returning to the start of the form to change the sequence
or the region searched).
"""

###
# CONFIGURATION DETAILS FOR LOCAL INSTALLATION CAN BE FOUND BELOW IN CONSTANTS
###


__version__ = '0.2'
__author__ = 'Paul-Michael Agapow (pma@agapow,net)'


### IMPORTS

import re
from cgi import parse_qs, escape
from exceptions import AssertionError

import config
import dbconn
import views.config as view_config
from analysis import match_seqs
from views import templates
from views import formbuilder


### CONSTANTS & DEFINES

SPACE_RE = re.compile (r'\s+')

# all comparison methods available
ALL_METHODS = [
	('fasta', 'FASTA'),
	('nj', 'Neighbour-joining')
]


### IMPLEMENTATION ###

def dev_show_traceback():
	if config.DEV_MODE:
		import sys, traceback
		print "Exception in user code:"
		print '-' * 60
		traceback.print_exc (file=sys.stdout)
		print '-' * 60
	
	
def make_db_connection():
	return dbconn.DbConn (
		protocol = config.DB_TYPE,
		host = config.DB_HOST,
		user = config.DB_USER, 
		passwd = config.DB_PASSWD,
		name = config.DB_NAME,
	)


### RENDER PAGES

def choose_regions_page (args):
	all_regions = make_db_connection().select_regions()
	rgn_tuples = [(r['id'], r['gene_region']) for r in all_regions]
	return formbuilder.select_region_form (args, ALL_METHODS, rgn_tuples)


def select_genes_page (args):
	sel_regions = make_db_connection().select_regions_by_ids (args['regions'])
	region_names = [r['gene_region'] for r in sel_regions]
	available_genes = make_db_connection().select_seqs_by_regions (region_names)
	return formbuilder.select_genes_form (args, available_genes)


def show_results_page (args):
	"""
	Perform the actions required for the search results page.
	
	That means do the match and render the results. 
	"""
	messages = []
	
	# gather the necessary arguments for everything
	sel_genes = make_db_connection().select_seqs_by_ids (args['refseqs'])
	messages.append (('note', '%s reference genes selected for matching' %
		len (sel_genes)))
	sel_regions = make_db_connection().select_regions_by_ids (args['regions'])
	region_names = [r['gene_region'] for r in sel_regions]
	messages.append (('note', 'reference genes selected from genomic regions %s' %
		', '.join(region_names)))
	target_seq = args.get ('seq', '')
	method = args['match_by']
	
	results = ''
	try:
		results = match_seqs (method, sel_genes, region_names, target_seq,
			config.EXEPATHS)
	except StandardError, err:
		dev_show_traceback()
		messages.append (('error', 'matching failed (%s)' % str(err)))
	except:
		messages.append (('error', 'matching failed (unknown problem)'))
		
	return messages, formbuilder.show_results_form (args, results)


### APPLICATION & PROCESSING

def get_args (environ):
	"""
	Unpack arguments in appropriate way.
	
	Note that the awkward code below is actually legit if you are using POST
	(which we are due to the size of what we're passing). If you use GET, it's
	a lot simpler.
	"""
	# In this idiom you must issue a list containing a default value.
	#age = d.get('age', [''])[0] # Returns the first age value.
	#hobbies = d.get('hobbies', []) # Returns a list of hobbies.
	
	# Always escape user input to avoid script injection
	#age = escape(age)
	#hobbies = [escape(hobby) for hobby in hobbies]
	
	# grab args as a dict containing lists as values
	try:
		request_body_size = int(environ.get('CONTENT_LENGTH', 0))
	except (ValueError):
		request_body_size = 0
	request_body = environ['wsgi.input'].read(request_body_size)
	args = {}
	for k, v in parse_qs(request_body).iteritems():
		if (k not in ['regions', 'refseqs']):
			if v in [[''], []]:
				v = None
			else:
				v = v[0]
		args[k] = v
		
	return args


def application(environ, start_response):

	# get the passed values, if any
	args = get_args (environ)
	
	# process request
	messages = []
	results = []
	
	# the form can be in 3 different states or stages as we progress thru form:
	# 1. opening page, enter seq, select region
	# 2. second page, select method & enter sequence
	# 3. third & final page, show results
	if (args.get ("submit", False) in [False, config.SUBMIT_SELECT_REGIONS]):
		# 1. if we are new to the form or have returned to the initial page
		form_body = choose_regions_page (args)
		
	elif args.get ("submit", False) == config.SUBMIT_SELECT_GENES:
		# 2. have just selected regions, if valid allow entry & selection of genes
		# otherwise return to region selection
		try:
			# validate number of regions
			assert (0 < len (args.get("regions", []))), \
				"need to select at least 1 region"
			# validate & clean seq
			target_seq = SPACE_RE.sub ('', args.get("seq", '')).upper()
			for i,x in enumerate (target_seq):
				assert x in 'GATCRYWSMKHBVDN-', \
					"nucleotide '%s' at position %s is not a legal symbol" % (x, i+1)
			args['seq'] = target_seq
			# validate method
			if args['match_by'] == 'fasta':
				assert target_seq, "Fasta comparison requires a target sequence"
			# okay, we're good, do the form
			messages.append (('note', 'gene regions and matching method selected'))
			form_body = select_genes_page (args)
		except StandardError, err:
			# TODO: handle err & append error message
			messages.append (('error', str(err)))
			form_body = choose_regions_page (args)

	elif args.get ("submit", False) == config.SUBMIT_MATCH_GENES:
		# 3. have entered & selected genes, if valid do match, otherwise return
		# gene entry page
		# TODO: check args
		try:
			# validate stuff
			assert (config.MIN_REFSEQS <= len (args.get("refseqs", []))), \
				"need to select at least %s reference sequences" % config.MIN_REFSEQS
			# okay, we're good, do the form
			messages.append (('note', 'genes selected'))
			msgs, form_body = show_results_page (args)
			messages.extend (msgs)
		except StandardError, err:
			dev_show_traceback()
			# TODO: handle err & append error message
			messages.append (('error', str(err)))
			form_body = select_genes_page (args)
			
	else:
		print "SHOULD NOT GET HERE", args
	
	# build page to return
	response_body = (templates.PAGE % {
		'MESSAGES': '\n'.join (['<p class="%s">%s</p>' % (r[0], r[1]) for r in messages]),
		'RESULTS': '\n'.join (results),
		'SCRIPT_NAME': __file__,
		'FORM': form_body,
	}).encode ('utf-8')
	response_headers = [
		('Content-Type', 'text/html; charset=UTF-8'),
		('Content-Length', str (len (response_body))),
	]
	start_response ('200 OK', response_headers)
	
	## Postconditions & return:
	# NOTE: have to do this as wsgi doesn't handle unicode but bytes
	# NOTE: must return an iterable
	return [response_body]


### TEST & DEBUG

# Run this file to make a server & run the application to test it
# for Windows kill in Task Manager, in Linux use Ctrl-C 
# Note that this doesn't serve static content.
if __name__ == '__main__':
	if config.DEV_MODE:
		# These settings are to be used solely for the purposes of testing and
		# development and should not be taken as advisory for production.
		
		TEST_ADDR = "158.119.147.40"
		TEST_PORT = 9123
		
		from wsgiref.simple_server import make_server
		httpd = make_server (TEST_ADDR, TEST_PORT, application)
		httpd.serve_forever()
	else:
		import wsgiref.handlers
		wsgiref.handlers.CGIHandler().run (application)


### END ###

