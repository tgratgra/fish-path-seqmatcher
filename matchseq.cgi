#!/usr/bin/python2.6
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


__version__ = '0.4'
__author__ = 'Paul-Michael Agapow (pma@agapow,net)'


### IMPORTS

import re
from cgi import parse_qs, escape
from exceptions import AssertionError
import sys
import types

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

import config
from utils import dbconn
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
	if config.DEV_MODE or config.CAPTURE_TRACEBACKS:
		import sys, traceback
		print >> sys.stderr, "Exception in user code:"
		print >> sys.stderr, '-' * 60
		traceback.print_exc (file=sys.stderr)
		print >> sys.stderr, '-' * 60


def debug_msg(s):
	if config.DEV_MODE:
		print >> sys.stderr, "DEBUG: %s" % s
	
	
def make_db_connection():
	return dbconn.DbConn (
		protocol = config.DB_TYPE,
		host = config.DB_HOST,
		user = config.DB_USER, 
		passwd = config.DB_PASSWD,
		name = config.DB_NAME,
	)


### RENDER PAGES

def choose_pathogen_page (args):
	all_pathogens = make_db_connection().select_pathogens()
	patho_tuples = [(r['id'], '%s: %s' % (r['title'], r.get('desc', '') or '')) for r in all_pathogens]
	return formbuilder.select_pathogen_form (args, ALL_METHODS, patho_tuples)


def choose_regions_page (args):
	all_regions = make_db_connection().select_regions (args['pathogen'])
	rgn_tuples = [(r['id'], r['gene_region']) for r in all_regions]
	return formbuilder.select_region_form (args, ALL_METHODS, rgn_tuples)


def select_genes_page (args):
	sel_regions = make_db_connection().select_regions_by_ids (args['pathogen'], args['regions'])
	region_names = [r['gene_region'] for r in sel_regions]
	available_genes = make_db_connection().select_seqs_by_regions (args['pathogen'], region_names)
	return formbuilder.select_genes_form (args, available_genes)


def show_results_page (args):
	"""
	Perform the actions required for the search results page.
	
	That means do the match and render the results. 
	"""
	messages = []
	
	# gather the necessary arguments for everything
	sel_genes = make_db_connection().select_seqs_by_ids (args['pathogen'], args['refseqs'])
	messages.append (('note', '%s reference genes selected for matching' %
		len (sel_genes)))
	sel_regions = make_db_connection().select_regions_by_ids (args['pathogen'], args['regions'])
	region_names = [r['gene_region'] for r in sel_regions]
	messages.append (('note', 'reference genes selected from genomic regions %s' %
		', '.join(region_names)))
	target_seq = args.get ('seq', '')
	method = args['match_by']
	outgroup = args['outgroup']
	
	results = ''
	try:
		# will return either a table of fasta matches or an etet2 tree, rerooted 
		results = match_seqs (method, sel_genes, region_names, target_seq,
			outgroup, config.EXEPATHS)
	except StandardError, err:
		dev_show_traceback()
		messages.append (('error', 'matching failed (%s)' % str(err)))
	except:
		messages.append (('error', 'matching failed (unknown problem)'))
		
	return messages, formbuilder.show_results_form (args, results)


### APPLICATION & PROCESSING

### HANDLING PASSED ARGUMENTS ###

class CgiArgs (object):

	def __init__ (self, env):
		self._args = get_args (env)

	def get (self, key, default=None):
		return self._args.get (key, default)

	def get_list (self, key):
		return self.get (key, [])

	def __getitem__ (self, key):
		return self.get (key)

	def set (self, key, val):
		self._args[key] = val

	def __setitem__ (self, key, val):
		self.set (key, val)

	def append_list (self, key, val):
		old_val = self.get_list (key)
		assert (type (old_val) == type ([]))
		old_val.append (val)
		self.set (key, old_val)

	def set_list (self, key, val_list):
		assert (type (val_list) in [types.ListType, types.Tuple.Type])
		self.set (List (val_list))


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
	if config.REQUEST_METHOD == 'POST':	
		# grab args as a dict containing lists as values
		try:
			request_body_size = int(environ.get('CONTENT_LENGTH', 0))
		except (ValueError):
			request_body_size = 0
		request_body = environ['wsgi.input'].read(request_body_size)
		raw_arg_dict = parse_qs (request_body)
	else:
		# method is 'GET'
		debug_msg("querystring=%s" % environ.get ('QUERY_STRING', ''))
		raw_arg_dict = parse_qs(environ.get ('QUERY_STRING', ''))
		debug_msg(raw_arg_dict)

	args = {}
	for k, v in raw_arg_dict.iteritems():
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
	
	# the form can be in 4 different states or stages as we progress thru form:
	# 0. opening page. select the pathogen
	# 1. enter seq, select region
	# 2. second page, select method & enter sequence
	# 3. third & final page, show results
	if (args.get ("submit", False) in [False, config.SUBMIT_SELECT_PATHOGEN]):
		# 0. if we are new to the form or have returned to the initial page
		debug_msg("submit=selectpathogen")
		form_body = choose_pathogen_page (args)
		
	if (args.get ("submit", False) == config.SUBMIT_SELECT_REGIONS):
		# 1. have selected pathogen
		# could popssibly route straight to here
		# must validate pathogen choice
		debug_msg("submit=selectregions")
		form_body = choose_regions_page (args)
		
	elif args.get ("submit", False) == config.SUBMIT_SELECT_GENES:
		# 2. have just selected regions, if valid allow entry & selection of genes
		# otherwise return to region selection
		debug_msg("submit=selectgenes")
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
			dev_show_traceback()
			messages.append (('error', str(err)))
			form_body = choose_regions_page (args)

	elif args.get ("submit", False) == config.SUBMIT_MATCH_GENES:
		# 3. have entered & selected genes, if valid do match, otherwise return
		# gene entry page
		# TODO: check args
		debug_msg("submit=matchgenes")
		try:
			## extract, validate, & cleanse arguments ...
			# ... outgroup may or may not be in refseqs, add if necessary
			ref_seqs = args.get("refseqs", [])
			outgroup = args.get ('outgroup', None)
			if outgroup and (outgroup not in ref_seqs):
				ref_seqs.append (outgroup)
				args['ref_seqs'] = ref_seqs
			# ... check min number of ref seqs
			assert (config.MIN_REFSEQS <= len (ref_seqs)), \
				"need to select at least %s reference sequences" % config.MIN_REFSEQS

			## okay, we're good, now process the page
			messages.append (('note', 'genes selected'))
			msgs, form_body = show_results_page (args)
			messages.extend (msgs)

		except StandardError, err:
			dev_show_traceback()
			# TODO: handle err & append error message
			messages.append (('error', str(err)))
			form_body = select_genes_page (args)
			
	else:
		sys.stderr.write ("%s\n" % args)
		print >> sys.stderr, "SHOULD NOT GET HERE (css fetch error in debug may cause)", args
	
	# build page to return
	response_body = (templates.PAGE % {
		'MESSAGES': '\n'.join (['<p class="%s">%s: %s</p>' % (r[0], r[0].title(), r[1]) for r in messages]),
		'RESULTS': '\n'.join (results),
		'SCRIPT_NAME': __file__,
		'METHOD': config.REQUEST_METHOD,
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
		from wsgiref.simple_server import make_server
		httpd = make_server (config.TEST_ADDR, config.TEST_PORT, application)
		httpd.serve_forever()
	else:
		import wsgiref.handlers
		wsgiref.handlers.CGIHandler().run (application)


### END ###

