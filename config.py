#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Various application-wide configuration details.

This should be the only place that settings (primarily db settings) have to be
adjusted.
"""


### CONSTANTS & DEFINES

DEV_MODE = False
CAPTURE_TRACEBACKS = True

# where is this running?
LOCALE = 'AN'
#LOCALE = 'DTU'

# IP for dev, DTU is live so doesn't need this
if LOCALE == 'AN':
	TEST_ADDR = "72.14.179.53"

TEST_PORT = 9123


# Db connection details

if LOCALE in [ 'AN', 'HPA']:
	DB_TYPE = 'sqlite'
	DB_HOST = "fish.sqlite"
	DB_USER = None
	DB_PASSWD = None
	DB_NAME = None
else:
	DB_TYPE = 'mysql'
	DB_HOST = "localhost"
	DB_USER = ""
	DB_PASSWD = ""
	DB_NAME = ""


# submission buttons for forms
#from views.config import *


# maximum sequences to return from region selection
SEQLIMIT = 150

# user must select at least this many reference sequences
MIN_REFSEQS = 3

# for phyogeny reconstruction
NUM_BOOTSTRAPS = 100


# how to write branchlengths on the trees
SUPPORT_FMT = "%1.2f"

SUBMIT_SELECT_PATHOGEN = '0: Select pathogen'
SUBMIT_SELECT_REGIONS = '1: Select regions'
SUBMIT_SELECT_GENES = '2: Select sequences'
SUBMIT_MATCH_GENES = '3: Match sequences'

TREE_SHAPE = 'square'
TREE_HT = 1000
TREE_WT = 1000


# executable paths
if LOCALE == 'HPA':
	EXEPATHS = {
		'mafft': '/usr/bin/mafft',
		'fasta': '/home/f0/paul/Installed/bin/fasta36',
		'qjoin': '/home/f0/paul/Bin/qjoin',
		'fasttree': '/usr/local/bin/FastTree',
	}
elif LOCALE == 'AN':
        EXEPATHS = {
                'mafft': '/usr/local/bin/mafft',
                'fasta': '/usr/local/bin/fasta36',
                'qjoin': '/usr/local/bin/qjoin',
		'fasttree': '/usr/local/bin/FastTree',
        }
elif LOCALE == 'DTU':
	EXEPATHS = {
                'mafft': '/home/tgray/apps/mafft/bin/mafft',
                'fasta': '/home/tgray/apps/fasta-35.4.11/bin/fasta35',
                'qjoin': '/home/tgray/apps/quick-join-1.0.10/qjoin',
		'fasttree': '/usr/local/bin/FastTree',
        }



	

REQUEST_METHOD = 'GET'
#REQUEST_METHOD = 'POST'


### END ######################################################################
