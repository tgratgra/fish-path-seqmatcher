#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Various application-wide configuration details.

This should be the only place that settings (primarily db settings) have to be
adjusted.
"""


### CONSTANTS & DEFINES

DEV_MODE = True

### Db connection details

if DEV_MODE:
	DB_TYPE = 'sqlite'
	DB_HOST = "fish.sqlite"
	DB_USER = None
	DB_PASSWD = None
	DB_NAME = None
else:
	DB_TYPE = 'mysql'
	DB_HOST = "localhost"
	DB_USER = "fishpathogens"
	DB_PASSWD = "PgDenv23"
	DB_NAME = "fishpathogens"


### submission buttons for forms
SUBMIT_SELECT_REGIONS = 'Select regions'
SUBMIT_RESELECT_REGIONS = 'Re-select regions'
SUBMIT_MATCH_GENES = 'Match genes'


### maximum sequences to return
SEQLIMIT = 150


### END ######################################################################
