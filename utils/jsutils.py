#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Various utility functions for working with javascript.

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

import re


__all__ = [
	'make_js_str',
]


### CONSTANTS & DEFINES ###

STR_ESC_RE = re.compile (r'([\n\'\"])')

### IMPLEMENTATION ###	
		
def make_js_str (str):
	r"""
	Make a string safe for use in javascript.

	:Parameters:
		str
			a string, including quotes and newlines

	:Returns:
		the string with all the appropriate characters escaped

	Javascripts lame way of doing multiline strings is to have every newlines
	backslash escaped. We also have the issue of quotes inside the string. So,
	we just escape everything. 

	For example::

		>>> print make_js_str ('foo')
		foo
		>>> print make_js_str ('foo"bar')
		foo\"bar
		>>> print make_js_str ('\n'.join (['foo', 'bar', 'baz']))
		foo\
		bar\
		baz

	"""
	# NOTE: a real arse writnig doctests that include newlines: make the
	# docstring a raw string, remember that \n in a rawstring is actually
	# \n (escaped n), not a newline
	return STR_ESC_RE.sub (r'\\\1', str)



### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	_doctest()


### END ######################################################################
