#!/usr/bin/env python
"""
The overall layout of the webapp page.

This includes whatever appears in the header, look-and-feel, page-title, etc.
Customize this (carefully) as you so wish.

"""

### INMPORTS

### CONSTANTS & DEFINES

PAGE = """
<html>
	<head>
		<!-- TODO: just for dev so wsgi doesn't  -->
		<title>TESTING TESTING</title>
		<!-- <link rel="stylesheet" type="text/css" href="./matchseq.css" /> -->
		<link rel="shortcut icon" href="http://www.example.com/my_empty_resource" />
	<head>
	<body>
	
		<!-- TODO: alter page title if need be -->
		<h1>Sequence Matching</h1>
	
		<p class="description">
			This webservice accepts an input sequence and attempts to match it
			against selection of database sequences, using either FASTA matching
			or by building a phylogeny using Neighbour-joining.
		</p>
	
		%(MESSAGES)s
	
		%(RESULTS)s
	
		<hr />
	
	   <form method="post" action="%(SCRIPT_NAME)s">
			%(FORM)s
		</form>
	
		<hr />
		
		<!-- TODO: alter footer info if need be -->
		<div class='footer'>
			
		</div>
		
	</body>
</html>
"""

# the layout for a form
FORM_BODY = """
<H2>%(TITLE)s</H2>
<p class="description">%(DESC)s</p>

%(FIELDS)s

<HR class="form_divider" />

<div class="form_controls">
%(CONTROLS)s
</div>
""" 

# the layout for a field in a form
FORM_FIELD = """
<div class='form_field'>
	<label>%(LABEL)s</label><br />
	%(INPUTS)s
	<p class='helptext'>%(HELPTEXT)s</p>
</div>
"""



### END ###
