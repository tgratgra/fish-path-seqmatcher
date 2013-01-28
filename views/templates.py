#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
		<title>FishPathogens.Eu Sequence Matcher</title>

		<link rel="shortcut icon" href="http://localhost/" />

		<link rel="stylesheet" type="text/css" href="./media/css/matchseq.css" />
		<script type="text/javascript" src="./media/js/raphael-min.js" ></script> 
		<script type="text/javascript" src="./media/js/jsphylosvg-min.js"></script>
		<script language='javascript'>
			function SetAllCheckBoxes(FormName, FieldName, CheckValue) {
				if(!document.forms[FormName])
					return;
				var objCheckBoxes = document.forms[FormName].elements[FieldName];
				if(!objCheckBoxes)
					return;
				var countCheckBoxes = objCheckBoxes.length;
				if(!countCheckBoxes)
					objCheckBoxes.checked = CheckValue;
				else
					// set the check value for all check boxes
					for(var i = 0; i < countCheckBoxes; i++)
						objCheckBoxes[i].checked = CheckValue;
			}
		</script>

	<head>
	<body>


		

<table id='ms'  width="100%%" >
<tr>
<td >
<p><a href="http://www.fishpathogens.eu" ><img  alt="Fishpathogens.eu logo" 
style="border:0px;" src="http://www.fishpathogens.eu/images/logo.jpg" ></a></p>
</td>
<td style="text-align:center;vertical-align:middle;">
<img alt="CRL logo"  style="border:0px" width="588px" src="http://www.fishpathogens.eu/images/eurl.png" ></div>
</td>
<td><img alt="Epizone logo" style="border:0px" src="http://www.fishpathogens.eu/images/epizone.gif" >
</td></tr></table>



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
	
	   <form method="%(METHOD)s" action="" id="dvifish" name="dvifish">
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

<div class="form_controls">
%(CONTROLS)s
</div>

<HR class="form_divider" />

%(FIELDS)s


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
