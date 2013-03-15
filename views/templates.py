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
		<title></title>

		<link rel="shortcut icon" href="http://localhost/" />

		<link rel="stylesheet" type="text/css" href="http://www.fishpathogens.eu/media/css/matchseq.css" />
		<script type="text/javascript" src="http://www.fishpathogens.eu/media/js/raphael-min.js" ></script> 
		<script type="text/javascript" src="http://www.fishpathogens.eu/media/js/jsphylosvg-min.js"></script>
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
		
		<script  src="/fp-common/jscript/script.js" type="text/javascript"></script>
<script src="/fp-common/jscript/menu.js" type="text/javascript"></script>
<script  src="/fp-common/jscript/ajax.js" type="text/javascript"></script>
<script  src="/fp-common/jscript/form.js" type="text/javascript"></script>
<script  src="/fp-common/jscript/geo.js" type="text/javascript"></script>
<LINK href="/fp-common/css/menu.css" rel="stylesheet" type="text/css">
<LINK href="/fp-common/css/style.css" rel="stylesheet" type="text/css">

 <script src="http://code.jquery.com/jquery.js"></script>
    <script src="/fp-common/bootstrap/js/bootstrap.min.js"></script>			
<link href="/fp-common/bootstrap/css/bootstrap-cerulean.min.css" rel="stylesheet" media="screen">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">

  
		

	<head>
	<body>


		
	<div class="navbar  navbar-fixed-top">
			
      <div class="navbar-inner">
        <div class="container-fluid">
          <button type="button" class="btn btn-navbar" data-toggle="collapse" data-target=".nav-collapse">
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
          </button>
      <!--    <a class="brand" ></a>-->
      
          <div class="nav-collapse collapse">
   
            <ul class="nav">
            
              <li class="dropdown"><a href="#" class="dropdown-toggle"  data-toggle="dropdown"><i class="icon-home icon-white"></i> Home</a>
              <ul class="dropdown-menu">
                <li><a href="/">Main</a></li>
  				  <li><a href="/fp-vhsv/">VHSV</a></li>
  				 <li><a href="/fp-ihnv/">IHNV</a></li>
                <li><a href="/fp-betanodavirus/">Betanodavirus</a></li>
   </ul>
   
            </li>
              <li><a href="%(SIDEBAR_PATHOGEN_STR)sindex.php#about">About</a></li>
              <li><a href="%(SIDEBAR_PATHOGEN_STR)sindex.php#contact"><i class="icon-envelope icon-white"></i> Contact</a></li>
            </ul>
          </div><!--/.nav-collapse -->
        </div>
      </div>
    </div>
			<br/><Br/><br/>
	<div class="container-fluid">
		
      <div class="row-fluid">
			
        <div class="span3">
			<p align='center'>
			<img  alt="Fishpathogens.eu logo" 
style="border:0px;"  src="/fp-common/images/logo.jpg" ></a></p>	
          <div class="well sidebar-nav">
            <ul class="nav nav-list">
              <li class="nav-header">Isolates</li>
              <li><a href="%(SIDEBAR_PATHOGEN_STR)sreports_browse.php?type=1">Isolate reports</a></li>
			<li><a href="%(SIDEBAR_PATHOGEN_STR)sreports_browse.php?type=2">Sequence reports</a></li>

			 <li class="nav-header">Search Reports</li>
			<li><a href="%(SIDEBAR_PATHOGEN_STR)ssearch.php">Search with a text query <i class='icon-search icon-white'></i></a></li>

<li><a href="%(SIDEBAR_PATHOGEN_STR)ssearch_blast.php">Search with a sequence</a></li>
<li><a href="/cgi-bin/matchseq/matchseq.cgi?submit=1%%3A+Select+regions&pathogen=%(SIDEBAR_PATHOGEN_ID)s&dbid=%(SIDEBAR_PATHOGEN_STR)s">Sequence matcher</a></li>
			<li class="nav-header">Add a Report</li>
<li><a href="%(SIDEBAR_PATHOGEN_STR)sinput_isolate.php">Add an Isolate report</a></li>
<li><a href="%(SIDEBAR_PATHOGEN_STR)sinput_sequence.php">Add a Sequence report</a></li>
              
             
             
              <li class="nav-header">Related Information</li>
              
<li><a href="%(SIDEBAR_PATHOGEN_STR)sg_links.php">Links to external resources</a></li>
<li><a href="/%(SIDEBAR_PATHOGEN_STR)sg_ncbi_resources.php">Records for VHSV at the NCBI</a></li>
            </ul>
          </div><!--/.well -->
		
		<img alt="CRL logo"  style="border:0px" width="294px" src="/fp-common/images/eurl.png" >
		<br/>
		<img alt="Epizone logo" style="border:0px" src="/fp-common/images/epizone.gif" >
		<br/>	
		
		</div><!--/span-->
	<div class="span9">


<h2>Sequence matcher</h2>
	
		<p class="description">
			This webservice accepts an input sequence and attempts to match it
			against selection of database sequences, using either FASTA matching
			or by building a phylogeny using Neighbour-joining.
		</p>
	
		%(MESSAGES)s
	
		%(RESULTS)s
	
		
	
	   <form method="%(METHOD)s" class='form' action="" id="dvifish" name="dvifish">
			%(FORM)s
		</form>
	
		
		
		 <hr>

      <footer>
        <p align='center'> </p>
      </footer>

    </div>
		
	</body>
</html>
"""

# the layout for a form
FORM_BODY = """
<H3>%(TITLE)s</H3>
<p class="description">%(DESC)s</p>

<div class="form_controls">
%(CONTROLS)s
</div>


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
