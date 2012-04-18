#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Abstracted database connection.

In the absence of being able to connect to the DVI db, I'll have to trust that
this works.

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


### CONSTANTS & DEFINES ###

### IMPLEMENTATION ###

class DbConn (object):
	def __init__ (self, type, host, user, passwd, name):
		# NOTE: 'name' is the name of the database
		# TODO: have to allow for types
		# TODO: need port?
		
		#import MySQLdb as DbMod
		pass
		
		#self.conn = DbMod.connect (
		#	host = host,
		#	user = user,
		#	passwd = passwd,
		#	db = name
		#)
		
	def select_regions (self):
		"""
		Return list of gene regions.
		"""
		return [
			('1', 'Foo'),
			('2', 'Bar'),
			('3', 'Baz'),
		]
		
	def select_refset (self):
		cursor = conn.cursor()
		# execute SQL statement
		cursor.execute("SELECT id, gene_region from vhsv_gene_regions order by gene_region")
		# get the resultset as a tuple
		results = cursor.fetchall()
		return [(r[0], r[1]) for r in results]
		
	def select_refseqs (self):
		return [
			SeqRecord (id='Foo', name="Foo seq",
				seq=Seq ("ACGTACGTACGTACGTACGCGTACGTACGTACGTACGTACGTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")),
			SeqRecord (id='Bar', name="Bar seq",
				seq=Seq ("ACGTACGCCGTGCGGTGCGTCGGTACGTACGTACGAAATACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACCCGTCGTCAAAAAATACGTACGTACGTACG")),
			SeqRecord (id='Baz', name="Baz seq",
				seq=Seq ("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGAAATACGTACGTACGTCGTACGTACGTACGTACGTACGTACACGTACGTACGTACGTACGTACGTACGTACGTACT")),
			SeqRecord (id='QuuxS', name="Quux seq",
				seq=Seq ("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGAAATACGTACGTACGTCGTACGTACGTACGTACGTACGTACACGTACGTACGTACGTACGTACGTACGTACGTACT")),
			SeqRecord (id='Plugh', name="Phlugh seq",
				seq=Seq ("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGAAATACGTACGTACGTCGTACGTACGTACGTACGTACGTACACGTACGTACGTACGTACGTACGTACGTACGTACT")),
		]



### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


### MAIN ###

if __name__ == '__main__':
	main()


### END ######################################################################
