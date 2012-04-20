#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Abstracted database connection.

In the absence of being able to connect to the DVI db, I'll have to trust that
this works.

"""

__docformat__ = 'restructuredtext en'


### IMPORTS ###

from exceptions import ValueError

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import config


### CONSTANTS & DEFINES ###

### IMPLEMENTATION ###

class DbConn (object):
	def __init__ (self, protocol, host, user, passwd, name):
		# NOTE: 'name' is the name of the database
		# TODO: need port?
		
		# NOTE: the whole types things is mainly so we can test on Sqlite
		self.protocol = protocol.lower()
		if self.protocol == 'mysql':
			import MySQLdb as DbMod
			self.conn = DbMod.connect (
				host = host,
				user = user,
				passwd = passwd,
				db = name
			)
		elif self.protocol == 'sqlite':
			import sqlite3
			self.conn = sqlite3.connect(host)
		else:
			raise ValueError, "unrecognised db protocol '%s'" % self.protocol

		
	def select_regions (self):
		"""
		Return list of all gene regions.
		"""
		cursor = self.conn.cursor()
		cursor.execute("SELECT id, gene_region from VHSV_gene_regions order by gene_region")
		results = [{'id': r[0], 'gene_region': r[1]} for r in cursor.fetchall()]
		return results
	
	def select_regions_by_ids (self, ids):
		cursor = self.conn.cursor()
		qry_tmpl = "SELECT id, gene_region from VHSV_gene_regions WHERE id IN (%s) order by gene_region"
		qry = qry_tmpl % ', '.join (["%s" % x for x in ids])
		cursor.execute (qry)
		results = [{'id': r[0], 'gene_region': r[1]} for r in cursor.fetchall()]
		return results
		
	def select_seqs_by_regions (self, regions):
		"""
		Return the sequence information associated with every region.
		
		Note it the the region name, not the reghion id
		"""
		# NOTE: regions is a list, must call for each region
		cursor = self.conn.cursor()
		fetched_fields = [
			'id',
			'isolate_name',
			'accession_number',
			'sequencing_method',
			'author',
		]
		field_str = ', '.join (fetched_fields)
		qry_tmpl = "SELECT %s FROM VHSV_sequence WHERE gene_region LIKE '%%;%s;;%%' LIMIT %s"
		
		# for each region do query
		records = []
		for gr in regions:
			qry = qry_tmpl % (
				field_str,
				gr,
				config.SEQLIMIT,
			)
			cursor.execute (qry)
			region_re = cursor.fetchall()
			records.extend (region_re)
			
		# process records to coherent form
		results = []
		for r in records:
			new_res = {}
			for i in range (len (fetched_fields)):
				k = fetched_fields[i]
				v = r[i]
				new_res[k] = v
			results.append (new_res)
			
		return results

	def select_seqs_by_ids (self, ids):
		"""
		Return sequence info for subsequent alignment.
		"""
		fetched_fields = [
			'id',
			'isolate_name',
			'accession_number',
			'nucleotide_sequence',
			'gene_region',
		]
		field_str = ', '.join (fetched_fields)
		qry_tmpl = "SELECT %s from VHSV_sequence  WHERE id IN (%s) LIMIT %s"
		qry = qry_tmpl % (
			field_str,
			', '.join (["%s" % x for x in ids]),
			config.SEQLIMIT,
		)
		print qry
		
		cursor = self.conn.cursor()
		cursor.execute (qry)
		records = cursor.fetchall()
		
		# process records to coherent form
		results = []
		for r in records:
			new_res = {}
			for i in range (len (fetched_fields)):
				k = fetched_fields[i]
				v = r[i]
				new_res[k] = v
			results.append (new_res)
			
		return results


### TEST & DEBUG ###

def _doctest ():
	import doctest
	doctest.testmod ()


def test():
	conn = DbConn('sqlite', 'fish.sqlite', None, None, None)
	print conn.select_regions()
	print conn.select_regions_by_ids ([15, 19, 11])
	print conn.select_seqs_by_regions (['G-gene ORF', 'P (M1)-gene ORF', ])
	print conn.select_seqs_by_ids ([166, 169, 28, 65, 196])


### MAIN ###

if __name__ == '__main__':
	test()


### END ######################################################################
