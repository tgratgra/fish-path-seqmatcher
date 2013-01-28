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

# map virus type to name


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

	def get_tables_by_pathogen_name (self, name):
		"""
		Map a virus name to the names of the two tables used.
		"""
		# NOTE: most virus names map straight to prefix
		# TODO: what does betanodavirus map to?
		clean_name = name.upper()
		prefix = {
			'BETANODAVIRUS': 'betanodavirus',
		}.get (clean_name, clean_name) 
		return (
			'%s_gene_regions' % prefix,
			'%s_sequence' % prefix,
		)

	def get_tables_by_pathogen_id (self, id):
		cursor = self.conn.cursor()
		cursor.execute("SELECT identifier from dbdesc WHERE id = %s" % id)
		name = cursor.fetchone()[0]
		return self.get_tables_by_pathogen_name (name)
		
	def select_pathogens (self):
		"""
		Return list of all pathogens.
		"""
		cursor = self.conn.cursor()
		cursor.execute("SELECT * FROM dbdesc where id in (1,3,6)")
		results = [{'id': r[0], 'title': r[1], 'desc': r[2]} for r in cursor.fetchall()]
		return results
		
	def select_regions (self, path_id):
		"""
		Return list of all gene regions.
		"""
		reg_table, seq_table = self.get_tables_by_pathogen_id (path_id)
		cursor = self.conn.cursor()
		cursor.execute("SELECT id, gene_region from %s order by gene_region" % reg_table)
		results = [{'id': r[0], 'gene_region': r[1]} for r in cursor.fetchall()]
		return results
	
	def select_regions_by_ids (self, path_id, ids):
		reg_table, seq_table = self.get_tables_by_pathogen_id (path_id)
		cursor = self.conn.cursor()
		qry_tmpl = "SELECT id, gene_region from %s WHERE id IN (%s) order by gene_region"
		qry = qry_tmpl % (reg_table, ', '.join (["%s" % x for x in ids]))
		cursor.execute (qry)
		results = [{'id': r[0], 'gene_region': r[1]} for r in cursor.fetchall()]
		return results
		
	def select_seqs_by_regions (self, path_id, regions):
		"""
		Return the sequence information associated with every region.
		
		Note it the the region name, not the reghion id
		"""
		reg_table, seq_table = self.get_tables_by_pathogen_id (path_id)
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
		qry_tmpl = "SELECT %s FROM %s WHERE gene_region LIKE '%%;%s;;%%' LIMIT %s" 
		
		# for each region do query
		records = []
		for gr in regions:
			qry = qry_tmpl % (
				field_str,
				seq_table,
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

	def select_seqs_by_ids (self, path_id, ids):
		"""
		Return sequence info for subsequent alignment.
		"""
		reg_table, seq_table = self.get_tables_by_pathogen_id (path_id)
		fetched_fields = [
			'id',
			'isolate_name',
			'accession_number',
			'nucleotide_sequence',
			'gene_region',
		]
		field_str = ', '.join (fetched_fields)
		qry_tmpl = "SELECT %s from %s WHERE id IN (%s) LIMIT %s"
		qry = qry_tmpl % (
			field_str,
			seq_table,
			', '.join (["%s" % x for x in ids]),
			config.SEQLIMIT,
		)
		
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
