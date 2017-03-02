#!/usr/bin/python

import time
import sqlite3
import sys

def energy_minimization(raw_molecule):
	rms = 17.1
	min_molecule = "FOO BAR"
        output = [rms,min_molecule]
	return output

def synthetic_feasibility(molecule):
	sigma = 0.7765
	return sigma

dbname = sys.argv[1]
rms_cutoff = 10.0
query1 = "SELECT compound_id,raw_structure FROM Compound WHERE minimized_structure IS NULL;"
query2 = "SELECT compound_id,minimized_structure FROM Compound WHERE synthetic_feasibility IS NULL AND root_mean_square < " + str(rms_cutoff) + ";"

while True:
	db = sqlite3.connect(sys.argv[1])
	q = db.execute(query1)
	result = q.fetchall()
	for row_id in result:
		print "Doing energy minimization on molecule",row_id[0]
		output = energy_minimization(row_id[1])
		query = "UPDATE Compound SET minimized_structure='" + output[1] + "',root_mean_square=" + str(output[0]) + " WHERE compound_id=" + str(row_id[0]) + ";"
		db.execute(query);
		db.commit()
	print "Done energy minimization"
	q = db.execute(query2)
	result = q.fetchall()
	for row_id in result:
		print "Doing synthetic feasibility on molecule",row_id[0]
		output = synthetic_feasibility(row_id[1])
		query = "UPDATE Compound SET synthetic_feasibility=" + str(output) + " WHERE compound_id=" + str(row_id[0]) + ";"
		db.execute(query);
		db.commit()
	db.close()
	print "Done synthetic feasibility, sleeping..."
	time.sleep(180)

