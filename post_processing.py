#!/usr/bin/env python2

import time
import sqlite3
import pybel
import math
import sys

def parse_molecule(molecule):
	# A function to extract the atomic coordinates from an MDL MOL string
	data = molecule.split("\n")
	k = 0
	natom = -1
	coords = []
	for line in data:
		line.strip()
		if (k == 3):
			temp = line.split()
			natom = int(temp[0])
		elif (k > 3):
			if (k <= (3 + natom)):
				temp = line.split()
				coords.append(float(temp[0]))
				coords.append(float(temp[1]))
				coords.append(float(temp[2]))
		k = k + 1
	assert len(coords) == 3*natom
	output = [natom,coords]
	return output

def energy_minimization(raw_molecule):
	output = parse_molecule(raw_molecule)
	natom = output[0]
	raw_coords = output[1]

	mol = pybel.readstring("mol",raw_molecule)
	mol.localopt('mmff94',2500)
	min_molecule = mol.write("mol")
	output = parse_molecule(min_molecule)
	min_coords = output[1]
	assert natom == output[0]    
	rms = 0.0
	for i in range(0,3*natom):
		delta = raw_coords[i] - min_coords[i]
		rms = rms + delta*delta
	rms = math.sqrt(rms/(3.0*natom))
	output = [rms,min_molecule]
	return output

def synthetic_feasibility(molecule):
	# This should be a function which computes the synthetic feasibility of a molecule from 
	# its energy-minimized form, with 0 corresponding to a molecule which is impossible to 
	# create in a lab and 1 denoting a compound that can be purchased off the shelf from a 
	# commercial supplier.
	sigma = 0.5
	return sigma

if (len(sys.argv) != 2):
	print "Usage: ./post-processing.py database file"
	sys.exit(0)
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
		update = "UPDATE Compound SET minimized_structure='" + output[1] + "',root_mean_square=" + str(output[0]) + " WHERE compound_id=" + str(row_id[0]) + ";"
		db.execute(update);
		db.commit()
	print "Done energy minimization"
	q = db.execute(query2)
	result = q.fetchall()
	for row_id in result:
		print "Doing synthetic feasibility on molecule",row_id[0]
		output = synthetic_feasibility(row_id[1])
		update = "UPDATE Compound SET synthetic_feasibility=" + str(output) + " WHERE compound_id=" + str(row_id[0]) + ";"
		db.execute(update);
		db.commit()
	db.close()
	print "Done synthetic feasibility, sleeping..."
	time.sleep(180)

