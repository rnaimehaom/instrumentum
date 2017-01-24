#!/usr/local/bin/python

import sys
import os

ffile = open(sys.argv[1],"r")
data = ffile.readlines()
ffile.close()
kount = 1
mol = ''
for i in range(0,len(data)):
    mol = mol + data[i]
    if (data[i] == '$$$$\n'):
       file2 = open("temp.mol","w")
       file2.write(mol)
       file2.close()
       os.system("babel temp.mol temp" + str(kount) + ".gpr");
       mol = ''
       kount = kount + 1
os.system("rm temp.mol")
 
