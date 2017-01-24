#!/usr/bin/python

import Tkinter
import FileDialog
import random
import os
import sys

# One point to note is that this GUI must be run on a monitor with
# resolution at least 800x600 for proper viewing of the complete interface. 

class masm:
    def __init__(self,master=None):
        self.master = master
        self.master.title('Molecular Assembler, version 4.0')
        self.master.geometry('800x500')
        self.master.resizable(0,0)
        if (os.name == 'nt'):
           self.hdir = "H:\\"
        elif (os.name == 'posix'):
           self.hdir = "~/"

        label1 = Tkinter.Label(text='Number of molecules to create:',wraplength=250,justify=Tkinter.LEFT)
        label2 = Tkinter.Label(text='Filename in which to store molecules:',wraplength=250,justify=Tkinter.LEFT)
        label5 = Tkinter.Label(text='Percentage of Initial Nodes to Delete:',wraplength=250,justify=Tkinter.LEFT)
        label6 = Tkinter.Label(text='Maximum Number of Initial Deletion Attempts:',wraplength=250,justify=Tkinter.LEFT)
        label7 = Tkinter.Label(text='Maximum Number of Secondary Deletion Attempts:',wraplength=250,justify=Tkinter.LEFT)
        label8 = Tkinter.Label(text='Radius of Steric Volume:',wraplength=250,justify=Tkinter.LEFT)
        label9 = Tkinter.Label(text='Number of Pharmacophoric Nodes:',wraplength=250,justify=Tkinter.LEFT)
        label14 = Tkinter.Label(text='Maximum Number of Quaternary Carbon Atoms for Ending Secondary Deletion:',wraplength=250,justify=Tkinter.LEFT)
        label15 = Tkinter.Label(text='Maximum Number of Four Ring Carbons for Ending Secondary Deletion:',wraplength=250,justify=Tkinter.LEFT)
        label16 = Tkinter.Label(text='Maximum Number of Rings for Ending Secondary Deletion:',wraplength=250,justify=Tkinter.LEFT)
        label17 = Tkinter.Label(text='Number of Initial Deletion Iterations:',wraplength=250,justify=Tkinter.LEFT)
        label18 = Tkinter.Label(text='Number of Secondary Deletion Iterations:',wraplength=250,justify=Tkinter.LEFT)
        label19 = Tkinter.Label(text='Number of Path Hardening Iterations:',wraplength=250,justify=Tkinter.LEFT)
        label20 = Tkinter.Label(text='Number of Demethylation Iterations:',wraplength=250,justify=Tkinter.LEFT)
        label21 = Tkinter.Label(text='Number of Desaturation/Heteroatom Substitution Iterations:',wraplength=250,justify=Tkinter.LEFT)
        label24 = Tkinter.Label(text='Internode Distance:',wraplength=250,justify=Tkinter.LEFT)
        label45 = Tkinter.Label(text='Parameter Filename:',wraplength=250,justify=Tkinter.LEFT)

        chem_label1 = Tkinter.Label(text='Percentage of Methyl Groups to Prune:',wraplength=250,justify=Tkinter.LEFT)
        chem_label2 = Tkinter.Label(text='Minimum Number of Rings for Desaturating Scaffold:',wraplength=250,justify=Tkinter.LEFT)
        chem_label3 = Tkinter.Label(text='Maximum Number of Rings for Desaturating Scaffold:',wraplength=250,justify=Tkinter.LEFT)
        
        self.percent_methyl = Tkinter.DoubleVar()
        self.min_rings = Tkinter.IntVar()
        self.max_rings = Tkinter.IntVar()
        chem_entry1 = Tkinter.Entry(width=7,textvariable=self.percent_methyl)
        chem_entry2 = Tkinter.Entry(width=7,textvariable=self.min_rings)
        chem_entry3 = Tkinter.Entry(width=7,textvariable=self.max_rings)

        self.percent_methyl.set(35.0)
        self.min_rings.set(1)
        self.max_rings.set(6)
        
        chem_label1.grid(row=0,column=3,sticky=Tkinter.W)
        chem_entry1.grid(row=0,column=4,sticky=Tkinter.W)
        chem_label2.grid(row=14,sticky=Tkinter.W)
        chem_entry2.grid(row=14,column=1,sticky=Tkinter.W)
        chem_label3.grid(row=15,sticky=Tkinter.W)
        chem_entry3.grid(row=15,column=1,sticky=Tkinter.W)
        
        self.nmol = Tkinter.IntVar()
        self.ofilename = Tkinter.StringVar()
        self.ofilename_display = Tkinter.StringVar() 
        self.parameter_fname = Tkinter.StringVar()
        self.parameter_fname_display = Tkinter.StringVar()

        self.percent = Tkinter.DoubleVar()
        self.nattempts1 = Tkinter.IntVar()
        self.nattempts2 = Tkinter.IntVar()
        self.pharm_radius = Tkinter.DoubleVar()
        self.blength = Tkinter.DoubleVar()
        self.npharm = Tkinter.IntVar()
        self.nc4 = Tkinter.IntVar()
        self.nrings = Tkinter.IntVar()
        self.nc4rings = Tkinter.IntVar()
        self.iteration1 = Tkinter.IntVar()
        self.iteration2 = Tkinter.IntVar()
        self.iteration3 = Tkinter.IntVar()
        self.iteration4 = Tkinter.IntVar()
        self.iteration5 = Tkinter.IntVar()
        self.oxy = Tkinter.IntVar()
        self.nit = Tkinter.IntVar()
        self.sul = Tkinter.IntVar()
        self.penta = Tkinter.IntVar()
        self.amide = Tkinter.IntVar()
        self.fgrp = Tkinter.IntVar()
        self.dbond = Tkinter.IntVar()
        self.tbond = Tkinter.IntVar()
        self.uniqueness = Tkinter.IntVar()
        self.ligand = Tkinter.IntVar()
        self.strip_axial_methyls = Tkinter.IntVar()
        entry1 = Tkinter.Entry(width=7,textvariable=self.nmol)
        self.nmol.set(50000)
        self.blength.set(1.52)
        self.pharm_radius.set(3.5)
        self.nattempts1.set(100)
        self.nattempts2.set(100)
        self.npharm.set(3)
        self.percent.set(50.0)
        self.nc4.set(2)
        self.nrings.set(4)
        self.nc4rings.set(0)
        self.oxy.set(1)
        self.nit.set(1)
        self.sul.set(0)
        self.amide.set(1)
        self.penta.set(1)
        self.fgrp.set(0)
        self.uniqueness.set(1)
        self.dbond.set(1)
        self.tbond.set(0)
        self.iteration1.set(3)
        self.iteration2.set(3)
        self.iteration3.set(4)
        self.iteration4.set(3)
        self.iteration5.set(5)
        self.strip_axial_methyls.set(1)
        
        entry5 = Tkinter.Entry(width=7,textvariable=self.percent)
        entry6 = Tkinter.Entry(width=7,textvariable=self.nattempts1)
        entry7 = Tkinter.Entry(width=7,textvariable=self.nattempts2)
        entry8 = Tkinter.Entry(width=7,textvariable=self.pharm_radius)
        entry10 = Tkinter.Entry(width=7,textvariable=self.npharm)
        entry17 = Tkinter.Entry(width=7,textvariable=self.nc4)
        entry18 = Tkinter.Entry(width=7,textvariable=self.nc4rings)
        entry19 = Tkinter.Entry(width=7,textvariable=self.nrings)
        entry20 = Tkinter.Entry(width=7,textvariable=self.iteration1)
        entry21 = Tkinter.Entry(width=7,textvariable=self.iteration2)
        entry22 = Tkinter.Entry(width=7,textvariable=self.iteration3)
        entry23 = Tkinter.Entry(width=7,textvariable=self.iteration4)
        entry25 = Tkinter.Entry(width=7,textvariable=self.iteration5)
        entry24 = Tkinter.Entry(width=7,textvariable=self.blength)
        entry45 = Tkinter.Entry(width=18,textvariable=self.parameter_fname)

        axial_methyl_check = Tkinter.Checkbutton(text='Strip axial methyls from rings?',variable=self.strip_axial_methyls)
        dbond_check = Tkinter.Checkbutton(text='Create double bonds?',variable=self.dbond)
        tbond_check = Tkinter.Checkbutton(text='Create triple bonds?',variable=self.tbond)
        oxy_check = Tkinter.Checkbutton(text='Substitute oxygen atoms?',variable=self.oxy)
        nit_check = Tkinter.Checkbutton(text='Substitute nitrogen atoms?',variable=self.nit)
        sul_check = Tkinter.Checkbutton(text='Substitute sulfur atoms?',variable=self.sul)
        penta_check = Tkinter.Checkbutton(text='Create penta-atomic rings?',variable=self.penta)
        fgrp_check = Tkinter.Checkbutton(text='Substitute functional groups?',variable=self.fgrp)
        amide_check = Tkinter.Checkbutton(text='Substitute amides, sulfonamides and esters?',variable=self.amide)
        
        self.hchoice = Tkinter.IntVar()
        label26 = Tkinter.Label(text='Initial Node for Path Hardening is a',wraplength=250,justify=Tkinter.LEFT)
        rbutton5 = Tkinter.Radiobutton(text='Random Interior Node',value=0,variable=self.hchoice)
        rbutton6 = Tkinter.Radiobutton(text='Pharmacophoric Node',value=1,variable=self.hchoice)
        self.hchoice.set(1)
        
        button1 = Tkinter.Button(text='Write Parameter File',command=self.startjob)
        button2 = Tkinter.Button(text='Exit',command=root.quit)
        self.button3 = Tkinter.Button(textvariable=self.ofilename_display,width=16,command=self.get_sdfile1)
        self.button0 = Tkinter.Button(textvariable=self.parameter_fname_display,width=16,command=self.get_sdfile0)

        label1.grid(row=0,sticky=Tkinter.W)
        entry1.grid(row=0,column=1,sticky=Tkinter.W)
        label45.grid(row=1,sticky=Tkinter.W)
        self.button0.grid(row=1,column=1,sticky=Tkinter.W)
        label2.grid(row=2,sticky=Tkinter.W)
        self.button3.grid(row=2,column=1,sticky=Tkinter.W)

        label24.grid(row=3,sticky=Tkinter.W)
        entry24.grid(row=3,column=1,sticky=Tkinter.W)
        label5.grid(row=4,sticky=Tkinter.W)
        entry5.grid(row=4,column=1,sticky=Tkinter.W)
        label6.grid(row=5,sticky=Tkinter.W)
        entry6.grid(row=5,column=1,sticky=Tkinter.W)
        label7.grid(row=6,sticky=Tkinter.W)
        entry7.grid(row=6,column=1,sticky=Tkinter.W)
        label8.grid(row=7,sticky=Tkinter.W)
        entry8.grid(row=7,column=1,sticky=Tkinter.W)
        label9.grid(row=8,sticky=Tkinter.W)
        entry10.grid(row=8,column=1,sticky=Tkinter.W)

        label26.grid(row=9,column=0,sticky=Tkinter.W)
        rbutton5.grid(row=9,column=1,sticky=Tkinter.W)
        rbutton6.grid(row=10,column=1,sticky=Tkinter.W)
        axial_methyl_check.grid(row=6,column=3,sticky=Tkinter.W)
        dbond_check.grid(row=7,column=3,sticky=Tkinter.W)
        tbond_check.grid(row=8,column=3,sticky=Tkinter.W)
        oxy_check.grid(row=9,column=3,sticky=Tkinter.W)
        sul_check.grid(row=10,column=3,sticky=Tkinter.W)
        nit_check.grid(row=11,column=3,sticky=Tkinter.W)
        penta_check.grid(row=12,column=3,sticky=Tkinter.W)
        amide_check.grid(row=13,column=3,sticky=Tkinter.W)
        fgrp_check.grid(row=14,column=3,sticky=Tkinter.W)

        label14.grid(row=11,sticky=Tkinter.W)
        entry17.grid(row=11,column=1,sticky=Tkinter.W)
        label15.grid(row=12,sticky=Tkinter.W)
        entry18.grid(row=12,column=1,sticky=Tkinter.W)
        label16.grid(row=13,sticky=Tkinter.W)
        entry19.grid(row=13,column=1,sticky=Tkinter.W)
        label17.grid(row=1,column=3,sticky=Tkinter.W)
        entry20.grid(row=1,column=4,sticky=Tkinter.W)
        label18.grid(row=2,column=3,sticky=Tkinter.W)
        entry21.grid(row=2,column=4,sticky=Tkinter.W)
        label19.grid(row=3,column=3,sticky=Tkinter.W)
        entry22.grid(row=3,column=4,sticky=Tkinter.W)
        label20.grid(row=4,column=3,sticky=Tkinter.W)
        entry23.grid(row=4,column=4,sticky=Tkinter.W)
        label21.grid(row=5,column=3,sticky=Tkinter.W)
        entry25.grid(row=5,column=4,sticky=Tkinter.W)
        button1.grid(row=31,column=0)
        button2.grid(row=31,column=3)

    def get_sdfile0(self):
        filename = FileDialog.LoadFileDialog(root).go(self.hdir,"*.txt")
        if not(filename is None):
           if (self.hdir == "H:\\"):
              temp = filename.split("\\")
              pure_name = temp[len(temp)-1].split(".")
              self.parameter_fname_display.set(pure_name[0])
           elif (self.hdir == "~/"):
              temp = filename.split("/")
              pure_name = temp[len(temp)-1].split(".")
              self.parameter_fname_display.set(pure_name[0])
        self.parameter_fname.set(filename)

    def get_sdfile1(self):
        filename = FileDialog.SaveFileDialog(root).go(self.hdir,"*.sdf")
        if not(filename is None):
           if (self.hdir == "H:\\"):
              temp = filename.split("\\")
              pure_name = temp[len(temp)-1].split(".")
              self.ofilename_display.set(pure_name[0])
           elif (self.hdir == "~/"):
              temp = filename.split("/")
              pure_name = temp[len(temp)-1].split(".")
              self.ofilename_display.set(pure_name[0])
        self.ofilename.set(filename)
        
    def startjob(self):
        para_filename = self.parameter_fname.get()
        pfile = open(para_filename,'w')
        pfile.write('output_file = ' + ofilename.get() + '\n')
        pfile.write('percent = ' + str(self.percent.get()) + '\n')
        pfile.write('max_attempts = ' + str(self.nattempts1.get()) + '\n')
        pfile.write('max_secondary = ' + str(self.nattempts2.get()) + '\n')
        pfile.write('bond_length = ' + str(self.blength.get()) + '\n')
        pfile.write('npharm = ' + self.npharm.get() + '\n')
        pfile.write('pharm_radius = ' + str(self.pharm_radius.get()) + '\n')
        pfile.write('nc4 = ' + str(self.nc4.get()) + '\n')
        pfile.write('nc4rings = ' + str(self.nc4rings.get()) + '\n')
        pfile.write('nrings = ' + str(self.nrings.get()) + '\n')
        pfile.write('subs_oxy = ' + str(self.oxy.get()) + '\n')
        pfile.write('subs_nit = ' + str(self.nit.get()) + '\n')
        pfile.write('subs_sul = ' + str(self.sul.get()) + '\n')
        pfile.write('subs_fun = ' + str(self.fgrp.get()) + '\n')
        pfile.write('create_double = ' + str(self.dbond.get()) + '\n')
        pfile.write('create_triple = ' + str(self.tbond.get()) + '\n')
        pfile.write('create_exotic = ' + str(self.amide.get()) + '\n')
        pfile.write('create_penta = ' + str(self.penta.get()) + '\n')
        pfile.write('strip_axial_methyls = ' + str(self.strip_axial_methyls.get()) + '\n')
        pfile.write('n_initial = ' + str(self.iteration1.get()) + '\n')
        pfile.write('n_path = ' + str(self.iteration2.get()) + '\n')
        pfile.write('n_secondary = ' + str(self.iteration3.get()) + '\n')
        pfile.write('n_demethylate = ' + str(self.iteration4.get()) + '\n')
        pfile.write('n_desaturate = ' + str(self.iteration5.get()) + '\n')
        pfile.write('pharm_hardening = ' + str(self.hchoice.get()) + '\n')
        pfile.write('percent_methyl = ' + str(self.percent_methyl.get()) + '\n')
        pfile.write('n_mols = ' + str(self.nmol.get()) + '\n')
        pfile.write('min_rings = ' + str(self.min_rings.get()) + '\n')
        pfile.write('max_rings = ' + str(self.max_rings.get()) + '\n')
        pfile.close()
        
root = Tkinter.Tk()
MA = masm(root)
root.mainloop()










