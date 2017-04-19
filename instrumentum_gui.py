#!/usr/bin/python

import os
import Tkinter
import FileDialog
import tkMessageBox

# One point to note is that this GUI must be run on a monitor with
# resolution at least 800x600 for proper viewing of the complete interface. 

class instrumentum:
    def __init__(self,master=None):
        self.master = master
        self.master.title('Instrumentum')
        self.master.geometry('800x600')
        self.master.resizable(0,0)

        gparams = Tkinter.LabelFrame(self.master,text="Global",padx=5,pady=5)
        #gparams.pack(padx=10,pady=10)

        label1 = Tkinter.Label(text='Number of molecules to create:',wraplength=250,justify=Tkinter.LEFT)
        label2 = Tkinter.Label(text='Database in which to store molecules:',wraplength=250,justify=Tkinter.LEFT)
        label5 = Tkinter.Label(text='Initial Percentage of Nodes to Delete:',wraplength=250,justify=Tkinter.LEFT)
        label6 = Tkinter.Label(text='Maximum Number of Initial Deletion Attempts:',wraplength=250,justify=Tkinter.LEFT)
        label7 = Tkinter.Label(text='Maximum Number of Secondary Deletion Attempts:',wraplength=250,justify=Tkinter.LEFT)
        label8 = Tkinter.Label(text='Pharmacophore Radius (in angstroms):',wraplength=250,justify=Tkinter.LEFT)
        label9 = Tkinter.Label(text='Number of Pharmacophoric Nodes:',wraplength=250,justify=Tkinter.LEFT)
        label14 = Tkinter.Label(text='Maximum Number of Quaternary Carbon Atoms for Ending Secondary Deletion:',wraplength=250,justify=Tkinter.LEFT)
        label15 = Tkinter.Label(text='Maximum Number of Four Ring Carbons for Ending Secondary Deletion:',wraplength=250,justify=Tkinter.LEFT)
        label16 = Tkinter.Label(text='Maximum Number of Rings for Ending Secondary Deletion:',wraplength=250,justify=Tkinter.LEFT)
        label17 = Tkinter.Label(text='Number of Initial Deletion Iterations:',wraplength=250,justify=Tkinter.LEFT)
        label18 = Tkinter.Label(text='Number of Secondary Deletion Iterations:',wraplength=250,justify=Tkinter.LEFT)
        label19 = Tkinter.Label(text='Number of Path Hardening Iterations:',wraplength=250,justify=Tkinter.LEFT)
        label20 = Tkinter.Label(text='Number of Demethylation Iterations:',wraplength=250,justify=Tkinter.LEFT)
        label21 = Tkinter.Label(text='Number of Desaturation/Heteroatom Substitution Iterations:',wraplength=250,justify=Tkinter.LEFT)
        label24 = Tkinter.Label(text='Bond Length (in angstroms):',wraplength=250,justify=Tkinter.LEFT)
        label45 = Tkinter.Label(text='Parameter Filename:',wraplength=250,justify=Tkinter.LEFT)

        chem_label1 = Tkinter.Label(text='Percentage of Methyl Groups to Prune:',wraplength=250,justify=Tkinter.LEFT)
        chem_label2 = Tkinter.Label(text='Minimum Number of Rings for Desaturating Scaffold:',wraplength=250,justify=Tkinter.LEFT)
        chem_label3 = Tkinter.Label(text='Maximum Number of Rings for Desaturating Scaffold:',wraplength=250,justify=Tkinter.LEFT)
        
        self.percent_methyl = Tkinter.DoubleVar()
        self.min_rings = Tkinter.IntVar()
        self.max_rings = Tkinter.IntVar()
        chem_entry1 = Tkinter.Entry(width=7,textvariable=self.percent_methyl)
        chem_entry2 = Tkinter.Entry(width=7,textvariable=self.min_rings)
        #chem_entry2 = Tkinter.Spinbox(width=7,from_=0,to=10,textvariable=self.min_rings)
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
        self.database = Tkinter.StringVar()
        self.database_display = Tkinter.StringVar() 
        self.parameter_filename = Tkinter.StringVar()
        self.parameter_filename_display = Tkinter.StringVar()

        self.percent = Tkinter.DoubleVar()
        self.rng_seed = Tkinter.IntVar()
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
        self.pharm_harden = Tkinter.BooleanVar()
        self.oxy = Tkinter.BooleanVar()
        self.nit = Tkinter.BooleanVar()
        self.sul = Tkinter.BooleanVar()
        self.penta = Tkinter.BooleanVar()
        self.amide = Tkinter.BooleanVar()
        self.fgrp = Tkinter.BooleanVar()
        self.dbond = Tkinter.BooleanVar()
        self.tbond = Tkinter.BooleanVar()
        self.strip_axial_methyls = Tkinter.BooleanVar()
        entry1 = Tkinter.Entry(width=7,textvariable=self.nmol)
        self.nmol.set(50000)
        self.rng_seed.set(0)
        self.blength.set(1.52)
        self.pharm_radius.set(3.5)
        self.nattempts1.set(100)
        self.nattempts2.set(100)
        self.npharm.set(3)
        self.percent.set(50.0)
        self.nc4.set(2)
        self.nrings.set(4)
        self.nc4rings.set(0)
        self.pharm_harden.set(1)
        self.oxy.set(1)
        self.nit.set(1)
        self.sul.set(0)
        self.amide.set(1)
        self.penta.set(1)
        self.fgrp.set(0)
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
        entry45 = Tkinter.Entry(width=18,textvariable=self.parameter_filename)

        axial_methyl_check = Tkinter.Checkbutton(text='Strip axial methyls from rings?',variable=self.strip_axial_methyls)
        dbond_check = Tkinter.Checkbutton(text='Create double bonds?',variable=self.dbond,command=self.dbond_change)
        self.tbond_check = Tkinter.Checkbutton(text='Create triple bonds?',variable=self.tbond)
        oxy_check = Tkinter.Checkbutton(text='Substitute oxygen atoms?',variable=self.oxy)
        nit_check = Tkinter.Checkbutton(text='Substitute nitrogen atoms?',variable=self.nit)
        sul_check = Tkinter.Checkbutton(text='Substitute sulfur atoms?',variable=self.sul)
        penta_check = Tkinter.Checkbutton(text='Create penta-atomic rings?',variable=self.penta)
        fgrp_check = Tkinter.Checkbutton(text='Substitute functional groups?',variable=self.fgrp)
        self.amide_check = Tkinter.Checkbutton(text='Substitute amides, sulfonamides and esters?',variable=self.amide)
        
        label26 = Tkinter.Label(text='Initial Node for Path Hardening is a',wraplength=250,justify=Tkinter.LEFT)
        rbutton5 = Tkinter.Radiobutton(text='Random Interior Node',value=0,variable=self.pharm_harden)
        rbutton6 = Tkinter.Radiobutton(text='Pharmacophoric Node',value=1,variable=self.pharm_harden)
        
        button1 = Tkinter.Button(text='Write Parameter File',command=self.write_parameters)
        button2 = Tkinter.Button(text='Exit',command=root.quit)
        self.button3 = Tkinter.Button(textvariable=self.database_display,width=16,command=self.get_dbase_file)
        self.button0 = Tkinter.Button(textvariable=self.parameter_filename_display,width=16,command=self.get_parameter_file)

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
        self.tbond_check.grid(row=8,column=3,sticky=Tkinter.W)
        oxy_check.grid(row=9,column=3,sticky=Tkinter.W)
        sul_check.grid(row=10,column=3,sticky=Tkinter.W)
        nit_check.grid(row=11,column=3,sticky=Tkinter.W)
        penta_check.grid(row=12,column=3,sticky=Tkinter.W)
        self.amide_check.grid(row=13,column=3,sticky=Tkinter.W)
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

    def dbond_change(self):
        if (self.dbond.get() == 1):
            self.tbond_check.config(state = Tkinter.NORMAL)
            self.amide_check.select()
            self.amide_check.config(state = Tkinter.NORMAL)
        else:
            self.tbond_check.deselect()
            self.tbond_check.config(state = Tkinter.DISABLED) 
            self.amide_check.deselect()
            self.amide_check.config(state = Tkinter.DISABLED)

    def get_parameter_file(self):
        filename = FileDialog.LoadFileDialog(root).go("*.txt")
        if not(filename is None):
            self.parameter_filename_display.set(os.path.basename(filename))
        self.parameter_filename.set(filename)

    def get_dbase_file(self):
        filename = FileDialog.SaveFileDialog(root).go("*.sqlite3")
        if not(filename is None):
            self.database_display.set(os.path.basename(filename))
        self.database.set(filename)
        
    def convert_boolean(self,bvalue):
        if (bvalue == 0):
            return "No"
        else:
            return "Yes"

    def write_parameters(self):
        # Various sanity checks are now required...
        if (self.parameter_filename.get() is None):
            tkMessageBox.showerror("Illegal Value","The parameter filename must not be empty!")
            return
        if (self.database.get() is None):
            tkMessageBox.showerror("Illegal Value","The database filename must not be empty!")
            return
        if (self.rng_seed.get() < 0):
            tkMessageBox.showerror("Illegal Value","The random number seed must be non-negative!")
            return
        if (self.percent.get() < 0 or self.percent.get() > 100):
            tkMessageBox.showerror("Illegal Value","The initial percentage parameter must be between 0 and 100!")
            return
        if (self.percent_methyl.get() < 0 or self.percent_methyl.get() > 100):
            tkMessageBox.showerror("Illegal Value","The methyl pruning parameter must be between 0 and 100!")
            return
        if not(self.blength.get() > 0):
            tkMessageBox.showerror("Illegal Value","The bond length must be positive!")
            return
        if (self.min_rings.get() < 0):
            tkMessageBox.showerror("Illegal Value","The minimum ring number must be non-negative!")
            return
        if (self.max_rings.get() < self.min_rings.get()):
            tkMessageBox.showerror("Illegal Value","The maximum ring number must not be less than the minimum ring number!")
            return
        if (self.npharm.get() <= 0):
            tkMessageBox.showerror("Illegal Value","The number of pharmacophores must be positive!")
            return
        if not(self.pharm_radius.get() > 0):
            tkMessageBox.showerror("Illegal Value","The pharmacophore radius must be positive!")
            return
        if (self.nmol.get() <= 0):
            tkMessageBox.showerror("Illegal Value","The number of molecules must be positive!")
            return
        if (self.nc4.get() < 0):
            tkMessageBox.showerror("Illegal Value","The number of quaternary carbons must be non-negative!")
            return
        if (self.nc4rings.get() < 0):
            tkMessageBox.showerror("Illegal Value","The number of four-ring carbons must be non-negative!")
            return
        if (self.nrings.get() < 0):
            tkMessageBox.showerror("Illegal Value","The number of rings must be non-negative!")
            return

        parameter_filename = self.parameter_filename.get()
        pfile = open(parameter_filename,'w')
        pfile.write('DatabaseFile = ' + self.database.get() + '\n')
        pfile.write('RandomSeed = ' + str(self.rng_seed) + '\n')
        pfile.write('InitialPercentage = ' + str(self.percent.get()) + '\n')
        pfile.write('MaximumAttempts = ' + str(self.nattempts1.get()) + '\n')
        pfile.write('MaximumSecondary = ' + str(self.nattempts2.get()) + '\n')
        pfile.write('BondLength = ' + str(self.blength.get()) + '\n')
        pfile.write('NumberPharmacophores = ' + self.npharm.get() + '\n')
        pfile.write('PharmacophoreRadius = ' + str(self.pharm_radius.get()) + '\n')
        pfile.write('NumberC4Atoms = ' + str(self.nc4.get()) + '\n')
        pfile.write('NumberC4Rings = ' + str(self.nc4rings.get()) + '\n')
        pfile.write('NumberRings = ' + str(self.nrings.get()) + '\n')
        pfile.write('SubstituteOxygen = ' + convert_boolean(self.oxy.get()) + '\n')
        pfile.write('SubstituteNitrogen = ' + convert_boolean(self.nit.get()) + '\n')
        pfile.write('SubstituteSulfur = ' + convert_boolean(self.sul.get()) + '\n')
        pfile.write('SubstituteFunctionalGroups = ' + convert_boolean(self.fgrp.get()) + '\n')
        pfile.write('CreateDoubleBonds = ' + convert_boolean(self.dbond.get()) + '\n')
        pfile.write('CreateTripleBonds = ' + convert_boolean(self.tbond.get()) + '\n')
        pfile.write('CreateExotic = ' + convert_boolean(self.amide.get()) + '\n')
        pfile.write('CreateFiveMemberRings = ' + convert_boolean(self.penta.get()) + '\n')
        pfile.write('StripAxialMethyls = ' + convert_boolean(self.strip_axial_methyls.get()) + '\n')
        pfile.write('NumberInitial = ' + str(self.iteration1.get()) + '\n')
        pfile.write('NumberPath = ' + str(self.iteration2.get()) + '\n')
        pfile.write('NumberSecondary = ' + str(self.iteration3.get()) + '\n')
        pfile.write('NumberDemethylate = ' + str(self.iteration4.get()) + '\n')
        pfile.write('NumberDesaturate = ' + str(self.iteration5.get()) + '\n')
        pfile.write('PharmacophoreHardening = ' + convert_boolean(self.pharm_harden.get()) + '\n')
        pfile.write('PercentMethyl = ' + str(self.percent_methyl.get()) + '\n')
        pfile.write('NumberMolecules = ' + str(self.nmol.get()) + '\n')
        pfile.write('MinimumRings = ' + str(self.min_rings.get()) + '\n')
        pfile.write('MaximumRings = ' + str(self.max_rings.get()) + '\n')
        pfile.close()
        
root = Tkinter.Tk()
gui = instrumentum(root)
root.mainloop()










