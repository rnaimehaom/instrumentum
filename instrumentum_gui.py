#!/usr/bin/env python3

import os
import tkinter
from tkinter.filedialog import askopenfilename
from tkinter import messagebox

# One point to note is that this GUI must be run on a monitor with
# resolution at least 800 x 675 for proper viewing of the complete interface. 

class instrumentum:
    def __init__(self,master=None):
        self.master = master
        self.master.title('Instrumentum')
        self.master.geometry('800x675')
        self.master.resizable(0,0)

        global_group = tkinter.LabelFrame(self.master,text="Global",padx=5,pady=5)
        pharma_group = tkinter.LabelFrame(self.master,text="Pharmacophore",padx=5,pady=5)
        skeleton_group = tkinter.LabelFrame(self.master,text="Carbon Skeleton",padx=5,pady=5)
        ratio_group = tkinter.LabelFrame(self.master,text="Rationalization",padx=5,pady=5)
        desat_group = tkinter.LabelFrame(self.master,text="Desaturation & Heteroatoms",padx=5,pady=5)
        masm_group = tkinter.LabelFrame(self.master,text="Molecular Assembly",padx=5,pady=5)

        global_group.grid(row=0,column=0,sticky="ns")
        pharma_group.grid(row=1,column=0,sticky="ns")
        skeleton_group.grid(row=2,column=0,rowspan=2,sticky="ns")

        ratio_group.grid(row=0,column=1,sticky="ns")
        desat_group.grid(row=1,column=1,rowspan=2,sticky="ns")
        masm_group.grid(row=3,column=1,sticky="ns")


        label1 = tkinter.Label(global_group,text='Number of molecules to create:',wraplength=250,justify=tkinter.LEFT)
        label2 = tkinter.Label(global_group,text='Database in which to store molecules:',wraplength=250,justify=tkinter.LEFT)
        label5 = tkinter.Label(skeleton_group,text='Initial Percentage of Nodes to Delete:',wraplength=250,justify=tkinter.LEFT)
        label6 = tkinter.Label(skeleton_group,text='Maximum Number of Initial Deletion Attempts:',wraplength=250,justify=tkinter.LEFT)
        label7 = tkinter.Label(skeleton_group,text='Maximum Number of Secondary Deletion Attempts:',wraplength=250,justify=tkinter.LEFT)
        label8 = tkinter.Label(pharma_group,text='Pharmacophore Radius (in angstroms):',wraplength=250,justify=tkinter.LEFT)
        label9 = tkinter.Label(pharma_group,text='Number of Pharmacophoric Nodes:',wraplength=250,justify=tkinter.LEFT)
        label14 = tkinter.Label(skeleton_group,text='Maximum Number of Quaternary Carbon Atoms for Ending Secondary Deletion:',wraplength=250,justify=tkinter.LEFT)
        label15 = tkinter.Label(skeleton_group,text='Maximum Number of Four Ring Carbons for Ending Secondary Deletion:',wraplength=250,justify=tkinter.LEFT)
        label16 = tkinter.Label(skeleton_group,text='Maximum Number of Rings for Ending Secondary Deletion:',wraplength=250,justify=tkinter.LEFT)
        label17 = tkinter.Label(masm_group,text='Number of Initial Deletion Iterations:',wraplength=250,justify=tkinter.LEFT)
        label18 = tkinter.Label(masm_group,text='Number of Secondary Deletion Iterations:',wraplength=250,justify=tkinter.LEFT)
        label19 = tkinter.Label(masm_group,text='Number of Path Hardening Iterations:',wraplength=250,justify=tkinter.LEFT)
        label20 = tkinter.Label(masm_group,text='Number of Demethylation Iterations:',wraplength=250,justify=tkinter.LEFT)
        label21 = tkinter.Label(masm_group,text='Number of Desaturation/Heteroatom Substitution Iterations:',wraplength=250,justify=tkinter.LEFT)
        label24 = tkinter.Label(global_group,text='Bond Length (in angstroms):',wraplength=250,justify=tkinter.LEFT)
        label45 = tkinter.Label(global_group,text='Parameter Filename:',wraplength=250,justify=tkinter.LEFT)
        label46 = tkinter.Label(global_group,text='Random Number Generator Seed',wraplength=250,justify=tkinter.LEFT)

        chem_label1 = tkinter.Label(ratio_group,text='Percentage of Methyl Groups to Prune:',wraplength=250,justify=tkinter.LEFT)
        chem_label2 = tkinter.Label(ratio_group,text='Minimum Number of Rings for Desaturating Scaffold:',wraplength=250,justify=tkinter.LEFT)
        chem_label3 = tkinter.Label(ratio_group,text='Maximum Number of Rings for Desaturating Scaffold:',wraplength=250,justify=tkinter.LEFT)
        
        self.nmol = tkinter.IntVar()
        self.database = tkinter.StringVar()
        self.parameter_filename = tkinter.StringVar()
        self.percent_methyl = tkinter.DoubleVar()
        self.min_rings = tkinter.IntVar()
        self.max_rings = tkinter.IntVar()
        self.percent = tkinter.DoubleVar()
        self.rng_seed = tkinter.IntVar()
        self.nattempts1 = tkinter.IntVar()
        self.nattempts2 = tkinter.IntVar()
        self.pharm_radius = tkinter.DoubleVar()
        self.blength = tkinter.DoubleVar()
        self.npharm = tkinter.IntVar()
        self.nc4 = tkinter.IntVar()
        self.nrings = tkinter.IntVar()
        self.nc4rings = tkinter.IntVar()
        self.iteration1 = tkinter.IntVar()
        self.iteration2 = tkinter.IntVar()
        self.iteration3 = tkinter.IntVar()
        self.iteration4 = tkinter.IntVar()
        self.iteration5 = tkinter.IntVar()
        self.pharm_harden = tkinter.BooleanVar()
        self.oxy = tkinter.BooleanVar()
        self.nit = tkinter.BooleanVar()
        self.sul = tkinter.BooleanVar()
        self.penta = tkinter.BooleanVar()
        self.amide = tkinter.BooleanVar()
        self.fgrp = tkinter.BooleanVar()
        self.dbond = tkinter.BooleanVar()
        self.tbond = tkinter.BooleanVar()
        self.strip_axial_methyls = tkinter.BooleanVar()

        self.clear_parameters()
        
        chem_entry1 = tkinter.Entry(ratio_group,width=7,textvariable=self.percent_methyl)
        chem_entry2 = tkinter.Entry(ratio_group,width=7,textvariable=self.min_rings)
        chem_entry3 = tkinter.Entry(ratio_group,width=7,textvariable=self.max_rings)
        entry1 = tkinter.Entry(global_group,width=7,textvariable=self.nmol)
        entry5 = tkinter.Entry(skeleton_group,width=7,textvariable=self.percent)
        entry6 = tkinter.Entry(skeleton_group,width=7,textvariable=self.nattempts1)
        entry7 = tkinter.Entry(skeleton_group,width=7,textvariable=self.nattempts2)
        entry8 = tkinter.Entry(pharma_group,width=7,textvariable=self.pharm_radius)
        entry10 = tkinter.Entry(pharma_group,width=7,textvariable=self.npharm)
        entry17 = tkinter.Entry(skeleton_group,width=7,textvariable=self.nc4)
        entry18 = tkinter.Entry(skeleton_group,width=7,textvariable=self.nc4rings)
        entry19 = tkinter.Entry(skeleton_group,width=7,textvariable=self.nrings)
        entry20 = tkinter.Entry(masm_group,width=7,textvariable=self.iteration1)
        entry21 = tkinter.Entry(masm_group,width=7,textvariable=self.iteration2)
        entry22 = tkinter.Entry(masm_group,width=7,textvariable=self.iteration3)
        entry23 = tkinter.Entry(masm_group,width=7,textvariable=self.iteration4)
        entry25 = tkinter.Entry(masm_group,width=7,textvariable=self.iteration5)
        entry24 = tkinter.Entry(global_group,width=7,textvariable=self.blength)
        entry45 = tkinter.Entry(global_group,width=18,textvariable=self.parameter_filename,state=tkinter.DISABLED)
        entry46 = tkinter.Entry(global_group,width=18,textvariable=self.database)
        entry47 = tkinter.Entry(global_group,width=7,textvariable=self.rng_seed)

        axial_methyl_check = tkinter.Checkbutton(desat_group,text='Strip axial methyls from rings?',variable=self.strip_axial_methyls)
        dbond_check = tkinter.Checkbutton(desat_group,text='Create double bonds?',variable=self.dbond,command=self.dbond_change)
        self.tbond_check = tkinter.Checkbutton(desat_group,text='Create triple bonds?',variable=self.tbond)
        oxy_check = tkinter.Checkbutton(desat_group,text='Substitute oxygen atoms?',variable=self.oxy)
        nit_check = tkinter.Checkbutton(desat_group,text='Substitute nitrogen atoms?',variable=self.nit)
        sul_check = tkinter.Checkbutton(desat_group,text='Substitute sulfur atoms?',variable=self.sul)
        penta_check = tkinter.Checkbutton(desat_group,text='Create penta-atomic rings?',variable=self.penta)
        fgrp_check = tkinter.Checkbutton(desat_group,text='Substitute functional groups?',variable=self.fgrp)
        self.amide_check = tkinter.Checkbutton(desat_group,text='Substitute amides, sulfonamides and esters?',variable=self.amide)
        
        label26 = tkinter.Label(skeleton_group,text='Initial Node for Path Hardening is a',wraplength=250,justify=tkinter.LEFT)
        rbutton5 = tkinter.Radiobutton(skeleton_group,text='Random Interior Node',value=0,variable=self.pharm_harden)
        rbutton6 = tkinter.Radiobutton(skeleton_group,text='Pharmacophoric Node',value=1,variable=self.pharm_harden)
        
        button1 = tkinter.Button(text='Save Parameters to Disk',command=self.save_parameters)
        button2 = tkinter.Button(text='Exit',command=root.quit)
        button3 = tkinter.Button(text='Load Parameters from Disk',command=self.load_parameters)
        button4 = tkinter.Button(text='Clear Parameters',command=self.clear_parameters)

        label1.grid(row=0,column=0,sticky=tkinter.W)
        entry1.grid(row=0,column=1,sticky=tkinter.W)
        label24.grid(row=1,column=0,sticky=tkinter.W)
        entry24.grid(row=1,column=1,sticky=tkinter.W)
        label46.grid(row=2,column=0,sticky=tkinter.W)
        entry47.grid(row=2,column=1,sticky=tkinter.W)
        label45.grid(row=3,column=0,sticky=tkinter.W)
        entry45.grid(row=3,column=1,sticky=tkinter.W)
        label2.grid(row=4,column=0,sticky=tkinter.W)
        entry46.grid(row=4,column=1,sticky=tkinter.W)

        label8.grid(row=0,column=0,sticky=tkinter.W)
        entry8.grid(row=0,column=1,sticky=tkinter.W)
        label9.grid(row=1,column=0,sticky=tkinter.W)
        entry10.grid(row=1,column=1,sticky=tkinter.W)

        chem_label1.grid(row=0,column=0,sticky=tkinter.W)
        chem_entry1.grid(row=0,column=1,sticky=tkinter.W)
        chem_label2.grid(row=1,column=0,sticky=tkinter.W)
        chem_entry2.grid(row=1,column=1,sticky=tkinter.W)
        chem_label3.grid(row=2,column=0,sticky=tkinter.W)
        chem_entry3.grid(row=2,column=1,sticky=tkinter.W)

        label17.grid(row=0,column=0,sticky=tkinter.W)
        entry20.grid(row=0,column=1,sticky=tkinter.W)
        label18.grid(row=1,column=0,sticky=tkinter.W)
        entry21.grid(row=1,column=1,sticky=tkinter.W)
        label19.grid(row=2,column=0,sticky=tkinter.W)
        entry22.grid(row=2,column=1,sticky=tkinter.W)
        label20.grid(row=3,column=0,sticky=tkinter.W)
        entry23.grid(row=3,column=1,sticky=tkinter.W)
        label21.grid(row=4,column=0,sticky=tkinter.W)
        entry25.grid(row=4,column=1,sticky=tkinter.W)

        label5.grid(row=0,column=0,sticky=tkinter.W)
        entry5.grid(row=0,column=1,sticky=tkinter.W)
        label6.grid(row=1,column=0,sticky=tkinter.W)
        entry6.grid(row=1,column=1,sticky=tkinter.W)
        label7.grid(row=2,column=0,sticky=tkinter.W)
        entry7.grid(row=2,column=1,sticky=tkinter.W)
        label14.grid(row=3,column=0,sticky=tkinter.W)
        entry17.grid(row=3,column=1,sticky=tkinter.W)
        label15.grid(row=4,column=0,sticky=tkinter.W)
        entry18.grid(row=4,column=1,sticky=tkinter.W)
        label16.grid(row=5,column=0,sticky=tkinter.W)
        entry19.grid(row=5,column=1,sticky=tkinter.W)
        label26.grid(row=6,column=0,sticky=tkinter.W)
        rbutton5.grid(row=6,column=1,sticky=tkinter.W)
        rbutton6.grid(row=7,column=1,sticky=tkinter.W)

        axial_methyl_check.grid(row=0,column=0,sticky=tkinter.W)
        dbond_check.grid(row=1,column=0,sticky=tkinter.W)
        self.tbond_check.grid(row=2,column=0,sticky=tkinter.W)
        self.amide_check.grid(row=3,column=0,sticky=tkinter.W)
        penta_check.grid(row=4,column=0,sticky=tkinter.W)
        fgrp_check.grid(row=5,column=0,sticky=tkinter.W)
        oxy_check.grid(row=6,column=0,sticky=tkinter.W)        
        nit_check.grid(row=7,column=0,sticky=tkinter.W)      
        sul_check.grid(row=8,column=0,sticky=tkinter.W)      

        button1.grid(row=4,column=0)
        button4.grid(row=4,column=1)
        button3.grid(row=5,column=0)
        button2.grid(row=5,column=1)

    def dbond_change(self):
        if (self.dbond.get() == 1):
            self.tbond_check.config(state = tkinter.NORMAL)
            self.amide_check.select()
            self.amide_check.config(state = tkinter.NORMAL)
        else:
            self.tbond_check.deselect()
            self.tbond_check.config(state = tkinter.DISABLED) 
            self.amide_check.deselect()
            self.amide_check.config(state = tkinter.DISABLED)

    def load_parameters(self):
        if (self.parameter_filename.get() == ""):
            filename = askopenfilename(filetypes=(("Text File", "*.txt"),("All Files","*.*")))
            if filename:
            	self.parameter_filename.set(os.path.basename(filename))
            else:
                return
            self.read_parameters()

    def save_parameters(self):
        if (self.parameter_filename.get() == ""):
            filename = tkinter.filedialog.SaveFileDialog(root).go("*.txt")
            if filename:
                self.parameter_filename.set(os.path.basename(filename))
            else:
                return
        self.write_parameters()

    def clear_parameters(self):
        self.percent_methyl.set(35.0)
        self.min_rings.set(1)
        self.max_rings.set(6)
        self.database.set("")
        self.parameter_filename.set("")
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
        
    def convert_boolean(self,bvalue):
        if (bvalue == 0):
            return "No"
        else:
            return "Yes"

    def convert_string(self,val):
        val.upper()
        if (val == "YES"):
            return 1
        else:
            return 0

    def read_parameters(self):
        parameter_filename = self.parameter_filename.get()
        pfile = open(parameter_filename,'r')
        for line in pfile:
            data = line.split('=')
            if (len(data) != 2):
                continue
            name = data[0].strip()
            value = data[1].strip()
            if (name == "NumberMolecules"):
                self.nmol.set(int(value))
            elif (name == "DatabaseFile"):
                self.database.set(value)
            elif (name == "MinimumRings"):
                self.min_rings.set(int(value))
            elif (name == "MaximumRings"):
                self.max_rings.set(int(value))
            elif (name == "InitialPercentage"):
                self.percent.set(100.0*float(value))
            elif (name == "RandomSeed"):
                self.rng_seed.set(int(value))
            elif (name == "BondLength"):
                self.blength.set(float(value))
            elif (name == "MaximumAttempts"):
                self.nattempts1.set(int(value))
            elif (name == "MaximumSecondary"):
                self.nattempts2.set(int(value))
            elif (name == "NumberC4Atoms"):
                self.nc4.set(int(value))
            elif (name == "NumberC4Rings"):
                self.nc4rings.set(int(value))
            elif (name == "NumberPharmacophores"):
                self.npharm.set(int(value))
            elif (name == "PharmacophoreRadius"):
                self.pharm_radius.set(float(value))
            elif (name == "NumberRings"):
                self.nrings.set(int(value))
            elif (name == "NumberInitial"):
                self.iteration1.set(int(value))
            elif (name == "NumberPath"):
                self.iteration2.set(int(value))
            elif (name == "NumberSecondary"):
                self.iteration3.set(int(value))
            elif (name == "NumberRationalize"):
                self.iteration4.set(int(value))
            elif (name == "NumberDesaturate"):
                self.iteration5.set(int(value))
            elif (name == "PharmacophoreHardening"):
                self.pharm_harden.set(self.convert_string(value))
            elif (name == "CreateDoubleBonds"):
                self.dbond.set(self.convert_string(value))
            elif (name == "CreateTripleBonds"):
                self.tbond.set(self.convert_string(value))
            elif (name == "CreateExotic"):
                self.amide.set(self.convert_string(value))
            elif (name == "CreateFiveMemberRings"):
                self.penta.set(self.convert_string(value))
            elif (name == "SubstituteOxygen"):
                self.oxy.set(self.convert_string(value))
            elif (name == "SubstituteNitrogen"):
                self.nit.set(self.convert_string(value))
            elif (name == "SubstituteSulfur"):
                self.sul.set(self.convert_string(value))
            elif (name == "SubstituteFunctionalGroups"):
                self.fgrp.set(self.convert_string(value))
            elif (name == "StripAxialMethyls"):
                self.strip_axial_methyls.set(self.convert_string(value))
            elif (name == "PercentMethyl"):
                self.percent_methyl.set(100.0*float(value))
        pfile.close()        

    def write_parameters(self):
        # Various sanity checks are now required...
        if (self.database.get() is None):
            messagebox.showerror("Illegal Value","The database filename must not be empty!")
            return
        if (self.rng_seed.get() < 0):
            messagebox.showerror("Illegal Value","The random number seed must be non-negative!")
            return
        if (self.percent.get() < 0 or self.percent.get() > 100):
            messagebox.showerror("Illegal Value","The initial percentage parameter must be between 0 and 100!")
            return
        if (self.percent_methyl.get() < 0 or self.percent_methyl.get() > 100):
            messagebox.showerror("Illegal Value","The methyl pruning parameter must be between 0 and 100!")
            return
        if not(self.blength.get() > 0):
            messagebox.showerror("Illegal Value","The bond length must be positive!")
            return
        if (self.min_rings.get() < 0):
            messagebox.showerror("Illegal Value","The minimum ring number must be non-negative!")
            return
        if (self.max_rings.get() < self.min_rings.get()):
            messagebox.showerror("Illegal Value","The maximum ring number must not be less than the minimum ring number!")
            return
        if (self.npharm.get() <= 0):
            messagebox.showerror("Illegal Value","The number of pharmacophores must be positive!")
            return
        if not(self.pharm_radius.get() > 0):
            messagebox.showerror("Illegal Value","The pharmacophore radius must be positive!")
            return
        if (self.nmol.get() <= 0):
            messagebox.showerror("Illegal Value","The number of molecules must be positive!")
            return
        if (self.nc4.get() < 0):
            messagebox.showerror("Illegal Value","The number of quaternary carbons must be non-negative!")
            return
        if (self.nc4rings.get() < 0):
            messagebox.showerror("Illegal Value","The number of four-ring carbons must be non-negative!")
            return
        if (self.nrings.get() < 0):
            messagebox.showerror("Illegal Value","The number of rings must be non-negative!")
            return
        if (self.nattempts1.get() <= 0):
            messagebox.showerror("Illegal Value","The maximum number of initial deletion attempts must be positive!")
            return
        if (self.nattempts2.get() <= 0):
            messagebox.showerror("Illegal Value","The maximum number of secondary deletion attempts must be positive!")
            return
        if (self.iteration1.get() <= 0):
            messagebox.showerror("Illegal Value","The number of initial deletion iterations must be positive!")
            return
        if (self.iteration2.get() <= 0):
            messagebox.showerror("Illegal Value","The number of secondary deletion iterations must be positive!")
            return
        if (self.iteration3.get() <= 0):
            messagebox.showerror("Illegal Value","The number of path hardening iterations must be positive!")
            return
        if (self.iteration4.get() <= 0):
            messagebox.showerror("Illegal Value","The number of rationalization iterations must be positive!")
            return
        if (self.iteration5.get() <= 0):
            messagebox.showerror("Illegal Value","The number of desaturation and heteroatom iterations must be positive!")
            return

        parameter_filename = self.parameter_filename.get()
        pfile = open(parameter_filename,'w')
        pfile.write('DatabaseFile = ' + self.database.get() + '\n')
        pfile.write('RandomSeed = ' + str(self.rng_seed.get()) + '\n')
        pfile.write('InitialPercentage = ' + str(self.percent.get()/100.0) + '\n')
        pfile.write('MaximumAttempts = ' + str(self.nattempts1.get()) + '\n')
        pfile.write('MaximumSecondary = ' + str(self.nattempts2.get()) + '\n')
        pfile.write('BondLength = ' + str(self.blength.get()) + '\n')
        pfile.write('NumberPharmacophores = ' + str(self.npharm.get()) + '\n')
        pfile.write('PharmacophoreRadius = ' + str(self.pharm_radius.get()) + '\n')
        pfile.write('NumberC4Atoms = ' + str(self.nc4.get()) + '\n')
        pfile.write('NumberC4Rings = ' + str(self.nc4rings.get()) + '\n')
        pfile.write('NumberRings = ' + str(self.nrings.get()) + '\n')
        pfile.write('SubstituteOxygen = ' + self.convert_boolean(self.oxy.get()) + '\n')
        pfile.write('SubstituteNitrogen = ' + self.convert_boolean(self.nit.get()) + '\n')
        pfile.write('SubstituteSulfur = ' + self.convert_boolean(self.sul.get()) + '\n')
        pfile.write('SubstituteFunctionalGroups = ' + self.convert_boolean(self.fgrp.get()) + '\n')
        pfile.write('CreateDoubleBonds = ' + self.convert_boolean(self.dbond.get()) + '\n')
        pfile.write('CreateTripleBonds = ' + self.convert_boolean(self.tbond.get()) + '\n')
        pfile.write('CreateExotic = ' + self.convert_boolean(self.amide.get()) + '\n')
        pfile.write('CreateFiveMemberRings = ' + self.convert_boolean(self.penta.get()) + '\n')
        pfile.write('StripAxialMethyls = ' + self.convert_boolean(self.strip_axial_methyls.get()) + '\n')
        pfile.write('NumberInitial = ' + str(self.iteration1.get()) + '\n')
        pfile.write('NumberPath = ' + str(self.iteration2.get()) + '\n')
        pfile.write('NumberSecondary = ' + str(self.iteration3.get()) + '\n')
        pfile.write('NumberRationalize = ' + str(self.iteration4.get()) + '\n')
        pfile.write('NumberDesaturate = ' + str(self.iteration5.get()) + '\n')
        pfile.write('PharmacophoreHardening = ' + self.convert_boolean(self.pharm_harden.get()) + '\n')
        pfile.write('PercentMethyl = ' + str(self.percent_methyl.get()/100.0) + '\n')
        pfile.write('NumberMolecules = ' + str(self.nmol.get()) + '\n')
        pfile.write('MinimumRings = ' + str(self.min_rings.get()) + '\n')
        pfile.write('MaximumRings = ' + str(self.max_rings.get()) + '\n')
        pfile.close()
        
root = tkinter.Tk()
gui = instrumentum(root)
root.mainloop()










