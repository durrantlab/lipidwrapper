#!/usr/bin/python

''' Copyright (c) 2014, Jacob D. Durrant
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.'''

import os

try:
    from Tkinter import *
    import tkFileDialog
    import tkMessageBox
except:
    print "\nERROR: Unfortunately, it seems you don't have the python module Tkinter"
    print "installed on your computer, so the graphical user interface is not"
    print "available. Consider using the command-line program, lipidwrapper.py,"
    print "instead.\n"
    sys.exit()
    
class LipidWrapper(Frame):
    def __init__(self,Master=None,**kw):
        # set up the GUI

        apply(Frame.__init__,(self,Master),kw)
        self.__Frame2 = Frame(self)
        self.__Frame2.pack(side='top')
        self.__Label1 = Label(self.__Frame2,background='#000000'
            ,foreground='White',padx=2,pady=2
            ,text='Initial Planar Bilayer Model',width=50)
        self.__Label1.pack(pady=5,side='top')
        self.__Frame1 = Frame(self,height=250,width=700)
        self.__Frame1.pack(side='top')
        self.__Frame8 = Frame(self)
        self.__Frame8.pack(pady=5,side='top')
        self.__Frame10 = Frame(self)
        self.__Frame10.pack(side='top')
        self.__Label5 = Label(self.__Frame10,background='#000000'
            ,foreground='White',padx=2,pady=2,text='Create the Surface Mesh'
            ,width=50)
        self.__Label5.pack(pady=5,side='top')
        self.__Frame11 = Frame(self)
        self.__Frame11.pack(side='top')
        self.__Frame6 = Frame(self)
        self.__Frame6.pack(pady=5,side='top')
        self.__Frame3 = Frame(self)
        self.__Frame3.pack(pady=5,side='top')
        self.__Label7 = Label(self.__Frame3,background='Black'
            ,foreground='White',padx=2,pady=2,text='Model Dimensions',width=50)
        self.__Label7.pack(side='top')
        self.__Frame12 = Frame(self,borderwidth=1,height=95,relief='sunken'
            ,width=725)
        self.__Frame12.pack(side='top')
        self.__Frame13 = Frame(self)
        self.__Frame13.pack(pady=5,side='top')
        self.__Frame41 = Frame(self)
        self.__Frame41.pack(side='top')
        self.__Label6 = Label(self.__Frame41,background='Black'
            ,foreground='White',padx=2,pady=2,text='Resolving Lipid Clashes'
            ,width=50)
        self.__Label6.pack(pady=5,side='top')
        self.__Frame42 = Frame(self,height=100,width=750)
        self.__Frame42.pack(padx=8,side='top')
        self.__Frame40 = Frame(self)
        self.__Frame40.pack(pady=5,side='top')
        self.__Frame68 = Frame(self)
        self.__Frame68.pack(side='top')
        self.__Label25 = Label(self.__Frame68,background='Black'
            ,foreground='White',text='Output/Etc',width=50)
        self.__Label25.pack(side='top')
        self.__Frame70 = Frame(self,relief='sunken')
        self.__Frame70.pack(ipadx=10,ipady=5,padx=10,side='top')
        self.__Frame69 = Frame(self)
        self.__Frame69.pack(side='top')
        self.__Button_runit = Button(self.__Frame69,text='Run LipidWrapper!', command=self.__on_Button_runit_ButRel_1)
        self.__Button_runit.pack(side='top')
        self.__Frame7 = Frame(self.__Frame1,borderwidth=1,height=25
            ,relief='sunken',width=340)
        self.__Frame7.pack(ipadx=10,ipady=10,padx=3,side='left')
        self.__Button_load_planar = Button(self.__Frame7
            ,text='Select Planar Bilayer Model...', command=self.__on_Button_load_planar_ButRel_1)
        self.__Button_load_planar.pack(padx=5,side='left')
        self.planar_bilayer_filename = StringVar()
        self.__Label_Planar_bilayer_filename = Label(self.__Frame7
            ,text='Select...',textvariable=self.planar_bilayer_filename)
        self.__Label_Planar_bilayer_filename.pack(side='left')
        self.__Button_help1 = Button(self.__Frame7,padx='1m',text='?', command=self.__on_Button_help1_ButRel_1)
        self.__Button_help1.pack(padx=5,side='right')
        self.__Label3 = Label(self.__Frame7)
        self.__Label3.pack(padx=15,side='right')
        self.__Frame5 = Frame(self.__Frame1,borderwidth=1,height=25
            ,relief='sunken',width=340)
        self.__Frame5.pack(ipadx=10,ipady=10,padx=3,side='left')
        self.__Label9 = Label(self.__Frame5)
        self.__Label9.pack(padx=5,side='left')
        self.__Label2 = Label(self.__Frame5,text='Lipid Headgroup Marker')
        self.__Label2.pack(side='left')
        self.lipid_headgroup_marker = StringVar()
        self.__Entry_lipid_heqdgroup_marker = Entry(self.__Frame5
            ,textvariable=self.lipid_headgroup_marker,width=15)
        self.__Entry_lipid_heqdgroup_marker.pack(side='left')
        self.__Button_help2 = Button(self.__Frame5,padx='1m',text='?', command=self.__on_Button_help2_ButRel_1)
        self.__Button_help2.pack(padx=5,side='right')
        self.__Label4 = Label(self.__Frame5)
        self.__Label4.pack(padx=15,side='right')
        self.__Frame21 = Frame(self.__Frame11,borderwidth=1,height=300
            ,relief='sunken',width=340)
        self.__Frame21.pack(ipadx=10,ipady=10,padx=3,side='left')
        self.__Frame9 = Frame(self.__Frame11,borderwidth=1,height=300
            ,relief='sunken',width=340)
        self.__Frame9.pack(ipadx=10,ipady=10,padx=3,side='left')
        self.__Frame17 = Frame(self.__Frame12)
        self.__Frame17.pack(side='left')
        self.__Frame18 = Frame(self.__Frame12)
        self.__Frame18.pack(side='left')
        self.__Frame44 = Frame(self.__Frame42,borderwidth=1,height=500
            ,relief='sunken',width=340)
        self.__Frame44.pack(anchor='w',ipadx=10,ipady=10,padx=3,side='left')
        self.__Frame43 = Frame(self.__Frame42,borderwidth=1,height=500
            ,relief='sunken',width=340)
        self.__Frame43.pack(anchor='e',ipadx=10,ipady=10,padx=3,side='left')
        self.__Frame71 = Frame(self.__Frame70,borderwidth=1,height=85
            ,relief='sunken',width=340)
        self.__Frame71.pack(ipadx=10,ipady=10,side='left')
        self.__Frame72 = Frame(self.__Frame70,borderwidth=1,height=106
            ,relief='sunken',width=360)
        self.__Frame72.pack(padx=6,side='left')
        self.__Frame24 = Frame(self.__Frame21,height=15,width=340)
        self.__Frame24.pack(anchor='ne',pady=5,side='top')
        self.__Label11 = Label(self.__Frame24)
        self.__Label11.pack(padx=5,side='right')
        self.mesh_type = StringVar()
        self.__Radiobutton_mesh_from_filename = Radiobutton(self.__Frame24
            ,text='Create Mesh from a File',value='filename'
            ,variable=self.mesh_type, command=self.__on_Radiobutton_mesh_from_filename_Button_1)
        self.__Radiobutton_mesh_from_filename.pack(side='right')
        self.__Frame14 = Frame(self.__Frame21,height=25,width=340)
        self.__Frame14.pack(pady=5,side='top')
        self.__Button_load_mesh = Button(self.__Frame14
            ,text='Select Mesh File...', command=self.__on_Button_load_mesh_ButRel_1)
        self.__Button_load_mesh.pack(padx=5,side='left')
        self.mesh_filename = StringVar()
        self.__Label_mesh_filename = Label(self.__Frame14,text='Select...'
            ,textvariable=self.mesh_filename)
        self.__Label_mesh_filename.pack(side='left')
        self.__Button_help3 = Button(self.__Frame14,padx='1m',text='?', command=self.__on_Button_help3_ButRel_1)
        self.__Button_help3.pack(anchor='e',side='right')
        self.__Label12 = Label(self.__Frame14)
        self.__Label12.pack(padx=15,side='right')
        self.__Frame15 = Frame(self.__Frame21,height=25,width=340)
        self.__Frame15.pack(pady=5,side='top')
        self.__Label_max_height = Label(self.__Frame15,state='disabled'
            ,text='Maximum Height, If Image:')
        self.__Label_max_height.pack(padx=5,side='left')
        self.max_height = StringVar()
        self.__Entry_max_height = Entry(self.__Frame15,state='disabled'
            ,textvariable=self.max_height,width=5)
        self.__Entry_max_height.pack(side='left')
        self.__Button_help5 = Button(self.__Frame15,padx='1m',state='disabled'
            ,text='?', command=self.__on_Button_Help5_ButRel_1)
        self.__Button_help5.pack(side='right')
        self.__Frame22 = Frame(self.__Frame9,height=15,width=340)
        self.__Frame22.pack(anchor='nw',pady=5,side='top')
        self.__Radiobutton_mesh_from_eqn = Radiobutton(self.__Frame22
            ,text='Create Mesh from an Equation',value='equation'
            ,variable=self.mesh_type, command=self.__on_Radiobutton_mesh_from_eqn_Button_1)
        self.__Radiobutton_mesh_from_eqn.pack(anchor='nw',side='left')
        self.__Frame4 = Frame(self.__Frame9,height=25,width=340)
        self.__Frame4.pack(pady=5,side='top')
        self.__Label_equation = Label(self.__Frame4,state='disabled'
            ,text='Equation: z =')
        self.__Label_equation.pack(padx=5,side='left')
        self.the_equation = StringVar()
        self.__Entry_equation = Entry(self.__Frame4,state='disabled'
            ,textvariable=self.the_equation,width=26)
        self.__Entry_equation.pack(side='left')
        self.__Button_help4 = Button(self.__Frame4,padx='1m',state='disabled'
            ,text='?', command=self.__on_Button_help4_Button_1)
        self.__Button_help4.pack(anchor='e',side='right')
        self.__Frame16 = Frame(self.__Frame9,height=25,width=340)
        self.__Frame16.pack(pady=5,side='top')
        self.__Frame20 = Frame(self.__Frame17,height=25,width=340)
        self.__Frame20.pack(padx=10,side='top')
        self.__Frame19 = Frame(self.__Frame17,height=25,width=340)
        self.__Frame19.pack(padx=10,side='top')
        self.__Frame28 = Frame(self.__Frame17,height=25,width=340)
        self.__Frame28.pack(padx=10,side='top')
        self.__Frame26 = Frame(self.__Frame18,height=25,width=340)
        self.__Frame26.pack(padx=10,side='top')
        self.__Frame27 = Frame(self.__Frame18,height=25,width=340)
        self.__Frame27.pack(padx=10,side='top')
        self.__Frame25 = Frame(self.__Frame18,height=25,width=340)
        self.__Frame25.pack(padx=10,side='top')
        self.__Frame46 = Frame(self.__Frame44)
        self.__Frame46.pack(anchor='e',pady=5,side='top')
        self.delete_clashing_lipids = IntVar()
        
        self.__Checkbutton_delete_clashing_lipids = Checkbutton(self.__Frame46
            ,text='Delete Clashing Lipids',variable=self.delete_clashing_lipids, command=self.__on_Checkbutton_delete_clashing_lipids_ButRel_1)
        
        self.__Checkbutton_delete_clashing_lipids.pack(anchor='w',side='left')
        self.__Button_help20 = Button(self.__Frame46,padx='1m',text='?', command=self.__on_Button_help20_ButRel_1)
        self.__Button_help20.pack(anchor='e',padx=5,side='right')
        self.__Label14 = Label(self.__Frame46)
        self.__Label14.pack(padx=100,side='right')
        self.__Frame49 = Frame(self.__Frame44)
        self.__Frame49.pack(side='top')
        self.__Frame48 = Frame(self.__Frame43)
        self.__Frame48.pack(anchor='w',pady=5,side='top')
        self.fill_lipid_holes = IntVar()
        self.__Checkbutton_fill_lipid_holes = Checkbutton(self.__Frame48
            ,text='Fill Lipid Holes After Deleting'
            ,variable=self.fill_lipid_holes, command=self.__on_Checkbutton_fill_lipid_holes_ButRel_1)
        self.__Checkbutton_fill_lipid_holes.pack(anchor='w',side='left')
        self.__Label17 = Label(self.__Frame48)
        self.__Label17.pack(padx=60,side='left')
        self.__Button_help21 = Button(self.__Frame48,padx='1m',text='?', command=self.__on_Button_help21_ButRel_1)
        self.__Button_help21.pack(anchor='e',padx=5,side='right')
        self.__Frame47 = Frame(self.__Frame43)
        self.__Frame47.pack(side='top')
        self.__Frame51 = Frame(self.__Frame43)
        self.__Frame51.pack(side='top')
        self.__Frame74 = Frame(self.__Frame71,height=25,width=340)
        self.__Frame74.pack(pady=5,side='top')
        self.__Button_output_dir = Button(self.__Frame74
            ,text='Choose Output Directory...', command=self.__on_Button_output_dir_ButRel_1)
        self.__Button_output_dir.pack(padx=5,side='left')
        self.output_dir_val = StringVar()
        self.__Label_output_dir_name = Label(self.__Frame74,text='Select...'
            ,textvariable=self.output_dir_val)
        self.__Label_output_dir_name.pack(padx=5,side='left')
        self.__Button_help50 = Button(self.__Frame74,padx='1m',text='?', command=self.__on_Button_help50_ButRel_1)
        self.__Button_help50.pack(anchor='e',side='right')
        self.__Frame73 = Frame(self.__Frame71,height=25,width=340)
        self.__Frame73.pack(pady=5,side='top')
        self.generate_grid_points_and_triangles = IntVar()
        self.__Checkbutton1 = Checkbutton(self.__Frame73
            ,text=' Generate Grid Point and Triangle Files'
            ,variable=self.generate_grid_points_and_triangles)
        self.__Checkbutton1.pack(side='left')
        self.__Button_help51 = Button(self.__Frame73,padx='1m',text='?', command=self.__on_Button_help51_ButRel_1)
        self.__Button_help51.pack(side='right')
        self.__Frame75 = Frame(self.__Frame71,height=25,width=340)
        self.__Frame75.pack(pady=5,side='top')
        self.compress_output = IntVar()
        self.__Checkbutton2 = Checkbutton(self.__Frame75,text='Compress Output'
            ,variable=self.compress_output)
        self.__Checkbutton2.pack(side='left')
        self.__Button_help52 = Button(self.__Frame75,padx='1m',text='?', command=self.__on_Button_help52_ButRel_1)
        self.__Button_help52.pack(side='right')
        self.__Frame77 = Frame(self.__Frame72,height=25,width=340)
        self.__Frame77.pack(pady=5,side='top')
        self.__Label27 = Label(self.__Frame77,text='Number of Processors:')
        self.__Label27.pack(padx=5,side='left')
        self.num_proc = StringVar()
        self.__Entry1 = Entry(self.__Frame77,textvariable=self.num_proc,width=5)
        self.__Entry1.pack(side='left')
        self.__Button_help53 = Button(self.__Frame77,padx='1m',text='?', command=self.__on_Button_help53_ButRel_1)
        self.__Button_help53.pack(side='right')
        self.__Frame76 = Frame(self.__Frame72,height=25,width=340)
        self.__Frame76.pack(pady=5,side='top')
        self.use_disk_not_memory = IntVar()
        self.__Checkbutton3 = Checkbutton(self.__Frame76
            ,text='Use Disk Instead of Memory',variable=self.use_disk_not_memory)
        self.__Checkbutton3.pack(side='left')
        self.__Button_help54 = Button(self.__Frame76,padx='1m',text='?', command=self.__on_Button_help54_ButRel_1)
        self.__Button_help54.pack(side='right')
        self.__Frame78 = Frame(self.__Frame72,height=25,width=340)
        self.__Frame78.pack(side='top')
        self.__Frame23 = Frame(self.__Frame20,height=25,width=120)
        self.__Frame23.pack(side='left')
        self.__Label8 = Label(self.__Frame23,state='disabled'
            ,text='Minimum X Value:')
        self.__Label8.pack(anchor='w',side='top')
        self.__Frame29 = Frame(self.__Frame20)
        self.__Frame29.pack(side='left')
        self.min_x = StringVar()
        self.__Entry_min_x = Entry(self.__Frame29,state='disabled'
            ,textvariable=self.min_x,width=23)
        self.__Entry_min_x.pack(side='left')
        self.__Button_help6 = Button(self.__Frame29,padx='1m',text='?', command=self.__on_Button_help6_ButRel_1)
        self.__Button_help6.pack(anchor='e',padx=5,side='right')
        self.__Frame30 = Frame(self.__Frame19,height=25,width=120)
        self.__Frame30.pack(side='left')
        self.__Label10 = Label(self.__Frame30,state='disabled'
            ,text='Minimum Y Value:')
        self.__Label10.pack(anchor='w',side='left')
        self.__Frame31 = Frame(self.__Frame19)
        self.__Frame31.pack(side='left')
        self.min_y = StringVar()
        self.__Entry_min_y = Entry(self.__Frame31,state='disabled'
            ,textvariable=self.min_y,width=23)
        self.__Entry_min_y.pack(side='left')
        self.__Button_help7 = Button(self.__Frame31,padx='1m',text='?', command=self.__on_Button_help6_ButRel_1)
        self.__Button_help7.pack(padx=5,side='right')
        self.__Frame32 = Frame(self.__Frame28,height=25,width=120)
        self.__Frame32.pack(side='left')
        self.__Label13 = Label(self.__Frame32,state='disabled',text='Step X:')
        self.__Label13.pack(anchor='w',side='top')
        self.__Frame33 = Frame(self.__Frame28)
        self.__Frame33.pack(side='left')
        self.step_x = StringVar()
        self.__Entry_step_x = Entry(self.__Frame33,state='disabled'
            ,textvariable=self.step_x,width=23)
        self.__Entry_step_x.pack(side='left')
        self.__Button_help8 = Button(self.__Frame33,padx='1m',text='?', command=self.__on_Button_help6_ButRel_1)
        self.__Button_help8.pack(padx=5,side='right')
        self.__Frame35 = Frame(self.__Frame26,height=25,width=120)
        self.__Frame35.pack(side='left')
        self.__Label15 = Label(self.__Frame35,state='disabled'
            ,text='Maximum X Value:')
        self.__Label15.pack(anchor='w',side='top')
        self.__Frame34 = Frame(self.__Frame26)
        self.__Frame34.pack(side='left')
        self.max_x = StringVar()
        self.__Entry_max_x = Entry(self.__Frame34,state='disabled'
            ,textvariable=self.max_x,width=23)
        self.__Entry_max_x.pack(side='left')
        self.__Button_help9 = Button(self.__Frame34,padx='1m',text='?', command=self.__on_Button_help6_ButRel_1)
        self.__Button_help9.pack(padx=5,side='right')
        self.__Frame36 = Frame(self.__Frame27,height=25,width=120)
        self.__Frame36.pack(side='left')
        self.__Label16 = Label(self.__Frame36,state='disabled'
            ,text='Maximum Y Value:')
        self.__Label16.pack(anchor='w',side='left')
        self.__Frame37 = Frame(self.__Frame27)
        self.__Frame37.pack(side='left')
        self.max_y = StringVar()
        self.__Entry_max_y = Entry(self.__Frame37,state='disabled'
            ,textvariable=self.max_y,width=23)
        self.__Entry_max_y.pack(side='left')
        self.__Button_help10 = Button(self.__Frame37,padx='1m',text='?', command=self.__on_Button_help6_ButRel_1)
        self.__Button_help10.pack(padx=5,side='right')
        self.__Frame38 = Frame(self.__Frame25,height=25,width=120)
        self.__Frame38.pack(side='left')
        self.__Label18 = Label(self.__Frame38,state='disabled',text='Step Y:')
        self.__Label18.pack(anchor='w',side='top')
        self.__Frame39 = Frame(self.__Frame25)
        self.__Frame39.pack(side='left')
        self.step_y = StringVar()
        self.__Entry_step_y = Entry(self.__Frame39,state='disabled'
            ,textvariable=self.step_y,width=23)
        self.__Entry_step_y.pack(side='left')
        self.__Button_help11 = Button(self.__Frame39,padx='1m',text='?', command=self.__on_Button_help6_ButRel_1)
        self.__Button_help11.pack(padx=5,side='right')
        self.__Frame53 = Frame(self.__Frame49)
        self.__Frame53.pack(side='left')
        self.__Frame54 = Frame(self.__Frame49)
        self.__Frame54.pack(padx=10,side='left')
        self.__Frame64 = Frame(self.__Frame47,height=25,width=145)
        self.__Frame64.pack(side='left')
        self.__Label23 = Label(self.__Frame64,text='Filling Exhaustivenses:')
        self.__Label23.pack(anchor='w',side='top')
        self.__Frame65 = Frame(self.__Frame47)
        self.__Frame65.pack(side='left')
        self.fill_hole_exhaustiveness = StringVar()
        self.__Entry_fill_hole_exhaustiveness = Entry(self.__Frame65
            ,textvariable=self.fill_hole_exhaustiveness,width=3)
        self.__Entry_fill_hole_exhaustiveness.pack(side='top')
        self.__Frame66 = Frame(self.__Frame51,height=25,width=145)
        self.__Frame66.pack(side='left')
        self.__Label24 = Label(self.__Frame66,text='Memory Optimization:')
        self.__Label24.pack(anchor='w',side='top')
        self.__Frame67 = Frame(self.__Frame51)
        self.__Frame67.pack(side='left')
        self.memory_optimization_factor = StringVar()
        self.__Entry_memory_optimization_factor = Entry(self.__Frame67
            ,textvariable=self.memory_optimization_factor,width=3)
        self.__Entry_memory_optimization_factor.pack(side='top')
        self.__Frame55 = Frame(self.__Frame53)
        self.__Frame55.pack(side='top')
        self.__Frame56 = Frame(self.__Frame53)
        self.__Frame56.pack(side='top')
        self.__Frame58 = Frame(self.__Frame54)
        self.__Frame58.pack(side='top')
        self.__Frame57 = Frame(self.__Frame54)
        self.__Frame57.pack(side='top')
        self.__Frame45 = Frame(self.__Frame55,height=25,width=125)
        self.__Frame45.pack(anchor='w',side='left')
        self.__Label19 = Label(self.__Frame45,justify='left'
            ,text='Steric-Clash Cutoff:')
        self.__Label19.pack(anchor='w',side='top')
        self.__Frame50 = Frame(self.__Frame55)
        self.__Frame50.pack(side='left')
        self.clash_cutoff = StringVar()
        self.__Entry_clash_cutoff = Entry(self.__Frame50
            ,textvariable=self.clash_cutoff,width=3)
        self.__Entry_clash_cutoff.pack(side='top')
        self.__Frame59 = Frame(self.__Frame56,height=25,width=125)
        self.__Frame59.pack(anchor='w',side='left')
        self.__Label21 = Label(self.__Frame59,justify='left'
            ,text='Distant-Lipids Cutoff:')
        self.__Label21.pack(anchor='w',side='top')
        self.__Frame52 = Frame(self.__Frame56)
        self.__Frame52.pack(side='left')
        self.very_distant_lipids_cutoff = StringVar()
        self.__Entry_very_distant_lipids_cutoff = Entry(self.__Frame52
            ,textvariable=self.very_distant_lipids_cutoff,width=3)
        self.__Entry_very_distant_lipids_cutoff.pack(side='top')
        self.__Frame89 = Frame(self.__Frame56)
        self.__Frame89.pack(side='left')
        self.__Frame60 = Frame(self.__Frame58,height=25,width=145)
        self.__Frame60.pack(side='left')
        self.__Label20 = Label(self.__Frame60,text='Triangle-Margin Width:')
        self.__Label20.pack(anchor='w',side='top')
        self.__Frame61 = Frame(self.__Frame58)
        self.__Frame61.pack(side='left')
        self.clashing_potential_margin = StringVar()
        self.__Entry_clashing_potential_margin = Entry(self.__Frame61
            ,textvariable=self.clashing_potential_margin,width=3)
        self.__Entry_clashing_potential_margin.pack(anchor='e',side='top')
        self.__Frame62 = Frame(self.__Frame57,height=25,width=145)
        self.__Frame62.pack(side='left')
        self.__Label22 = Label(self.__Frame62,text='Triangle-Center Cutoff:')
        self.__Label22.pack(anchor='w',side='top')
        self.__Frame63 = Frame(self.__Frame57)
        self.__Frame63.pack(side='left')
        self.triangle_center_proximity_cutoff_distance = StringVar()
        self.__Entry_triangle_center_proximity_cutoff_distance = Entry(self.__Frame63
            ,textvariable=self.triangle_center_proximity_cutoff_distance,width=3)
        self.__Entry_triangle_center_proximity_cutoff_distance.pack(anchor='e'
            ,side='top')

        # some frames should just be the specified size, not expanding or shrinking
        self.__Frame5.pack_propagate(0)
        self.__Frame7.pack_propagate(0)
        self.__Frame22.pack_propagate(0)
        self.__Frame24.pack_propagate(0)
        self.__Frame14.pack_propagate(0)
        self.__Frame4.pack_propagate(0)
        self.__Frame15.pack_propagate(0)
        self.__Frame16.pack_propagate(0)
        self.__Frame12.pack_propagate(0)
        self.__Frame23.pack_propagate(0)
        self.__Frame30.pack_propagate(0)
        self.__Frame32.pack_propagate(0)
        self.__Frame35.pack_propagate(0)
        self.__Frame36.pack_propagate(0)
        self.__Frame38.pack_propagate(0)
        self.__Frame19.pack_propagate(0)
        self.__Frame20.pack_propagate(0)
        self.__Frame28.pack_propagate(0)
        self.__Frame25.pack_propagate(0)
        self.__Frame26.pack_propagate(0)
        self.__Frame27.pack_propagate(0)
        self.__Frame42.pack_propagate(0)
        self.__Frame43.pack_propagate(0)
        self.__Frame44.pack_propagate(0)
        self.__Frame45.pack_propagate(0)
        self.__Frame59.pack_propagate(0)
        self.__Frame60.pack_propagate(0)
        self.__Frame62.pack_propagate(0)
        self.__Frame64.pack_propagate(0)
        self.__Frame66.pack_propagate(0)
        self.__Frame73.pack_propagate(0)
        self.__Frame74.pack_propagate(0)
        self.__Frame75.pack_propagate(0)
        self.__Frame76.pack_propagate(0)
        self.__Frame77.pack_propagate(0)
        self.__Frame78.pack_propagate(0)
        self.__Frame71.pack_propagate(0)
        self.__Frame72.pack_propagate(0)

        # set up additional variables
        self.planar_bilayer_filename.set("Select...")
        self.mesh_filename.set("Select...")
        self.output_dir_val.set("Select...")
        self.lipid_headgroup_marker.set("_P,CHL1_O3")
        self.mesh_type.set("filename")
        self.the_equation.set("250*numpy.sin(x*x/60000 +y*y/60000)")
        self.max_height.set("15")
        self.min_x.set("0")
        self.min_y.set("0")
        self.step_x.set("25")
        self.max_x.set("100")
        self.max_y.set("100")
        self.step_y.set("25")
        self.delete_clashing_lipids.set(1)
        self.fill_lipid_holes.set(1)
        self.clash_cutoff.set("2.0")
        self.clashing_potential_margin.set("25.0")
        self.very_distant_lipids_cutoff.set("50.0")
        self.triangle_center_proximity_cutoff_distance.set("50.0")
        self.memory_optimization_factor.set("1")
        self.fill_hole_exhaustiveness.set("1")
        self.num_proc.set("1")
        self.generate_grid_points_and_triangles.set(0)
        self.compress_output.set(0)
        self.use_disk_not_memory.set(0)
        
        self.mesh_filename_reg_string = ""
        self.mesh_filename_type = ""
        self.planar_bilayer_model_filename = ""
        self.output_dir_full_path = ""

    # functions to toggle the visibility of various GUI-element groups
    def __toggle_enabled_equation(self,val):
        self.__Label_equation['state'] = val
        self.__Entry_equation['state'] = val
        self.__Button_help4['state'] = val
        
    def __toggle_enabled_mesh_filename(self,val):
        self.__Button_load_mesh['state'] = val
        self.__Label_mesh_filename['state'] = val
        self.__Button_help3['state'] = val
        
    def __toggle_enabled_mesh_max_height(self,val):
        self.__Label_max_height['state'] = val
        self.__Entry_max_height['state'] = val
        self.__Button_help5['state'] = val
        
    def __toggle_enabled_dimensions(self,val):
        self.__Label8['state'] = val
        self.__Label10['state'] = val
        self.__Label13['state'] = val
        self.__Label15['state'] = val
        self.__Label16['state'] = val
        self.__Label18['state'] = val
        self.__Entry_max_x['state'] = val
        self.__Entry_max_y['state'] = val
        self.__Entry_step_x['state'] = val
        self.__Entry_step_y['state'] = val
        self.__Entry_min_x['state'] = val
        self.__Entry_min_y['state'] = val

    def __toggle_enabled_delete_lipids(self,val):
        self.__Label19['state'] = val
        self.__Label20['state'] = val
        self.__Label21['state'] = val
        self.__Label22['state'] = val
        self.__Button_help20['state'] = val
        self.__Entry_triangle_center_proximity_cutoff_distance['state'] = val
        self.__Entry_very_distant_lipids_cutoff['state'] = val
        self.__Entry_clashing_potential_margin['state'] = val
        self.__Entry_clash_cutoff['state'] = val
        
    def __toggle_enabled_fill_holes(self,val):
        self.__Label23['state'] = val
        self.__Label24['state'] = val
        self.__Button_help21['state'] = val
        self.__Entry_fill_hole_exhaustiveness['state'] = val
        self.__Entry_memory_optimization_factor['state'] = val
        
    #Start of event handler methods
    # help buttons
    def __on_Button_Help5_ButRel_1(self,Event=None):
        tkMessageBox.showinfo('Help', 'If the bilayer model is to be generated from an image, the maximum height of the model (at those locations where the image is white) must be defined; black regions are assigned a height of 0. This feature is only available if the python PIL module has been installed on your system.')

    def __on_Button_help1_ButRel_1(self,Event=None):
        tkMessageBox.showinfo('Help', 'A PDB file containing an all-atom model of a planar lipid bilayer. LipidWrapper will wrap this lipid around the user-generated mesh.')

    def __on_Button_help20_ButRel_1(self,Event=None):
        tkMessageBox.showinfo('Help',"It's common for lipids to sterically clash at the interface of two adjacent surface-mesh tessellated triangles. LipidWrapper can delete clashing lipids.\n\nSteric-Clash Cutoff: This parameter determines how close two atoms must be (in Angstroms) to constitute a steric clash.\n\nDistant-Lipid Cutoff: Sometimes two lipids are so distant from each other that it's obvious there are no clashes. This is the distance cutoff used to determine whether or not the pair-wise distance comparison can be skipped.\n\nTriangle-Margin Width: Lipid clashes occur at the edges of adjacent tessellated triangles. Often it's faster to only check for clashes near the triangle edges. This variable specifies how far from the edges, in Angstroms, that LipidWrapper should look for clashes and (optionally) holes.\n\nTriangle-Center Cutoff: Lipid steric clashes/holes typically occur between lipids that belong to adjacent tessellated triangles. However, if tessellated triangles are small enough, clashes are possible between lipids that belong to non-adjacent triangles as well. Consequently, in addition to checking for adjacency, LipidWrapper also checks the distance between the triangle centers, using this user-specified value as a cutoff.")

    def __on_Button_help21_ButRel_1(self,Event=None):
        tkMessageBox.showinfo('Help',"Deleting lipids often leaves holes in the membrane. LipidWrapper can optionally try to fill these holes.\n\nFilling Exhaustiveness: Essentially, how long LipidWrapper should try to fill the holes. I recommend not changing this value.\n\nMemory Optimization: When the tessellated triangles are very large and consequently contain many individual lipids, the extensive pairwise distance comparisons required can result in memory errors. This parameter tells lipid Wrapper to divide the list of atoms being compared into smaller chunks. The pairwise distance comparison is performed piecewise on each chunk-chunk pair and so uses less memory, albeit at the expensive of speed. Only increase the value of this parameter if you run into memory errors. I recommend leaving it unchanged.")

    def __on_Button_help2_ButRel_1(self,Event=None):
        tkMessageBox.showinfo('Help','A comma-separated list of atom specifications (RESNAME_ATOMNAME) indicating lipid headgroups. By default, LipidWrapper looks for any atom named "P" (_P) or any cholesterol atom named "O3" (CHL1_O3).')

    def __on_Button_help3_ButRel_1(self,Event=None):
        tkMessageBox.showinfo('Help', 'If the file name ends in PDB, a surface mesh is generated from the coordinates of the PDB atoms. Otherwise, the file is assumed to be a gray-scale image, where black represents regions that are topologically low, and white represents regions that are topologically high.')

    def __on_Button_help4_Button_1(self,Event=None):
        tkMessageBox.showinfo('Help', 'Generate a surface mesh from a python-formatted equation defining z, given x and y. Python functions from the math, numpy, and scipy modules can be used.')

    def __on_Button_help50_ButRel_1(self,Event=None):
        tkMessageBox.showinfo('Help', 'If an output directory is specified, all LipidWrapper output files, as well as additional files representing the intermediate steps required to build the final bilayer, will be saved in that directory.\n\nWARNING TO WINDOWS USERS: Because of quirks in the way Windows interfaces with Tkinter, we recommend typing the path to your desired output directory, rather than relying on the Windows picker.')

    def __on_Button_help51_ButRel_1(self,Event=None):
        tkMessageBox.showinfo('Help', "Aside from producing PDB coordinates for lipid atoms, additional coordinates can be appended to the bottom of the output containing \"atoms\" named \"X\" that specify the location of the surface mesh points. Additionally, a separate file named \"triangles.tcl\" can be generated containing a tcl script that can be run in VMD to visualize the mesh surface.")

    def __on_Button_help52_ButRel_1(self,Event=None):
        tkMessageBox.showinfo('Help', 'Depending on the user options selected, LipidWrapper output can require a lot of disk space. The output can be automatically compressed using the gzip algorithm (Lempel-Ziv coding LZ77). The files can be uncompressed with the UNIX gunzip utility, or similar Windows-based packages.')

    def __on_Button_help53_ButRel_1(self,Event=None):
        tkMessageBox.showinfo('Help', 'Using multiple processors can significantly increase the speed of the LipidWrapper algorithm. Note: multi-processing is not available on Windows.')

    def __on_Button_help54_ButRel_1(self,Event=None):
        tkMessageBox.showinfo('Help', 'For very large systems, storing the growing model in memory can be problematic. In these cases, the growing model can be stored on the hard disk instead. Expect longer execution times if you choose to use this option.')

    def __on_Button_help6_ButRel_1(self,Event=None):
        tkMessageBox.showinfo('Help', 'The maximum and minimum values of X and Y must be defined if you\'re creating a mesh from an equation or an image. The distance between mesh points in both the X and Y directions (Step X and Step Y) must also be defined.')

    # load-file buttons
    def fix_filename(self, filename):
        return filename.replace('\\','/')

    def __on_Button_load_mesh_ButRel_1(self,Event=None):
         afile = tkFileDialog.askopenfilename(parent=self.master,title='Choose a File')
         if afile != "" and os.path.exists(afile):
             afile = self.fix_filename(afile)
             self.mesh_filename.set(os.path.basename(afile))
             self.mesh_filename_reg_string = afile
             if "PDB" in afile.upper(): self.mesh_filename_type = "PDB"
             elif "DAE" in afile.upper(): self.mesh_filename_type = "DAE"
             else: self.mesh_filename_type = "IMG"
             self.__on_Radiobutton_mesh_from_filename_Button_1()

    def __on_Button_load_planar_ButRel_1(self,Event=None):
         afile = tkFileDialog.askopenfilename(parent=self.master,title='Choose a Directory')
         if afile != "" and os.path.exists(afile):
             afile = self.fix_filename(afile)
             self.planar_bilayer_filename.set(os.path.basename(afile))
             self.planar_bilayer_model_filename = afile
             
    def __on_Button_output_dir_ButRel_1(self,Event=None):
         adir = tkFileDialog.askdirectory(parent=self.master,title='Choose a Directory')
         if adir != "": # and os.path.exists(afile):
             adir = self.fix_filename(adir)
             self.output_dir_val.set('...' + os.sep + os.path.basename(adir))
             self.output_dir_full_path = adir

    # run the program button
    def __on_Button_runit_ButRel_1(self,Event=None):
        # check a few things
        
        if self.planar_bilayer_model_filename == "":
            tkMessageBox.showinfo('Error',"You need to specify a planar bilayer model before proceeding.")
            return
        
        if self.__Label_mesh_filename['state'] == "normal" and self.mesh_filename_reg_string == "":
            tkMessageBox.showinfo('Error',"If you with to create your mesh from a file rather than an equation, you need to specify the mesh filename. Choose a PNG, PDB, or collada DAE file.")
            return
        
        if self.output_dir_full_path == "":
            tkMessageBox.showinfo('Error',"You need to specify an output directory before proceeding. Note that the command-line version of LipidWrapper does not require an output directory. It can simply write the final PDB coordinates to the screen.")
            return
        
        # run LipidWrapper!
        argv = [] # construct a fake commandline-arguments list
        
        # first, ones that will always be specified
        argv.append('lipidwrapper.py')
        argv.append('--lipid_pdb_filename')
        argv.append(self.planar_bilayer_model_filename)
        argv.append('--lipid_headgroup_marker')
        argv.append(self.lipid_headgroup_marker.get().replace(" ",''))
        argv.append('--delete_clashing_lipids')
        if self.delete_clashing_lipids.get() == 1: argv.append("TRUE")
        else: argv.append("FALSE")
        
        argv.append('--number_of_processors')
        argv.append(self.num_proc.get())
        
        argv.append('--show_grid_points')
        if self.generate_grid_points_and_triangles.get() == 1: argv.append("TRUE")
        else: argv.append("FALSE")
               
        argv.append('--create_triangle_tcl_file')
        if self.generate_grid_points_and_triangles.get() == 1: argv.append("TRUE")
        else: argv.append("FALSE")
        
        argv.append('--use_disk_instead_of_memory')
        if self.use_disk_not_memory.get() == 1: argv.append("TRUE")
        else: argv.append("FALSE")

        argv.append('--compress_output')
        if self.compress_output.get() == 1: argv.append("TRUE")
        else: argv.append("FALSE")
        
        # though optional for the commandline program, the below will be required in the GUI
        argv.append('--output_directory')
        argv.append(self.output_dir_full_path)

        # now ones that are optional
        if self.__Entry_max_height['state'] == "normal": 
            argv.append('--max_height')
            argv.append(self.max_height.get())
             
        if self.__Entry_equation['state'] == "normal": 
            argv.append('--surface_equation')
            argv.append("z = " + self.the_equation.get())
             
        if self.__Entry_min_x['state'] == "normal": 
            argv.append('--min_x')
            argv.append(self.min_x.get())
             
        if self.__Entry_max_x['state'] == "normal": 
            argv.append('--max_x')
            argv.append(self.max_x.get())
             
        if self.__Entry_min_y['state'] == "normal": 
            argv.append('--min_y')
            argv.append(self.min_y.get())
             
        if self.__Entry_max_y['state'] == "normal": 
            argv.append('--max_y')
            argv.append(self.max_y.get())
             
        if self.__Entry_step_x['state'] == "normal": 
            argv.append('--step_x')
            argv.append(self.step_x.get())
             
        if self.__Entry_step_y['state'] == "normal": 
            argv.append('--step_y')
            argv.append(self.step_y.get())
             
        if self.__Entry_fill_hole_exhaustiveness['state'] == "normal": 
            argv.append('--fill_hole_exhaustiveness')
            argv.append(self.fill_hole_exhaustiveness.get())

        if self.__Entry_memory_optimization_factor['state'] == "normal": 
            argv.append('--memory_optimization_factor')
            argv.append(self.memory_optimization_factor.get())
             
        if self.__Entry_clash_cutoff['state'] == "normal": 
            argv.append('--clash_cutoff')
            argv.append(self.clash_cutoff.get())
             
        if self.__Entry_very_distant_lipids_cutoff['state'] == "normal": 
            argv.append('--very_distant_lipids_cutoff')
            argv.append(self.very_distant_lipids_cutoff.get())

        if self.__Entry_clashing_potential_margin['state'] == "normal": 
            argv.append('--clashing_potential_margin')
            argv.append(self.clashing_potential_margin.get())

        if self.__Entry_triangle_center_proximity_cutoff_distance['state'] == "normal": 
            argv.append('--triangle_center_proximity_cutoff_distance')
            argv.append(self.triangle_center_proximity_cutoff_distance.get())

        if self.__Checkbutton_fill_lipid_holes['state'] == "normal": 
            argv.append('--fill_holes')
            if self.fill_lipid_holes.get() == 1: argv.append("TRUE")
            else: argv.append("FALSE")

        if self.mesh_filename_reg_string != "" and self.__Label_mesh_filename['state'] == "normal": 
            argv.append('--surface_filename')
            argv.append(self.mesh_filename_reg_string)

        #print argv
        
        tkMessageBox.showinfo('Starting...',"The lipid-bilayer model will now be generated. This could take some time depending on the options you selected. Check the console/terminal window to monitor the progress.")

        import lipidwrapper
        lipidwrapper.run_program(argv)
        
        tkMessageBox.showinfo('Done!',"LipidWrapper has finished running. The generated bilayer model is located in " + self.output_dir_full_path)
        
        sys.exit(0)
    
    # code to control other GUI element
    def __on_Checkbutton_delete_clashing_lipids_ButRel_1(self,Event=None):
        if self.delete_clashing_lipids.get() == 0: 
            self.__toggle_enabled_delete_lipids('disabled')
            self.__toggle_enabled_fill_holes('disabled')
            self.__Checkbutton_fill_lipid_holes['state'] = 'disabled'
        else: 
            self.__toggle_enabled_delete_lipids('normal')
            self.__on_Checkbutton_fill_lipid_holes_ButRel_1()
            self.__Checkbutton_fill_lipid_holes['state'] = 'normal'

    def __on_Checkbutton_fill_lipid_holes_ButRel_1(self,Event=None):
        if self.fill_lipid_holes.get() == 0: self.__toggle_enabled_fill_holes('disabled')
        else: self.__toggle_enabled_fill_holes('normal')

    def __on_Radiobutton_mesh_from_eqn_Button_1(self,Event=None):
        self.__toggle_enabled_equation('normal')
        self.__toggle_enabled_mesh_filename('disabled')
        self.__toggle_enabled_mesh_max_height('disabled')
        self.__toggle_enabled_dimensions('normal')

    def __on_Radiobutton_mesh_from_filename_Button_1(self,Event=None):
        self.__toggle_enabled_equation('disabled')
        self.__toggle_enabled_mesh_filename('normal')
        
        if self.mesh_filename_type == "IMG": self.__toggle_enabled_mesh_max_height("normal")
        else: self.__toggle_enabled_mesh_max_height("disabled")
        
        if self.mesh_filename_type == "IMG": self.__toggle_enabled_dimensions("normal")
        else: self.__toggle_enabled_dimensions("disabled")

if __name__ == '__main__':

    Root = Tk()
    App = LipidWrapper(Root)
    App.pack(expand='yes',fill='both')

    Root.geometry('750x693+10+10')
    Root.resizable(0,0)
    Root.title('LipidWrapper')
    Root.mainloop()
