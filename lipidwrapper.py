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

version = "1.14"

import sys
import math
import numpy
import scipy
from scipy.spatial import Delaunay
from scipy.spatial.distance import cdist
from scipy.spatial.distance import pdist
import getopt
import textwrap
import time
import random
import gc
import multiprocessing
import os
import shutil
import cPickle as pickle
import string
import platform
import gzip

################## MULTIPROCESSING ##################

class multi_threading():
    """A class for running calculations on multiple processors"""
    
    results = []
    
    def __init__(self, inputs, num_processors, task_class, params, progress_bar_prefix=''):
        """Launches a calculation on multiple processors
            
        Arguments:
        inputs -- A list, containing all the input required for the calculation
        num_processors -- An integer, the requested number of processors to use
        task_class -- An class, the class governing what calculations will be run on a given thread
        progress_bar_prefix -- An optional string, the prefix to append to the progress bar
        
        Returns:
        Nothing, though the objects self.results list is populated with the calculation results
        
        """
        
        # set up the progress bar
        self.results = []
        indices_to_star = []
        if len(inputs) < 50:
            indices_to_star = range(len(inputs))
        else:
            while len(indices_to_star) < 50:
                indx_to_add = random.choice(range(len(inputs)))
                if not indx_to_add in indices_to_star: indices_to_star.append(indx_to_add)

        if progress_bar_prefix != '':
            toadd = 78 - len(progress_bar_prefix) - 50
            progress_bar_prefix = progress_bar_prefix + (" "*toadd)
            sys.stdout.write(progress_bar_prefix)
        
        if num_processors == 1: # so just running on 1 processor, perhaps under windows
            single_thread = task_class()
            single_thread.total_num_tasks = len(inputs)
            single_thread.indices_to_star = indices_to_star

            single_thread.results = []
            for item in inputs: single_thread.value_func(item, None)
            
            self.results = single_thread.results
            
        else: # so it actually is running on multiple processors

            cpu_count = 1
            cpu_count = multiprocessing.cpu_count()
    
            # first, if num_processors <= 0, determine the number of processors to use programatically
            if num_processors <= 0: num_processors = cpu_count
    
            # reduce the number of processors if too many have been specified
            if len(inputs) < num_processors: num_processors = len(inputs)
            
            if len(inputs) == 0: # if there are no inputs, there's nothing to do.
                self.results = []
                return
            
            # now, divide the inputs into the appropriate number of processors
            inputs_divided = {}
            for t in range(num_processors): inputs_divided[t] = []
    
            for t in range(0, len(inputs), num_processors):
                for t2 in range(num_processors):
                    index = t + t2
                    if index < len(inputs): inputs_divided[t2].append(inputs[index])

            # now, run each division on its own processor
            running = multiprocessing.Value('i', num_processors)
            mutex = multiprocessing.Lock()
    
            arrays = []
            threads = []
            for i in range(num_processors):
                athread = task_class()
                athread.total_num_tasks = len(inputs)
                athread.indices_to_star = indices_to_star
                
                threads.append(athread)
                arrays.append(multiprocessing.Array('i',[0, 1]))
    
            results_queue = multiprocessing.Queue() # to keep track of the results
            
            processes = []
            for i in range(num_processors):
                p = multiprocessing.Process(target=threads[i].runit, args=(running, mutex, results_queue, inputs_divided[i]))
                p.start()
                processes.append(p)
    
            while running.value > 0: is_running = 0 # wait for everything to finish

            # compile all results into one list
            for thread in threads:
                chunk = results_queue.get()
                self.results.extend(chunk)
                
        if progress_bar_prefix != '': print # because the progress bar is now done

class general_task:
    """A parent class of others that governs what calculations are run on each thread"""
    
    results = []
    
    def print_star_if_appropriate(self, current_index):
        """If appropriate, prints an asterix as part of the progress bar
            
        Arguments:
        current_index -- An integer, the index of the current calculation
        
        """

        # if the current index is one of the ones that you should star, write out the asterix
        if current_index in self.indices_to_star: sys.stdout.write("*")
    
    def runit(self, running, mutex, results_queue, items):
        """Launches the calculations on this thread
            
        Arguments:
        running -- A multiprocessing.Value object
        mutex -- A multiprocessing.Lock object
        results_queue -- A multiprocessing.Queue() object for storing the calculation output
        items -- A list, the input data required for the calculation
        
        """

        for item in items: self.value_func(item, results_queue)
        
        mutex.acquire()
        running.value -= 1
        mutex.release()
        results_queue.put(self.results)

################## FUNCTIONS AND CLASSES TO EXTEND NUMPY ##################

class Quaternion:
    """A class supporting quaternion arithmetic"""
    
    def __init__(self, s, x, y, z):
        self.v = numpy.empty(4)
        self.v[0] = s
        self.v[1] = x
        self.v[2] = y
        self.v[3] = z
    
    def __str__(self):
        """String containing quaternion information in the form of s x y z
            
        Returns:
        A string, containing all information about this quaternion
            
        """
        
        return "" + str(self.v[0]) + "\t" + str(self.v[1]) + "\t" + str(self.v[2]) + "\t" + str(self.v[3])
    
    def copy_of(self):
        """Returns a copy of self"""
        
        return Quaternion(self.v[0], self.v[1], self.v[2], self.v[3])
    
    def load_from_mat(self, m):
        """Converts a rotation matrix that is pure orthogonal (det(matrix)=1) into a Quaternion. Adapted from http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm 
                
        Arguments:
        m -- A 2D numpy.array representing a pure orthogonal matrix
                
        """

        #Make sure m is a 3x3 array
        if m.shape[0] != 3 or m.shape[1] != 3:
            print "Could not load quaternion from matrix...size is not (3x3)"
            return
        
        #Check that matrix is orthogonal. m_T = m_inv
        if not numpy.array_equal(numpy.transpose(m),numpy.linalg.inv(m)):
            print "Load Quaternion error. Matrix is not orthogonal"
            return
        
        #Need to make sure that the matrix is special orthogonal
        if math.fabs(1-numpy.linalg.det(m)) > 0.000001: #Done for rounding errors
            print "Load Quaternion error.  Determinant is not 1"
            return
        
        #First calculate the sum of the diagonal elements
        t = m.trace()
        
        if t > 0:
            S = math.sqrt(t + 1.0) * 2
            self.v[0] = .25 * S
            self.v[1] = (m[2,1] - m[1,2]) / S
            self.v[2] = (m[0,2] - m[2,0]) / S
            self.v[3] = (m[1,0] - m[0,1]) / S
        elif m[0,0] > m[1,1] and m[0,0] > m[2,2]:
            S = math.sqrt(1.0 + m[0,0] - m[1,1] - m[2,2]) * 2
            self.v[0] = (m[2,1] - m[1,2]) / S
            self.v[1] = .25 * S
            self.v[2] = (m[0,1] + m[1,0]) / S
            self.v[3] = (m[0,2] + m[2,0]) / S
        elif m[1,1] > m[2,2]:
            S = math.sqrt(1.0 + m[1,1] - m[0,0] - m[2,2]) * 2
            self.v[0] = (m[0,2] - m[2,0]) / S
            self.v[1] = (m[0,1] + m[1,0]) / S
            self.v[2] = .25 * S
            self.v[3] = (m[2,1] + m[1,2]) / S
        else:
            S = math.sqrt(1.0) * 2
            self.v[0] = (m[1,0] - m[0,1]) / S
            self.v[1] = (m[0,2] + m[2,0]) / S
            self.v[2] = (m[2,1] + m[1,2]) / S
            self.v[3] = .25 * S
    
    def rep_as_44_matrix(self):
        """Creates a 4x4 matrix representation of the Quaternion.
            
        Returns:
        A 4x4 numpy array
        
        """
        
        n = self.normalize()
        qw = n.v[0]
        qx = n.v[1]
        qy = n.v[2]
        qz = n.v[3]
        
        return numpy.array([[qw,qx,qy,qz],[-qx,qw,-qz,qy],[-qy,qz,qw,-qx],[-qz,-qy,qx,qw]])

    def to_matrix(self):
        """Converts to a normalized 3x3 matrix
            
        Returns:
        A 3x3 numpy.array, corresponding to the quaternion
        
        """
        
        #First normalize
        n = self.normalize()
        qw = n.v[0]
        qx = n.v[1]
        qy = n.v[2]
        qz = n.v[3]
        return numpy.array(
                           [[1.0 - 2.0 * qy * qy - 2.0 * qz * qz, 2.0 * qx * qy - 2.0 * qz * qw, 2.0 * qx * qz + 2.0 * qy * qw],
                            [2.0 * qx * qy + 2.0 * qz * qw, 1.0 - 2.0 * qx * qx - 2.0 * qz * qz, 2.0 * qy * qz - 2.0 * qx * qw],
                            [2.0 * qx * qz - 2.0 * qy * qw,2.0 * qy * qz + 2.0 * qx * qw, 1.0 - 2.0 * qy * qy - 2.0 * qx * qx]])
                
    def add(self, q2):
        """Adds two quaternions
            
        Arguments:
        q2 -- A quaternion, to be added to self
        
        Returns:
        A Quaternion, with the values corresponding to self + q2
        
        """
        
        return Quaternion(self.v[0] + q2.v[0], self.v[1] + q2.v[1], self.v[2] + q2.v[2], self.v[3] + q2.v[3])
    
    def invert(self):
        """Takes the inverse of the quaternion for "division"
        
        Returns:
        A Quaternion, with the values corresponding to self^-1
        
        """
        
        return Quaternion(self.v[0], -1 * self.v[1], -1 * self.v[2], -1 * self.v[3])

    def minus(self, q2):
        """Multiplies two quaternions
        
        Arguments:
        q2 -- A quaternion, to be subtracted from self
        
        Returns:
        A Quaternion, with the values corresponding to self - q2
        
        """
        
        return Quaternion(self.v[0] - q2.v[0], self.v[1] - q2.v[1], self.v[2] - q2.v[2], self.v[3] - q2.v[3])
    
    def multiply(self, q2):
        """Multiplies two quaternions
            
        Arguments:
        q2 -- A quaternion, to be multiplied with self
        
        Returns:
        A Quaternion, with the values corresponding to self * q2
        
        """
        
        return Quaternion(self.v[0] * q2.v[0] - self.v[1] * q2.v[1] - self.v[2] * q2.v[2] - self.v[3] * q2.v[3],
                          self.v[1] * q2.v[0] + self.v[0] * q2.v[1] + self.v[2] * q2.v[3] - self.v[3] * q2.v[2],
                          self.v[0] * q2.v[2] - self.v[1] * q2.v[3] + self.v[2] * q2.v[0] + self.v[3] * q2.v[1],
                          self.v[0] * q2.v[3] + self.v[1] * q2.v[2] - self.v[2] * q2.v[1] + self.v[3] * q2.v[0])
    
    def normalize(self):
        """Normalizes the quaternion
            
        Returns:
        A normalized Quaternion
        
        """
        
        #First normalize
        n = math.sqrt(math.pow(self.v[0],2) + math.pow(self.v[1],2) + math.pow(self.v[2],2) + math.pow(self.v[3],2))
        
        return Quaternion(self.v[0] / n, self.v[1] / n, self.v[2] / n, self.v[3] / n)
                
    def scale(self, scalar):
        """Scales a quaternion
                
        Arguments:
        scalar -- the value to scale the quaternion by
        
        Returns:
        A Quaternion, with the values corresponding to self * scalar
            
        """
    
        return Quaternion(self.v[0] * scalar, self.v[1] * scalar, self.v[2] * scalar, self.v[3] * scalar)

def get_numpy_slice(numpy_array, indices):
    """A replacement for numpy's numpy_array[indices] functionality, where numpy_array and indices are both matrices.
    Unlike numpy's default behavior, this function returns an empty array if b is empty, rather than throwing an error.
        
    Arguments:
    numpy_array -- A numpy array
    indices -- A numpy array, the indices of the elements of numpy_array to return.
    
    Returns:
    A numpy array, numpy_array[indices] if indices is not empty, an empty numpy array otherwise.
    
    """

    try: return numpy_array[indices]
    except:
        if len(indices) == 0: return numpy.array([]) # why in the world isn't this numpy's default behavior?
        else: print "Error!"

################## MODIFYING AND MANIPULATING MOLECULAR MODELS ##################

class Molecule:
    """Loads, saves, and manupulates molecuar models."""
    
    def __init__ (self):
        """Initializes the variables of the Molecule class."""
        
        self.in_triangle_margin = True
        self.in_triangle_submargin = False
        self.headgroup_index = None
    
    def get_headgroup_index(self, lipid_headgroup_marker):
        """Get's the indices of the current molecule's headgroup
            
        Arguments:
        lipid_headgroup_marker -- A tuple of the form (chain, resname, resid, atomname) specifying the headgroup
        
        Returns:
        An integer, the index of this molecule's headgroup
        
        """

        if self.headgroup_index == None: self.headgroup_index = self.get_indices_of_mask_match(lipid_headgroup_marker)[0] # so calculate it only if it's never been calculated before
        return self.headgroup_index
    
    def load_pdb(self, filename):
        """Loads a PDB file into the current Molecule object from a file
        
        Arguments:
        filename -- a string, the name of the file to load

        """

        # Now load the file into a list
        file = open(filename,"r")
        lines = file.readlines()
        file.close()
        
        # load the molecule from the list
        self.load_pdb_from_lines(lines)
        
    def load_pdb_from_lines(self, lines):
        """Loads a PDB file into the current Molecule object from a list of PDB lines
        
        Arguments:
        lines -- a list, containing the PDB lines to load into the current object

        """

        self.__init__()

        gc.disable() # because appending objects slows down code if garbage collection turned on
        
        # set up the numpy arrays to store the data
        self.atom_inf_string_vals = numpy.empty((len(lines), 4), dtype='|S9') # chain, resname, atomname, id_keys 
        self.atom_inf_resids = numpy.empty(len(lines), dtype='|S4')
        self.all_atoms_numpy = numpy.empty((len(lines), 3))
        
        # read in the data from the lines
        count = 0
        for t in range(0,len(lines)):
            line=lines[t]
            if len(line) >= 7:
                if line[0:4]=="ATOM" or line[0:6]=="HETATM": # Load atom data (coordinates, etc.)
                    count = count + 1
                    
                    self.all_atoms_numpy[t][0] = float(line[30:38])
                    self.all_atoms_numpy[t][1] = float(line[38:46])
                    self.all_atoms_numpy[t][2] = float(line[46:54])
                    
                    resname = line[16:21].strip()
                    atomname = line[11:16].strip()

                    try: resid = line[22:26].strip()
                    except: resid = "0"

                    self.atom_inf_string_vals[t][0] = line[21:22].strip() # chain
                    self.atom_inf_string_vals[t][1] = resname # resname
                    self.atom_inf_string_vals[t][2] = atomname # atomname
                    self.atom_inf_string_vals[t][3] = resname + "_" + atomname # id_keys 
                
                    self.atom_inf_resids[t] = resid

        gc.enable()

        # now resize the array, cutting out bottom parts that were never populated
        self.atom_inf_string_vals = self.atom_inf_string_vals[:count]
        self.atom_inf_resids = self.atom_inf_resids[:count]
        self.all_atoms_numpy = self.all_atoms_numpy[:count]

    def save_pdb(self, filename):
        """Saves data to a PDB file.
        
        Arguments:
        filename -- A string, the filename to be written.
        
        """
        
        toprint = ""

        file = open(filename,"w")
        for index in range(len(self.all_atoms_numpy)): file.write(self.create_pdb_line(index) + "\n")
        file.close()
        
    def set_undo_point(self): # you can restore all atom positions to some undo point. This sets that point.
        """Sets ("saves") the undo point of all atoms. Any subsequent manipulations of atomic coordinates can be "undone" by reseting to this configuration."""
        
        self.all_atoms_numpy_undo = numpy.copy(self.all_atoms_numpy)

    def undo(self):
        """Resets the coordinates of all atoms to those saved using the set_undo_point function."""
        
        self.all_atoms_numpy = numpy.copy(self.all_atoms_numpy_undo)
        
    def rotate_mol_quat(self, rot_quat):
        """Support function that rotates a molecule according to a rotation quaternion
            
        Arguments:
        mol -- A molecule to be rotated in matrix format
        rot_quat -- Quaternion to rotate molecule
        
        """

        rot_mat = rot_quat.to_matrix()
        self.all_atoms_numpy = numpy.dot(self.all_atoms_numpy,rot_mat)
            
    def baseN(self, num, b, numerals="0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"):
        """Return the value of a number in another base
            
        Arguments:
        num -- An integer, the number in base 10.
        b -- An integer, the number of the new base
        numerals -- An optional string, containing the numerals to use in the new base
        
        Returns:
        A string, the representation of the original integer, now in the specified base
        
        """

        return ((num == 0) and numerals[0]) or (self.baseN(num // b, b, numerals).lstrip(numerals[0]) + numerals[num % b])
    
    def create_pdb_line(self, index, output_index=None, output_resid=None):
        """Create a string formatted according to the PDB standard from the atomic information contained in this atom class.
        
        Arguments:
        index -- An integer, the index of the atom in the Molecule object.
        output_index -- An optional integer, the index to use in the PDB-line output. If not specified, index is used.
        output_resid -- An optional integer, the resid to use in the PDB-line output. If not specified, the existing resid is used.
        
        Returns:
        A string, formatted according to the PDB standard.
        
        """
        
        # use user-specified index if provided
        if output_index is None: output_index = str(index)
        else: output_index = str(output_index)
        
        # PDB format is fixed column, so if the index is too big just turn it into stars
        if len(output_index) >= 7: output_index = "******"
        
        # use the user-specified resid if provided
        if output_resid is None: output_resid = self.atom_inf_resids[index]
        else: output_resid = str(output_resid)
        
        # PDB format is fixed column, so if the resid is too big, switch over to a string identifier that is unique to each residue
        if len(output_resid) >= 5: # you need to start using letters in the resid after 9999
            # 2383280 is "a001" in base 62.
            # so map 10000 to 2383280 and convert to base 62.
            output_resid = self.baseN(int(output_resid) + 2373280, 62)
            # max using this method is 35999 residues
        
        # create the PDB line
        output = "ATOM "
        output = output + str(output_index).rjust(6) + self.atom_inf_string_vals[index][2].rjust(5) + self.atom_inf_string_vals[index][1].rjust(5) + self.atom_inf_string_vals[index][0].rjust(1) + output_resid.rjust(4)
        output = output + ("%.3f" % self.all_atoms_numpy[index][0]).rjust(12)
        output = output + ("%.3f" % self.all_atoms_numpy[index][1]).rjust(8)
        output = output + ("%.3f" % self.all_atoms_numpy[index][2]).rjust(8)
        
        return output

    def copy_of(self):
        """Create a copy of the current molecule
        
        Returns:
        A Molecule object, a copy of the current one
        
        """
        
        new = Molecule()

        new.atom_inf_string_vals=self.atom_inf_string_vals.copy()
        new.atom_inf_resids=self.atom_inf_resids.copy()
        new.all_atoms_numpy = self.all_atoms_numpy.copy()
        new.headgroup_index = self.headgroup_index
        
        return new

    def portion_of(self, list_of_indices):
        """Get a portion of the current molecule
        
        Arguments:
        list_of_indices -- A numpy array of integers, corresponding to the atoms in the current Molecule that should be in the return portion
        
        Returns:
        A Molecule object, containing only the atoms specified in list_of_indices
        
        """

        new = Molecule()
        
        new.atom_inf_string_vals=self.atom_inf_string_vals[list_of_indices]
        new.atom_inf_resids=self.atom_inf_resids[list_of_indices]
        new.all_atoms_numpy=self.all_atoms_numpy[list_of_indices]
        
        return new

    def get_indices_of_mask_match(self, masks):
        """Get the indices of atoms that match queries (masks)
        
        Arguments:
        masks -- A list of tuples. Each tuple is a mask/query of the format (chain, resname, resid, atomname). '', None, or -9999 all mean wildcard. For example: [('A', 'CHL1', 57, 'O3'), ('B', '', 783, 'P')]
        
        Returns:
        A numpy array of integers, containing the indices of the atoms that match the mask
        
        """

        indices = numpy.array([], dtype=int)
        
        for mask in masks:
            chain = mask[0]
            resname = mask[1]
            resid = mask[2]
            atomname = mask[3]
            
            # get all the indices of the ones that have the same resname
            if chain == "" or chain is None or chain == -9999: indices_of_ones_with_this_chain = numpy.array(range(len(self.atom_inf_string_vals))) # so it can be anything
            else: indices_of_ones_with_this_chain = numpy.nonzero(self.atom_inf_string_vals[:,0] == chain)[0]

            if resname == "" or resname is None or resname == -9999: indices_of_ones_with_this_resname = numpy.array(range(len(self.atom_inf_string_vals))) # so it can be anything
            else: indices_of_ones_with_this_resname = numpy.nonzero(self.atom_inf_string_vals[:,1] == resname)[0]
            
            if resid == "" or resid is None or resid == -9999 or resid == "-9999": indices_of_ones_with_this_resid = numpy.array(range(len(self.atom_inf_resids))) # so it can be anything
            else: indices_of_ones_with_this_resid = numpy.nonzero(self.atom_inf_resids == resid)[0]

            if atomname == "" or atomname is None or atomname == -9999: indices_of_ones_with_this_atomname = numpy.array(range(len(self.atom_inf_string_vals))) # so it can be anything
            else: indices_of_ones_with_this_atomname = numpy.nonzero(self.atom_inf_string_vals[:,2] == atomname)[0]
            
            # the intersection is the one that has both
            indices_in_all = numpy.intersect1d(indices_of_ones_with_this_chain, indices_of_ones_with_this_resname, assume_unique=True)
            indices_in_all = numpy.intersect1d(indices_in_all, indices_of_ones_with_this_atomname, assume_unique=True)
            indices_in_all = numpy.intersect1d(indices_in_all, indices_of_ones_with_this_resid, assume_unique=True)
            indices = numpy.union1d(indices, indices_in_all)
            
        return indices

class Triangle():
    """A class describing a triangle in three dimensions"""
    
    def __init__(self, pts):
        """Create a Triangle object.
        
        Arguments:
        pts -- A 3x3 numpy array containing the points of the triangle
        
        """

        self.points = pts

    def __getitem__(self, index):
        """Get one of the triangel points
        
        Arguments:
        index -- An integer, the index of the point
        
        Returns:
        A (1x3) numpy array, the coordinates of the requested point
        
        """

        return self.points[index]

    def center(self):
        """Get the triangle center
        
        Returns:
        A (1x3) numpy array, the coordinates of the triangle center
        
        """

        try: return self.center_pt
        except:
            self.center_pt = numpy.average(self.points, 0)
            return self.center_pt
        
    def radii(self):
        """Get the set of distances from the trangle center to each of its points (i.e., the "radii" of the triangle)
        
        Returns:
        A 1x3 numpy array, containing the distances/radii
        
        """

        try: return self.radii_lengths
        except:
            self.radii_lengths = cdist(numpy.array([self.center()]), self.points, 'euclidean')
            return self.radii_lengths

    def max_radius(self):
        """Returns a float, the maximum distance from the triangle center to any of the contituent points"""

        try: return self.radius_length
        except:
            self.radius_length = numpy.amax(self.radii())
            return self.radius_length
        
    def project_points_onto_triangle(self, pts):
        """Projects a series of points onto the triangle plane
        
        Arguments:
        pts -- An nx3 numpy array, containing the points to be projected
        
        Returns:
        An nx3 numpy array, the coordinates of the requested point now projected onto the triangle plane
        
        """

        # define the triangle plane
        AB = self.points[1]-self.points[0]
        AC = self.points[2]-self.points[0]
        normal_to_plane = numpy.cross(AB, AC)
        normal_to_plane = normal_to_plane/numpy.linalg.norm(normal_to_plane) # so it's a unit vector now
        
        # project the headgroups onto the triangle plane
        return pts - ((numpy.transpose([numpy.dot(pts - self.points[1], normal_to_plane)])) * normal_to_plane)

    def get_indices_of_points_within_triangle_boundaries(self, pts):
        """For a set of points that are coplanar with the current triangle, identify the indices of the points that are within the triangle boundaries
        
        Arguments:
        pts -- An nx3 numpy array, containing the points to be evaluated
        
        Returns:
        A numpy array, containing the indices of the points that fall within the triangle boundaries
        
        """
        
        # if the triangle is really just a line or a point, return an empty list
        if numpy.array_equal(self.points[0], self.points[1]) or numpy.array_equal(self.points[0], self.points[2]) or numpy.array_equal(self.points[1], self.points[2]): return numpy.array([])
        
        # some error has occurred previously, return an empty list
        if numpy.isnan(self.points).any(): return numpy.array([])
        
        # get bounding box
        tri_min = numpy.min(self.points,0) - numpy.array([1e-6, 1e-6, 1e-6]) # subtract a little to avoid rounding errors
        tri_max = numpy.max(self.points,0) + numpy.array([1e-6, 1e-6, 1e-6])
        
        # identify points that couldn't possibly be in the triangle because they're outside the box
        pts_not_in_triangle = (pts < tri_min).any(1)
        pts_not_in_triangle = numpy.logical_or(pts_not_in_triangle, (pts > tri_max).any(1))
        pts_potentially_in_triangle = numpy.logical_not(pts_not_in_triangle)
        
        # get the indices of the ones that could possibly be inside the triangle
        indices_of_pts_potentially_in_triangle = numpy.nonzero(pts_potentially_in_triangle)[0]
        
        # verify which ones really are in the triangle
        indices_to_keep = []
        for t in indices_of_pts_potentially_in_triangle:
            
            # calculate three vectors from the triangle verticies to the projection
            t_v1 = self.points[0] - pts[t]
            t_v2 = self.points[1] - pts[t]
            t_v3 = self.points[2] - pts[t]
            
            # get the appropriate angles
            angle1 = angle_between(t_v1, t_v2)
            angle2 = angle_between(t_v1, t_v3)
            angle3 = angle_between(t_v2, t_v3)
            
            # sometimes, if a triangle is small and the comparison point is very far away,
            # two of the vectors can end up being the same, especially after normalization.
            # Inevitably, in this case the point is not in the triangle.
            # we should account for that.
            if angle1 == "NORMALIZED VECTORS EQUAL!" or angle2 == "NORMALIZED VECTORS EQUAL!" or angle3 == "NORMALIZED VECTORS EQUAL!": continue
                
            if math.fabs(angle1 + angle2 + angle3 - 2 * math.pi) < 0.01: # it's inside the triangle
                indices_to_keep.append(t)

        return numpy.array(indices_to_keep)
    
    def new_triangle_expanded_by_margin(self, margin):
        """Return a triangle that is larger or smaller than the current one. The triangle is not simply scaled. A new triangle is carefully constructed such that each of its edges and the edges of the original triangle are a user-specified distance apart.
        
        Arguments:
        margin -- A float, the distance between the edges of the new triangle and the corresponding edges of the original triangle
        
        Returns:
        A Triangle object, the new triangle
        
        """
        
        # first, if this triangle is already a single point and the margin is negative, you can't collapse it further
        if margin < 0.0:
            if numpy.array_equal(self.points[0], self.points[1]) and numpy.array_equal(self.points[0], self.points[2]):
                return Triangle(self.points)
        
        # get the centers of each side
        side_center_1 = numpy.average(self.points[[0,1]],0)
        side_center_2 = numpy.average(self.points[[1,2]],0)
        side_center_3 = numpy.average(self.points[[0,2]],0)
        side_centers = numpy.array([side_center_1, side_center_2, side_center_3])
        
        # get the vectors from the triangle center to each side center
        center_to_center_vec = side_centers - self.center()
        old_lengths = numpy.apply_along_axis(numpy.linalg.norm, 1, center_to_center_vec)
        new_lengths = old_lengths + margin
        
        # sanity check. None of the new lengths should be negative.
        # if any of them are, just set the whole triangle to a point
        # at the old triangle's center
        if new_lengths[0] < 0.0 or new_lengths[1] < 0.0 or new_lengths[2] < 0.0: return Triangle(numpy.array([self.center(), self.center(), self.center()]))
        
        # extend (or contract) these vectors farther out
        center_to_center_vec_normalized = (center_to_center_vec.T/old_lengths).T
        center_to_center_vec_new_length = self.center() + (center_to_center_vec_normalized.T * new_lengths).T
        
        # get an additional point on the extended side of the triangle
        second_points_on_new_side = center_to_center_vec_new_length + self.points - side_centers
        
        def seg_intersect(a1,a2, b1,b2):
            """Identify (or approximate) the point where two lines in 3D space intersect
            
            Arguments:
            a1 -- A 3x1 numpy array, a point on the first line
            a2 -- A 3x1 numpy array, a second point on the first line
            b1 -- A 3x1 numpy array, a point on the second line
            b2 -- A 3x1 numpy array, a second on the second line
            
            Returns:
            A 1x3 numpy array, the point (or an approximation thereof when an exact solution is not possible) where the defined lines intersect
            
            """
        
            # first, define the lines from the provided points
            pt1 = a1
            vec1 = a2-a1
            
            pt2 = b1
            vec2 = b2-b1
            
            # now get the points on the lines that are closest to each other
            coeffs = numpy.vstack((vec2, -vec1)).T
            best_sol_all = numpy.linalg.lstsq(coeffs, pt1-pt2)
            best_sol = best_sol_all[0]
            
            if best_sol_all[1][0] == 0.0: # an exact solution because the lines intersect
                return vec1 * best_sol[1] + pt1
            else: # return the average pt of the two points that are closest to each other
                close_pt1 = vec1 * best_sol[1] + pt1
                close_pt2 = vec2 * best_sol[0] + pt2
                
                return (close_pt1 + close_pt2) * 0.5 # return the average pt

        # get the corners of the new triangle
        pt1 = seg_intersect( center_to_center_vec_new_length[0],second_points_on_new_side[0], center_to_center_vec_new_length[1],second_points_on_new_side[1])
        pt2 = seg_intersect( center_to_center_vec_new_length[0],second_points_on_new_side[0], center_to_center_vec_new_length[2],second_points_on_new_side[2])
        pt3 = seg_intersect( center_to_center_vec_new_length[2],second_points_on_new_side[2], center_to_center_vec_new_length[1],second_points_on_new_side[1])
        
        new_pts = numpy.vstack((pt1, pt2, pt3))
        
        return Triangle(new_pts)
    
    def near_other_triangle(self, other_triangle, params):
        """Determines if another triangle is near this one
        
        Arguments:
        other_triangle -- A Triangle object, the other triangle to be evaluated
        params -- A dictionary, the user-specified command-line parameters
        
        Returns:
        A boolean, True of the two triangles are near each other, False otherwise
        
        """
        
        # check if the two triangles share a corner
        dists_between_triangle_pts = cdist(self.points, other_triangle.points)
        if True in (dists_between_triangle_pts == 0.0): return True # so they are adjacent
        
        # check if the distance to their centers is too close as well
        dist_between_center_points = numpy.linalg.norm(self.center() - other_triangle.center())
        if dist_between_center_points < params['triangle_center_proximity_cutoff_distance']: return True

############### GET COMMANDLINE PARAMETERS ###############

def get_commandline_parameters(argv):
    """Get the user-defined command-line parameters
    
    Returns:
    A dictionary, the user-specified command-line parameters
    
    """

    # first check if the user has requested the help file
    if '--help' in [t.lower() for t in argv]:
        print
        print "The initial lipid model"
        print "======================="
        print
        print textwrap.fill("--lipid_pdb_filename: This parameter specifies a PDB file containing an all-atom model of a planar lipid bilayer. LipidWrapper will wrap this lipid around the user-generated mesh. Example: --lipid_pdb_filename lipid.pdb", 70, subsequent_indent = "      ")
        print textwrap.fill("--lipid_headgroup_marker: A unique atom representing the headgroup of each lipid residue must be specified. The --lipid_headgroup_marker accepts a comma-separated lists of atom specifications (RESNAME_ATOMNAME). If either RESNAME or ATOMNAME is omitted, any value will be accepted. By default, LipidWrapper identifies lipid headgroups by looking for any atom named \"P\" (_P) or any atom named \"O3\" belonging to a cholesterol molecule (CHL1_O3). Example: --lipid_headgroup_marker \"_P,CHL1_O3\"", 70, subsequent_indent = "      ")
        print
        print "Methods for creating a surface mesh"
        print "==================================="
        print
        print textwrap.fill("--surface_equation: Generate a surface mesh from a python-formatted equation defining z, given x and y. The --min_x, --max_x, --min_y, and --max_y parameters are used to specify the region over which the function should be evaluated. The --step_x and --step_y parameters define the x-y distance between adjacent points. Python functions from the math, numpy, and scipy modules can be used. Example: --surface_equation \"z = 250*numpy.sin(x*x/60000 +y*y/60000)\"", 70, subsequent_indent = "      ")
        print textwrap.fill("--surface_filename: If this parameter specifies a file with the PDB extension, a surface mesh is generated from the coordinates of the PDB atoms. Example: --surface_filename mymesh.pdb", 70, subsequent_indent = "      ")
        print textwrap.fill("--surface_filename: If this parameter specifies a file with the DAE extension, the mesh points and triangulations will be taken from the file. Example: --surface_filename mymodel.dae", 70, subsequent_indent = "      ")
        print textwrap.fill("--surface_filename: If this parameter specifies a file that does not have the PDB extension, the file is assumed to be a gray-scale image, where black represents regions that are topologically low, and white represents regions that are topologically high. The --min_x, --max_x, --min_y, and --max_y parameters are used to specify the region where the mesh should be generated. The --step_x and --step_y parameters define the x-y distance between adjacent points. The --max_height parameter determines the height of the bilayer model at those locations where the image is white; black regions are assigned a height of 0. This feature is only available if the python PIL module has been installed on your system. Example: --surface_filename mymesh.png", 70, subsequent_indent = "      ")
        print
        print "Methods for resolving lipid clashes"
        print "==================================="
        print
        print textwrap.fill("--delete_clashing_lipids: It's common for lipids to sterically clash at the interface of two adjacent surface-mesh tessellated triangles. If this parameter is set to TRUE, any clashing lipids are deleted. Example: --delete_clashing_lipids TRUE", 70, subsequent_indent = "      ")
        print textwrap.fill("--clash_cutoff: If you do choose to delete clashing lipids, this parameter determines how close two atoms must be (in Angstroms) to constitute a steric clash. Example: --clash_cutoff 2.0", 70, subsequent_indent = "      ")
        print textwrap.fill("--clashing_potential_margin: Lipid clashes occur at the edges of adjacent tessellated triangles. If these triangles are very large, it's faster to only check for clashes and holes near the triangle edges. This variable specifies how far from the edges, in Angstroms, that LipidWrapper should look for clashes and holes. Example: --clashing_potential_margin 25.0", 70, subsequent_indent = "      ")
        print textwrap.fill("--fill_holes: Deleting lipids often leaves holes in the membrane. If this parameter is set to TRUE, LipidWrapper tries to fill the hole. Example: --fill_holes TRUE", 70, subsequent_indent = "      ")
        print textwrap.fill("--very_distant_lipids_cutoff: LipidWrapper determines if two lipids clash by comparing the distance between every atom in the first lipid with every atom in the second lipid. This can be computationally expensive. However, sometimes two lipids are so distant from each other, that it's obvious there are no clashes, making the pair-wise comparison unnecessary. Before performing this expensive pair-wise comparison, LipidWrapper calculates the distance between one atom of each lipid. If this distance is greater than this user-specified cutoff, the program will simply assume there are no clashes. WARNING: Remember to consider the width of your lipid bilayer when choosing this value. Adjacent lipids on opposite sides of the bilayer can seem distant when considering the distance between their headgroups, for example. Example: --very_distant_lipids_cutoff 50.0", 70, subsequent_indent = "      ")
        print textwrap.fill("--triangle_center_proximity_cutoff_distance: Lipid steric clashes/holes typically occur between lipids that belong to adjacent tessellated triangles. However, if tessellated triangles are small enough, clashes are possible between lipids that belong to non-adjacent triangles as well. Consequently, in addition to checking for adjacency, LipidWrapper also checks the distance between the triangle centers, using this user-specified value as a cutoff. Example: --triangle_center_proximity_cutoff_distance 50.0", 70, subsequent_indent = "      ")
        print textwrap.fill("--fill_hole_exhaustiveness: Essentially, how long LipidWrapper should try to fill the holes. Example: --fill_hole_exhaustiveness 10", 70, subsequent_indent = "      ")
        print textwrap.fill("--memory_optimization_factor: When the tessellated triangles are very large and consequently contain many individual lipids, the extensive pairwise distance comparisons required can result in memory errors. This parameter tells lipid Wrapper to divide the list of atoms being compared into smaller chunks. The pairwise distance comparison is performed piecewise on each chunk-chunk pair and so uses less memory, albeit at the expensive of speed. Only increase the value of this parameter if you run into memory errors. Example: --memory_optimization_factor 1", 70, subsequent_indent = "      ")
        print
        print "Additional options"
        print "=================="
        print
        print textwrap.fill("--number_of_processors: Using multiple processors can significantly increase the speed of the LipidWrapper algorithm. Example: --number_of_processors 8", 70, subsequent_indent = "      ")
        print textwrap.fill("--show_grid_points: Aside from producing PDB coordinates for lipid atoms, additional coordinates will be appended to the bottom of the output containing \"atoms\" named \"X\" that specify the location of the surface mesh points. Example: --show_grid_points TRUE", 70, subsequent_indent = "      ")
        print textwrap.fill("--create_triangle_tcl_file: A separate file named \"triangles.tcl\" will be generated containing a tcl script that can be run in VMD to visualize the mesh surface. Example: --create_triangle_tcl_file TRUE", 70, subsequent_indent = "      ")
        print textwrap.fill("--output_directory: If an output directory is specified, all LipidWrapper output files, as well as additional files representing the intermediate steps required to build the final bilayer, will be saved in that directory. Example: --output_directory ./my_output/", 70, subsequent_indent = "      ")
        print textwrap.fill("--use_disk_instead_of_memory: For very large systems, storing the growing model in memory can be problematic. If this parameter is set to TRUE, the growing model will be stored on the hard disk instead. However, expect longer execution times if this parameter is set to TRUE. Example: --use_disk_instead_of_memory TRUE", 70, subsequent_indent = "      ")
        print textwrap.fill("--compress_output: Depending on the user options selected, LipidWrapper output can require a lot of disk space. If this parameter is set to TRUE, the output will be automatically compressed using the gzip algorithm (Lempel-Ziv coding LZ77). The files can be uncompressed with the UNIX gunzip utility, or similar Windows-based packages. Example: --compress_output TRUE", 70, subsequent_indent = "      ")
        print
        print "Example"
        print "======="
        print
        print textwrap.fill('python lipidwrapper.py --surface_equation "z = 250*numpy.sin(x*x/60000 +y*y/60000) * (-numpy.sqrt(x*x+y*y)/(560 * numpy.sqrt(2)) + 1)" --min_x 500 --max_x 1000 --min_y 500 --max_y 1000 --step_x 25 --step_y 25 --lipid_pdb_filename lipid.pdb --lipid_headgroup_marker "_P,CHL1_O3" --delete_clashing_lipids TRUE --clash_cutoff 1.0 --fill_holes TRUE --fill_hole_exhaustiveness 10 > lipid_model.pdb', 70, subsequent_indent = "      ")
        print
        sys.exit(0)
    
    # defaults
    params = {}
    params['surface_filename'] = '' # could be a PDB or image file, depending on surface_source value
    params['surface_equation'] = 'z = 100*numpy.sin(x*x/60000 +y*y/60000) * (-numpy.sqrt(x*x+y*y)/(560 * numpy.sqrt(2)) + 1)' # used if surface_source is set to "EQUATION"
    params['min_x'] = 500 # used if surface_source is PNG or EQUATION
    params['max_x'] = 750 # used if surface_source is PNG or EQUATION
    params['min_y'] = 500 # used if surface_source is PNG or EQUATION
    params['max_y'] = 750 # used if surface_source is PNG or EQUATION
    params['step_x'] = 30 # used if surface_source is PNG or EQUATION
    params['step_y'] = 30 # used if surface_source is PNG or EQUATION
    params['max_height'] = 0 # used if surface_source is PNG
    params['lipid_pdb_filename'] = '' # the filename containing the small, planar lipid model
    params['lipid_headgroup_marker'] = '_P,CHL1_O3' # by default, any phosphate atom is considered a marker for the lipid headgroup, and also any O3 atom belonging to a cholesterol
    params['show_grid_points'] = 'FALSE'
    params['create_triangle_tcl_file'] = 'FALSE'
    params['delete_clashing_lipids'] = 'FALSE'
    params['use_disk_instead_of_memory'] = 'FALSE'
    params['clash_cutoff'] = 2.0
    params['fill_holes'] = 'TRUE'
    params['output_directory'] = ''
    params['fill_hole_exhaustiveness'] = 10
    params['number_of_processors'] = 1
    params['clashing_potential_margin'] = 25.0
    params['triangle_center_proximity_cutoff_distance'] = 50.0
    params['memory_optimization_factor'] = 1
    params['very_distant_lipids_cutoff'] = 50.0
    params['compress_output'] = "FALSE"    

    # get commandline parameters
    options, remainder = getopt.getopt(argv[1:], '', [ 'surface_filename=', 'surface_equation=', 'min_x=', 'max_x=', 'min_y=', 'max_y=', 'step_x=', 'step_y=', 'max_height=', 'lipid_pdb_filename=', 'lipid_headgroup_marker=', 'show_grid_points=', 'create_triangle_tcl_file=', 'delete_clashing_lipids=', 'clash_cutoff=', 'fill_holes=', 'fill_hole_exhaustiveness=', 'output_directory=', 'number_of_processors=', 'use_disk_instead_of_memory=', 'clashing_potential_margin=', 'triangle_center_proximity_cutoff_distance=', 'memory_optimization_factor=', 'very_distant_lipids_cutoff=', 'compress_output='])
    
    # set parameters to variables
    params_string = ['compress_output', 'surface_filename', 'surface_equation', 'lipid_pdb_filename', 'lipid_headgroup_marker', 'show_grid_points', 'create_triangle_tcl_file', 'delete_clashing_lipids', 'fill_holes', 'output_directory', 'use_disk_instead_of_memory']
    params_floats = ['very_distant_lipids_cutoff', 'memory_optimization_factor', 'triangle_center_proximity_cutoff_distance', 'clashing_potential_margin', 'min_x', 'max_x', 'min_y', 'max_y', 'step_x', 'step_y', 'max_height', 'clash_cutoff', 'fill_hole_exhaustiveness', 'number_of_processors']
    
    for opt, arg in options:
        opt = opt.replace('-','')
        if opt in params_floats: arg = float(arg)
        params[opt] = arg
    
    # some parameters should be integers
    params['fill_hole_exhaustiveness'] = int(params['fill_hole_exhaustiveness'])
    params['number_of_processors'] = int(params['number_of_processors'])
    params['memory_optimization_factor'] = int(params['memory_optimization_factor'])
    
    # directories should end in / or \ (depending on os)
    if params['output_directory'] != '' and params['output_directory'][-1:] != "/": params['output_directory'] = params['output_directory'] + "/"
    
    # check if running windows. If so, you can only use one processor
    if platform.system().lower() == "windows" and params['number_of_processors'] > 1:
        print "REMARK WARNING: Use of multiple processors is only supported on Linux and OS X."
        params['number_of_processors'] = 1

    # Print out header
    toprint = []
    #toprint.append("REMARK This lipid model was created using LipidWrapper")
    toprint.append("REMARK Parameters: (use the --help command-line parameter for further explanation)")
    for param in params.keys(): toprint.append("REMARK      " + param + ": " + str(params[param]))
    toprint.append("")
    
    # create the output directory if necessary, and write the parameters used to a file
    if params['output_directory'] == "": print "\n".join(toprint)
    else:
        try: os.mkdir(params['output_directory'])
        except: pass
        
        f = open(params['output_directory'] + 'parameters.input', 'w')
        f.write("\n".join(toprint))
        f.close()
    
    # in the case of the lipid_headgroup_marker, split it by the comma.
    params['lipid_headgroup_marker'] = [t.strip() for t in params['lipid_headgroup_marker'].split(',')]
    params['lipid_headgroup_marker'] = [(None, t.split('_')[0], None, t.split('_')[1]) for t in params['lipid_headgroup_marker']] # not that chain and resid are set to nothing, so any chain or resid will do

    # TRUE/FALSE answers need to be in caps.
    params['show_grid_points'] = params['show_grid_points'].upper()
    params['create_triangle_tcl_file'] = params['create_triangle_tcl_file'].upper()
    params['delete_clashing_lipids'] = params['delete_clashing_lipids'].upper()
    params['fill_holes'] = params['fill_holes'].upper()
    params['use_disk_instead_of_memory'] = params['use_disk_instead_of_memory'].upper()
    params['compress_output'] = params['compress_output'].upper()

    # specify the location of the temporary directory programatically
    params['memory_store_dir'] = params['output_directory'] + 'store_in_memory.tmp/'
    
    # now check a few of the parameters to make sure they're valid
    if not os.path.exists(params['lipid_pdb_filename']):
        print "ERROR: The file specified by the --lipid_pdb_filename parameter (" + params['lipid_pdb_filename'] + ") does not exist.\n"
        sys.exit(0)
    
    return params

################## FILE HANDLING FUNCTIONS ##################

def save_pickle(item, params, an_id=''):
    """Save an object to a pickle file
    
    Arguments:
    item -- The object to be saved
    item -- A dictionary, the user-specified comand-line parameters
    an_id -- An optional string, the id of the pickled object
    
    Returns:
    A string, the id of the current pickle. If the user-specified an_id is '', then a random an_id is generated.
    
    """

    # make an initial guess at hte pickle name
    filename = params['memory_store_dir'] + an_id + ".pickle"
    
    # need to make up an id, change filename, if an_id is not specified
    if an_id == '':
        an_id = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(12))
        filename = params['memory_store_dir'] + an_id + ".pickle"
        
        # keep trying until you come up with a unique an_id. Almost certainly on first try.
        while os.path.exists(filename):
            an_id = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(12))
            filename = params['memory_store_dir'] + an_id + ".pickle"
    
    while True: # keep trying to write until you succeed
        try:
            pickle.dump(item, gzip.open(filename, "wb" ), protocol=2)
            break
        except: pass
    
    return an_id # let user know what an_id you settled on

def load_pickle(an_id, params):
    """Save an object to a pickle file
    
    Arguments:
    an_id -- A string, the id of the pickled object to load
    params -- A dictionary, the user-specified comand-line parameters
    
    Returns:
    A python object, loaded from the pickle file
    
    """
    
    while True: # keep trying to return the pickle until you succeed.
        try: return pickle.load(gzip.open(params['memory_store_dir'] + an_id + ".pickle", "rb"))
        except: pass

def openfile(filename, mode, params):
    if params['compress_output'] == "TRUE": # open a gzip file
        return gzip.open(filename + ".gz", mode)
    else: # open a regular file
        return open(filename, mode)

################## LOAD AND POSITION LIPIDS ONTO THE 3D SURFACE ##################

def load_mesh_points_and_triangulations(params):
    """Load mesh points and perform a tesselation/triangulation if neceesary
    
    Arguments:
    params -- A dictionary, the user-specified comand-line parameters
    
    Returns:
    A list of Triangle objects, the tesselation/triangulation
    
    """

    # load the mesh points from whatever source the user specifried
    pts = Molecule()
    all_triangles = None
    
    # could be from a PDB file
    if params['surface_filename'][-3:].upper() == 'PDB': pts.load_pdb(params['surface_filename'])
    
    # could be from a blender-exported DAE file
    elif params['surface_filename'][-3:].upper() == 'DAE': # this is a Collada mesh generated by blender
        f = open(params['surface_filename'], 'r')
        while True:
            line = f.readline()
            if len(line) == 0: break # EOF
            if "<float_array" in line and "mesh-positions-array" in line: # so this is the line containing points
                pts_str = line.split(">")[1].split("<")[0].strip()
                while "  " in pts_str: pts_str = pts_str.replace('  ',' ')
                pts_float = [float(t) for t in pts_str.split(" ")]
                pts_list = [] # this is going to be so small that using python list is ok
                for t in range(0,len(pts_float),3): pts_list.append([pts_float[t], pts_float[t+1], pts_float[t+2]])
                pts.all_atoms_numpy = numpy.array(pts_list)
                
            if "<polylist" in line:
                # now figure out how many inputs there are and which one is the VERTEX
                line = f.readline()
                count_input = 0
                vertex_index = -1
                while "<input" in line:
                    count_input = count_input + 1
                    if "VERTEX" in line: vertex_index = count_input - 1
                    line = f.readline()
                
                # so the next line should be vcounts
                vcounts = line.split(">")[1].split("<")[0].strip()
                while "  " in vcounts: vcounts = vcounts.replace('  ',' ')
                vcounts = [int(t) for t in vcounts.split(" ")]
                all_threes = True
                for t in vcounts:
                    if t != 3:
                        all_threes = False
                        break
                if all_threes == False:
                    print "This mesh has not been triangulated. We recommend using blender. Press Ctrl-T in Edit Mode with the mesh selected."
                    sys.exit(0)
                
                # the next line has the triangles
                line = f.readline()
                verts = line.split(">")[1].split("<")[0].strip()
                while "  " in verts: verts = verts.replace('  ',' ')
                verts = [int(t) for t in verts.split(" ")]
                all_triangles = []
                for t in range(0,len(verts),3*count_input):
                    pt1_index = verts[t + vertex_index]
                    pt2_index = verts[t+count_input + vertex_index]
                    pt3_index = verts[t+count_input*2 + vertex_index]
                    
                    pt1 = pts.all_atoms_numpy[pt1_index]
                    pt2 = pts.all_atoms_numpy[pt2_index]
                    pt3 = pts.all_atoms_numpy[pt3_index]
                    
                    all_triangles.append([pt1, pt2, pt3])
                all_triangles = numpy.array(all_triangles)
                    
        f.close()
    
    # could be from some image
    elif params['surface_filename'] != '': # so it must be an image
        
        width = params['max_x'] - params['min_x']
        height = params['max_y'] - params['min_y']
    
        try: from PIL import Image
        except ImportError:
            print "Sorry, but to use an image as the surface source, PIL must be installed..."
            sys.exit(0)
    
        pic = Image.open(params['surface_filename'])
        pic = pic.resize((int(width), int(height)), Image.NEAREST)
        pic = numpy.array(pic)
    
        pts_list = []
        
        for x in numpy.arange(0, width, params['step_x']):
            for y in numpy.arange(0, height, params['step_y']):
                z = params['max_height'] * pic[x,y,0]/255.0 # 0 because it's R, G, B, alpha, and images should be greyscale
                pts_list.append(numpy.array([x + params['min_x'], y + params['min_y'], z]))
        pts.all_atoms_numpy = numpy.array(pts_list)
    
    # could be from an equation
    else: # so derive it from an equation
        pts_list = []
        for x in numpy.arange(params['min_x'], params['max_x'], params['step_x']):
            for y in numpy.arange(params['min_y'], params['max_y'], params['step_y']):
                z = 0.0
                exec params['surface_equation'] # to get the z value
                if not math.isnan(z): pts_list.append([x,y,z])
        pts.all_atoms_numpy = numpy.array(pts_list)
    
    # for everything but the DAE input, a tesselation/triangulation must also be performed
    if all_triangles == None: # so you need to get the triangulation
    
        # project the mesh onto the x-y plane (so it's important the it be oriented so that positive z is up)
        flatten = pts.all_atoms_numpy.copy()
        flatten = flatten[:,0:2]
        
        # now tesselate the 2D points
        tri1 = Delaunay(flatten)
        
        # project the points back onto the mesh surface (3d trinagles)
        all_triangles = []
        for ia, ib, ic in tri1.vertices: all_triangles.append([pts.all_atoms_numpy[ia], pts.all_atoms_numpy[ib], pts.all_atoms_numpy[ic]])
        all_triangles = numpy.array(all_triangles)
        
    # convert this list of triangle points into a list of Triangle objects
    gc.disable()
    all_triangles_obj = []
    for tri in all_triangles:
        tri2 = Triangle(tri)
        all_triangles_obj.append(tri2)
    gc.enable()
    
    return all_triangles_obj

def load_lipid_model(params):
    """Load the user-provided planar bilayer model
    
    Arguments:
    params -- A dictionary, the user-specified command-line parameters
    
    Returns:
    A Molecule object containing the planar bilayer model, a 1x3 numpy array representing the minimum corner of a bounding box, and a 1x3 numpy array representing the maximum corner of a bounding box.
    
    """

    # Now load the lipid pdb, which should also be oriented so Z is up
    # The lipid must have unique resids too. This is standard.
    lipid = Molecule()
    lipid.load_pdb(params['lipid_pdb_filename'])
    
    # get the coordinates of the lipid headgroups
    headgroup_markers_coors = lipid.all_atoms_numpy[lipid.get_indices_of_mask_match(params['lipid_headgroup_marker'])]
    
    # center the lipids and the headgroups so the x-y plane bisects the bilayer
    delta_z = numpy.array([0,0,numpy.average(headgroup_markers_coors[:,2])])
    lipid.all_atoms_numpy = lipid.all_atoms_numpy - delta_z
    headgroup_markers_coors = headgroup_markers_coors - delta_z
    
    # get the coordinates of a bounding box that encompasses the entire planar-bilayer model
    min_headgroups = numpy.amin(headgroup_markers_coors,0) # [x,y,z]
    max_headgroups = numpy.amax(headgroup_markers_coors,0) # [x,y,z]
    
    return lipid, min_headgroups, max_headgroups

def angle_between(v1, v2):
    """Calculates the angle between two vectors
    
    Arguments:
    v1 -- A 1x3 numpy array, the first vector
    v2 -- A 1x3 numpy array, the second vector
    
    Returns:
    A float, the angle between the two vectors in radians.
    If the two vectors are equal, the string "NORMALIZED VECTORS EQUAL!" is returned instead.
    
    """

    # get the unit vectors
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    
    # if the unit vectors are the same, you'll get an error.
    # return this string instead.
    if numpy.array_equal(v1_u, v2_u): return "NORMALIZED VECTORS EQUAL!"
    if numpy.linalg.norm(v1_u - v2_u) < 1e-7: return "NORMALIZED VECTORS EQUAL!"
    
    # if two vectors are pointing in the opposite directory, just return pi
    # This check is needed because sometimes numpy.dot(v1_u, v2_u) is actually slightly more than -1.0, giving an error
    if numpy.array_equal(v1_u, -v2_u): return numpy.pi
    if numpy.linalg.norm(v1_u + v2_u) < 1e-7: return numpy.pi
    
    # calculate the angle
    angle = numpy.arccos(numpy.dot(v1_u, v2_u))

    # if there's an error, modify the output
    if math.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return numpy.pi
    
    return angle

def unit_vector(vector):
    """Take a vector and scales it so its length is 1.0
    
    Arguments:
    vector -- A 1x3 numpy array, the vector to be scaled
    
    Returns:
    A 1x3 numpy array, the scaled vector
    
    """

    return vector/numpy.linalg.norm(vector)

def get_transformation_data(tethers1_static, tethers2_dynamic):
    """Determines how to translate and rotate one set of points onto another
    
    Arguments:
    tethers1_static -- An nx3 numpy array, the points that will remain fixed
    tethers2_dynamic -- An nx3 numpy array, the points that will be moved on to the points in tethers1_static. Note that both inputs must have the same point order.
    
    Returns:
    A tuple containing a 1x3 numpy array representing the center of the tethers2_dynamic points, a quarterion describing the required translation/rotation, and a 1x3 numpy array representing the center of the tethers1_static points
    
    """

    # Get the transformation matrix to move the dynamic_template onto the guide_static_template
    center_tethers1_pdb = numpy.mean(tethers1_static, 0)
    center_tethers2_pdb = numpy.mean(tethers2_dynamic, 0)

    # Translate com of tether molecule to origin
    tethers1_static = tethers1_static - center_tethers1_pdb
    tethers2_dynamic = tethers2_dynamic - center_tethers2_pdb

    # Get optimal rotation
    M = numpy.dot(numpy.transpose(tethers1_static), tethers2_dynamic)
    
    # Create symmetric 4x4 matrix K from M
    K = numpy.array([[M[0,0] + M[1,1] + M[2,2], M[1,2] - M[2,1], M[2,0] - M[0,2], M[0,1] - M[1,0]],
                     [M[1,2] - M[2,1], M[0,0] - M[1,1] - M[2,2], M[1,0] + M[0,1], M[2,0] + M[0,2]],
                     [M[2,0] - M[0,2], M[1,0] + M[0,1], M[1,1] - M[0,0] - M[2,2], M[1,2] + M[2,1]],
                     [M[0,1] - M[1,0], M[2,0] + M[0,2], M[1,2] + M[2,1], M[2,2] - M [0,0] - M[1,1]]])
            
    # Find eigenvector associated with the most positive eigenvalue of K.  Multiple quaternions can
    E,V = numpy.linalg.eig(K)
    index = numpy.argmax(E)
    eigenvector = V[:,index]
    rot_quat = Quaternion(eigenvector[0], eigenvector[1], eigenvector[2], eigenvector[3])
    
    return (center_tethers2_pdb, rot_quat, center_tethers1_pdb)

def apply_transformation(mol, transform_data):
    """Translate and rotate a molecule
    
    Arguments:
    mol -- A Molecule object, to be translated and rotated
    transform_data -- A tuple containing the output of the get_transformation_data() function

    """

    center_dynamic_pdb = transform_data[0]
    rot_quat = transform_data[1]
    center_static_pdb = transform_data[2]

    mol.all_atoms_numpy = mol.all_atoms_numpy - center_dynamic_pdb # move to center
    mol.rotate_mol_quat(rot_quat) # rotate
    mol.all_atoms_numpy = mol.all_atoms_numpy + center_static_pdb # move to new location

def position_lipid_model_on_triangulated_tiles(params, lipid, all_triangles, min_headgroups, max_headgroups):
    """Position portions of the user-specified planar bilayer onto the tessellated/triangled surface

    Arguments:
    params -- A dictionary, the user-specified command-line parameters
    lipid -- A Molecule object containing the initial user-specified planar bilayer model
    all_triangles -- A list of Triangle objects, the tesselation/triangulation
    min_headgroups -- A 1x3 numpy array representing the minimum corner of a bounding box
    max_headgroups -- A 1x3 numpy array representing the maximum corner of a bounding box

    Returns:
    A list of tuples, where each tuple contains a Triangle object and a list of lipid molecules (Molecule objects) that belong to that triangle
    
    """

    # position the lipid bilayers on the tiles
    lipid.set_undo_point()
    
    class position_lipids_multiprocessing(general_task):
        """A class for positioning lipid models onto a 3D mesh"""

        def value_func(self, item, results_queue): # so overwriting this function
            """Position lipid models onto a 3D mesh

            Arguments:
            item -- A list or tuple, the input data required for the calculation
            results_queue -- A multiprocessing.Queue() object for storing the calculation output
            
            """

            params = item[0]
            lipid = item[1]
            tile_triangle = item[2]
            min_headgroups = item[3]
            max_headgroups = item[4]
            current_index = item[5]
            
            self.print_star_if_appropriate(current_index)
            
            # get the tessellated 3D triangle center and radii (r1, r2, r3, overall radius)
            tile_r1 = tile_triangle.radii()[0][0]
            tile_r2 = tile_triangle.radii()[0][1]
            tile_r3 = tile_triangle.radii()[0][2]
        
            # get the interior angles of the tessellated 3D triangle
            u = tile_triangle[0] - tile_triangle.center()
            v = tile_triangle[1] - tile_triangle.center()
            r1_r2_angle = angle_between(u, v)
            
            u = tile_triangle[0] - tile_triangle.center()
            v = tile_triangle[2] - tile_triangle.center()
            r1_r3_angle = angle_between(u, v)
            
            # calculate the appropriate location of an identical triangle on the user-defined planar-bilayer model. ignoring z for now (since lipid is oriented along x-y plane)
            lipid_triangle_pt1 = numpy.array([0.0, tile_r1, 0.0])
            lipid_triangle_pt2 = numpy.array([-tile_r2 * numpy.cos(r1_r2_angle - (numpy.pi/2.0)), -tile_r2 * numpy.sin(r1_r2_angle - (numpy.pi/2.0)), 0.0])
            lipid_triangle_pt3 = numpy.array([tile_r3 * numpy.cos(r1_r3_angle - (numpy.pi/2.0)), -tile_r3 * numpy.sin(r1_r3_angle - (numpy.pi/2.0)), 0.0])
            
            # rotate the identical triangle randomly about z axis
            random_theta = numpy.random.uniform(0, 2.0 * numpy.pi)
            rot_matrix = numpy.array([[numpy.cos(random_theta), numpy.sin(random_theta), 0.0], [-numpy.sin(random_theta), numpy.cos(random_theta), 0.0], [0.0, 0.0, 1.0]])
            lipid_triangle_pt1 = numpy.dot(rot_matrix, lipid_triangle_pt1)
            lipid_triangle_pt2 = numpy.dot(rot_matrix, lipid_triangle_pt2)
            lipid_triangle_pt3 = numpy.dot(rot_matrix, lipid_triangle_pt3)
            
            # now pick an appropriate location, and move the identical triangle to that location along x-y plane
            lipid_triangle_center_x = numpy.random.uniform(min_headgroups[0]+tile_triangle.max_radius(), max_headgroups[0]-tile_triangle.max_radius())
            lipid_triangle_center_y = numpy.random.uniform(min_headgroups[1]+tile_triangle.max_radius(), max_headgroups[1]-tile_triangle.max_radius())
            lipid_triangle_center = numpy.array([lipid_triangle_center_x, lipid_triangle_center_y, 0.0])
            lipid_triangle_pt1 = lipid_triangle_pt1 + lipid_triangle_center
            lipid_triangle_pt2 = lipid_triangle_pt2 + lipid_triangle_center
            lipid_triangle_pt3 = lipid_triangle_pt3 + lipid_triangle_center
            
            lipid_triangle = numpy.array([lipid_triangle_pt1, lipid_triangle_pt2, lipid_triangle_pt3])
            
            # now get the transformation matrix to transform the identical triangle onto the 3D tessellated triangle
            # and apply the transformation
            transform_data = get_transformation_data(tile_triangle.points, lipid_triangle)
            lipid.undo()
            apply_transformation(lipid,transform_data)
    
            # project the headgroups of the planar bilayer onto the identical triangle
            headgroup_markers_indices = lipid.get_indices_of_mask_match(params['lipid_headgroup_marker'])
            headgroup_markers_coors = lipid.all_atoms_numpy[headgroup_markers_indices]
            headgroup_marker_proj_coors = tile_triangle.project_points_onto_triangle(headgroup_markers_coors)
            
            # identify the indices of the head groups that fall within the triangle when projected
            headgroup_indices_to_keep = get_numpy_slice(headgroup_markers_indices, tile_triangle.get_indices_of_points_within_triangle_boundaries(headgroup_marker_proj_coors))
            
            # there are legitamte circumstances where no headgroups should be retained (i.e., no headgroups fall within the triangle boundaries)
            # for example, sometimes the scipy tesselation could produce "triangles" that are actually lines.
            # so abort efforts if this is the case.
            if len(headgroup_indices_to_keep) == 0 and len(headgroup_marker_proj_coors) !=0:
                molecules_in_this_triangle = []
                if params['use_disk_instead_of_memory'] == "TRUE":
                    an_id = save_pickle(molecules_in_this_triangle, params)
                    self.results.append((tile_triangle, an_id))
                else: self.results.append((tile_triangle, molecules_in_this_triangle))
                return

            # identify the indices of the head groups that fall within a smaller triangle, that are too far from the edges to worry about steric clashes in subsequent steps.
            # this speeds up subsequent steps significantly
            smaller_tri = tile_triangle.new_triangle_expanded_by_margin(-params['clashing_potential_margin'])
            todel = smaller_tri.get_indices_of_points_within_triangle_boundaries(headgroup_marker_proj_coors)
            headgroup_indices_not_in_triangle_margin = get_numpy_slice(headgroup_markers_indices, todel)
            
            # identify the indices of the head groups that fall within an even smaller triangle, called the submargin.
            # identifying these headgroups will speed up the lipid-filling step later. Lipids will be placed in the margin, and the steric clashes are checked with the lipids of the submargin (as well as the margins of neighboring triangles)
            smaller_still_tri = smaller_tri.new_triangle_expanded_by_margin(-params['clashing_potential_margin'])
            todel = smaller_still_tri.get_indices_of_points_within_triangle_boundaries(headgroup_marker_proj_coors)
            headgroup_indices_not_in_triangle_margin_or_submargin = get_numpy_slice(headgroup_markers_indices, todel)

            # now get the entire lipids corresponding to these headgroups
            molecules_in_this_triangle = []

            # first, identify all residues
            atom_counts = len(lipid.atom_inf_resids)
            
            all_ids = numpy.core.defchararray.add(lipid.atom_inf_string_vals[:,0], numpy.array(['_'] * atom_counts))
            all_ids = numpy.core.defchararray.add(all_ids, lipid.atom_inf_resids)
            all_ids = numpy.core.defchararray.add(all_ids, numpy.array(['_'] * atom_counts))
            all_ids = numpy.core.defchararray.add(all_ids, lipid.atom_inf_string_vals[:,1])
            
            # now identify all residues to keep
            atom_counts = len(headgroup_indices_to_keep)
            if atom_counts != 0:
                hg_ids = numpy.core.defchararray.add(lipid.atom_inf_string_vals[headgroup_indices_to_keep,0], numpy.array(['_'] * atom_counts))
                hg_ids = numpy.core.defchararray.add(hg_ids, lipid.atom_inf_resids[headgroup_indices_to_keep])
                hg_ids = numpy.core.defchararray.add(hg_ids, numpy.array(['_'] * atom_counts))
                hg_ids = numpy.core.defchararray.add(hg_ids, lipid.atom_inf_string_vals[headgroup_indices_to_keep,1])
            else:
                hg_ids = numpy.array([])
                
            # now identify all residues to never delete (inside inner triangle)
            atom_counts = len(headgroup_indices_not_in_triangle_margin)
            if atom_counts != 0:
                hg_ids_not_in_triangle_margin = numpy.core.defchararray.add(lipid.atom_inf_string_vals[headgroup_indices_not_in_triangle_margin,0], numpy.array(['_'] * atom_counts))
                hg_ids_not_in_triangle_margin = numpy.core.defchararray.add(hg_ids_not_in_triangle_margin, lipid.atom_inf_resids[headgroup_indices_not_in_triangle_margin])
                hg_ids_not_in_triangle_margin = numpy.core.defchararray.add(hg_ids_not_in_triangle_margin, numpy.array(['_'] * atom_counts))
                hg_ids_not_in_triangle_margin = numpy.core.defchararray.add(hg_ids_not_in_triangle_margin, lipid.atom_inf_string_vals[headgroup_indices_not_in_triangle_margin,1])
            else:
                hg_ids_not_in_triangle_margin = numpy.array([])
            
            # now identify all residues that are even interior to the submargin
            atom_counts = len(headgroup_indices_not_in_triangle_margin_or_submargin)
            if atom_counts != 0:
                hg_ids_not_in_triangle_margin_or_submargin = numpy.core.defchararray.add(lipid.atom_inf_string_vals[headgroup_indices_not_in_triangle_margin_or_submargin,0], numpy.array(['_'] * atom_counts))
                hg_ids_not_in_triangle_margin_or_submargin = numpy.core.defchararray.add(hg_ids_not_in_triangle_margin_or_submargin, lipid.atom_inf_resids[headgroup_indices_not_in_triangle_margin_or_submargin])
                hg_ids_not_in_triangle_margin_or_submargin = numpy.core.defchararray.add(hg_ids_not_in_triangle_margin_or_submargin, numpy.array(['_'] * atom_counts))
                hg_ids_not_in_triangle_margin_or_submargin = numpy.core.defchararray.add(hg_ids_not_in_triangle_margin_or_submargin, lipid.atom_inf_string_vals[headgroup_indices_not_in_triangle_margin_or_submargin,1])
            else:
                hg_ids_not_in_triangle_margin_or_submargin = numpy.array([])
            
            # remove the lipids that are beyond the bounds of this triangle to speed up subsequent searching
            # Find the indices of elements of all_ids that are in hg_ids.
            iall_ids = numpy.in1d(all_ids.ravel(), hg_ids).reshape(all_ids.shape)
            indices_of_lipids_to_keep = numpy.where(iall_ids)[0]
            lipid = lipid.portion_of(indices_of_lipids_to_keep)
            all_ids = all_ids[indices_of_lipids_to_keep]
            
            # find where to split this lipid model into individual lipids
            all_ids_offset = all_ids.copy()
            all_ids_offset = numpy.append(all_ids_offset, all_ids_offset[0])
            all_ids_offset = all_ids_offset[numpy.arange(1, len(all_ids_offset), dtype=int)]
            indices_of_each_lipid = 1 + numpy.nonzero(numpy.logical_not(all_ids == all_ids_offset))[0]
            indices_of_each_lipid = numpy.insert(indices_of_each_lipid,0,0)
            
            # now move each individual lipid into its own object and append
            gc.disable() # because appending objects
            for t in range(len(indices_of_each_lipid)-1):
                start_index = indices_of_each_lipid[t]
                end_index = indices_of_each_lipid[t+1]
                atom_range = numpy.arange(start_index, end_index, dtype=int)
                single_lipid = lipid.portion_of(atom_range)
                single_lipid.in_triangle_submargin = False
                single_lipid.in_triangle_margin = True
                
                if all_ids[start_index] in hg_ids_not_in_triangle_margin:
                    single_lipid.in_triangle_margin = False
                    
                    if not all_ids[start_index] in hg_ids_not_in_triangle_margin_or_submargin:
                        single_lipid.in_triangle_submargin = True
                
                molecules_in_this_triangle.append(single_lipid)
            
            if params['use_disk_instead_of_memory'] == "TRUE":
                an_id = save_pickle(molecules_in_this_triangle, params)
                self.results.append((tile_triangle, an_id))
            else: self.results.append((tile_triangle, molecules_in_this_triangle))

            gc.enable()
    
    input_array = [(params, lipid, tile_triangle, min_headgroups, max_headgroups, index+1) for index, tile_triangle in enumerate(all_triangles)]
    tmp = multi_threading(input_array, params['number_of_processors'], position_lipids_multiprocessing, params, "REMARK ")
    
    molecules_by_triangle = tmp.results
    
    return molecules_by_triangle

################## REMOVE STATIC CLASHES ##################

def remove_steric_clashes(molecules_by_triangle, params):
    """Remove lipids that have steric clashes

    Arguments:
    molecules_by_triangle -- A list of tuples, where each tuple contains a Triangle object and a list of lipid molecules (Molecule objects) that belong to that triangle
    params -- A dictionary, the user-specified command-line parameters

    """

    class remove_clashes_multiprocessing(general_task):
        """A class for identifying lipids that have steric clashes"""
        
        def value_func(self, item, results_queue): # so overwriting this function
            """Identify lipids that have steric clashes

            Arguments:
            item -- A list or tuple, the input data required for the calculation
            results_queue -- A multiprocessing.Queue() object for storing the calculation output
            
            """

            triangle_index1 = item[0]
            triangle_index2 = item[1]
            triangle1 = item[2]
            triangle1_lipids = item[3]
            triangle2 = item[4]
            triangle2_lipids = item[5]
            params = item[6]
            
            self.print_star_if_appropriate(item[7])
            
            if triangle1.near_other_triangle(triangle2, params): # so only the lipids of proximate triangles are considered
                
                if params['use_disk_instead_of_memory'] == "TRUE": triangle1_lipids = load_pickle(triangle1_lipids, params)
                if params['use_disk_instead_of_memory'] == "TRUE": triangle2_lipids = load_pickle(triangle2_lipids, params)
                
                clash_map = {}
                
                # now generate a numpy array containing all the headgroups of the lipids of each triangle
                triangle1_headgroups = numpy.empty((len(triangle1_lipids), 3))
                for lipid_index, lipid in enumerate(triangle1_lipids):
                    headgroup_loc = lipid.all_atoms_numpy[lipid.get_headgroup_index(params['lipid_headgroup_marker'])]
                    triangle1_headgroups[lipid_index][0] = headgroup_loc[0]
                    triangle1_headgroups[lipid_index][1] = headgroup_loc[1]
                    triangle1_headgroups[lipid_index][2] = headgroup_loc[2]
                    
                triangle2_headgroups = numpy.empty((len(triangle2_lipids), 3))
                for lipid_index, lipid in enumerate(triangle2_lipids):
                    headgroup_loc = lipid.all_atoms_numpy[lipid.get_headgroup_index(params['lipid_headgroup_marker'])]
                    triangle2_headgroups[lipid_index][0] = headgroup_loc[0]
                    triangle2_headgroups[lipid_index][1] = headgroup_loc[1]
                    triangle2_headgroups[lipid_index][2] = headgroup_loc[2]
                    
                # get the indices of all the lipids in the margin of the first lipid
                indices_in_margin1 = []
                for idx, lip in enumerate(triangle1_lipids):
                    if lip.in_triangle_margin == True: indices_in_margin1.append(idx)
                indices_in_margin1 = numpy.array(indices_in_margin1)
                
                if len(indices_in_margin1) > 0: # so there are some lipids in the margin
                    
                    # get the indices of all the lipids in the margin of the second lipid
                    indices_in_margin2 = []
                    for idx, lip in enumerate(triangle2_lipids):
                        if lip.in_triangle_margin == True: indices_in_margin2.append(idx)
                    indices_in_margin2 = numpy.array(indices_in_margin2)
                    
                    if len(indices_in_margin2) > 0: # so there are some lipids in the margin
                    
                        # now, look at distances between all headgroups in margin to identify ones that are close enough to potentially clash
                        dists = cdist(triangle1_headgroups[indices_in_margin1], triangle2_headgroups[indices_in_margin2])
                        dists = dists < params['very_distant_lipids_cutoff']
                        indices_to_look_at1 = indices_in_margin1[numpy.nonzero(dists)[0]]
                        indices_to_look_at2 = indices_in_margin2[numpy.nonzero(dists)[1]]
                        
                        # now do a pairwise comparison of all lipids, looking for clashes
                        for t in range(len(indices_to_look_at1)):
                            lipid_index1 = indices_to_look_at1[t]
                            lipid1 = triangle1_lipids[lipid_index1]
                            if lipid1.in_triangle_margin == True:
                                lipid_index2 = indices_to_look_at2[t]
                                lipid2 = triangle2_lipids[lipid_index2]
                                if lipid2.in_triangle_margin == True:
                                    if two_lipids_clash(lipid1, lipid2, params['clash_cutoff'], 1, params, False):
                                        
                                        # there's a clash. update the clash map
                                        
                                        try:
                                            clash_map[(lipid_index1, triangle_index1)].append((lipid_index2, triangle_index2))
                                        except:
                                            clash_map[(lipid_index1, triangle_index1)] = []
                                            clash_map[(lipid_index1, triangle_index1)].append((lipid_index2, triangle_index2))
                
                                        try:
                                            clash_map[(lipid_index2, triangle_index2)].append((lipid_index1, triangle_index1))
                                        except:
                                            clash_map[(lipid_index2, triangle_index2)] = []
                                            clash_map[(lipid_index2, triangle_index2)].append((lipid_index1, triangle_index1))
                        self.results.append(clash_map)

    # generate a clash map, whcih specifies which lipids clash with each other
    some_input = []
    gc.disable() # because appending complex objects to a list
    t = 0
    for triangle_index1 in range(len(molecules_by_triangle)-1):
        for triangle_index2 in range(triangle_index1 + 1, len(molecules_by_triangle)):
            
            t = t + 1

            triangle1 = molecules_by_triangle[triangle_index1][0]
            triangle1_lipids = molecules_by_triangle[triangle_index1][1]
            
            triangle2 = molecules_by_triangle[triangle_index2][0]
            triangle2_lipids = molecules_by_triangle[triangle_index2][1]
            
            some_input.append((triangle_index1, triangle_index2, triangle1, triangle1_lipids, triangle2, triangle2_lipids, params, t))
    gc.enable()
    
    tmp = multi_threading(some_input, params['number_of_processors'], remove_clashes_multiprocessing, params, "REMARK            (step 1) ")
    
    # now combine all the clash maps into one
    clash_map = {}
    for amap in tmp.results:
        for akey in amap.keys():
            try:
                clash_map[akey].extend(amap[akey])
            except:
                clash_map[akey] = []
                clash_map[akey].extend(amap[akey])

    # how we actually delete molecules
    # keep deleting molecules that clash until everything's resolved.
    # this part runs on one processor, but clash detection ran on multiple ones
    
    # start eliminating molecules until all clashes in the clashmap are resolved, starting with the molecule that has the most clashes
    gc.disable()
    lipids_to_delete = {}
    while len(clash_map) > 0: # keep going until all clashes are resolved
        
        # identify the molecule that makes the most clashes
        most_clashes = 0
        most_clashes_mol_index = -1
        for mol_index in clash_map.keys():
            num_clashes = len(clash_map[mol_index])
            if num_clashes > most_clashes:
                most_clashes = num_clashes
                most_clashes_mol_index = mol_index
        
        # now go through each of the ones it clashes with and remove it from their lists
        for other_clasher in clash_map[most_clashes_mol_index]:
            clash_map[other_clasher].remove(most_clashes_mol_index)
            if len(clash_map[other_clasher]) == 0: del clash_map[other_clasher]
            
        # now assign its value in the lipids lists to None, to be deleted later
        try: lipids_to_delete[most_clashes_mol_index[1]].append(most_clashes_mol_index[0])
        except:
            lipids_to_delete[most_clashes_mol_index[1]] = []
            lipids_to_delete[most_clashes_mol_index[1]].append(most_clashes_mol_index[0])
            
        # now remove this one from clash_map as well
        del clash_map[most_clashes_mol_index]
    gc.enable()

    # now delete the lipids that have been marked for deletion by being assigned a value of None
    class remove_clashes2_multiprocessing(general_task):
        """A class for removing lipids that have steric clashes"""
        
        def value_func(self, item, results_queue): # so overwriting this function
            """Remove lipids that have steric clashes

            Arguments:
            item -- A list or tuple, the input data required for the calculation
            results_queue -- A multiprocessing.Queue() object for storing the calculation output
            
            """

            triangle_index = item[0]
            lipid_indices_to_delete = item[1]
            lipids_list = item[2]
            params = item[4]
            
            self.print_star_if_appropriate(item[3])
            
            if params['use_disk_instead_of_memory'] == "TRUE": lipids = load_pickle(lipids_list, params)
            else: lipids = lipids_list
                
            for lipid_index in lipid_indices_to_delete: lipids[lipid_index] = None
            while None in lipids: lipids.remove(None)
            
            if params['use_disk_instead_of_memory'] == "TRUE": save_pickle(lipids, params, lipids_list)
            else:
                gc.disable()
                self.results.append((triangle_index, lipids))
                gc.enable()

    # now actually go through and remove all the "lipids" that have been asigned None
    some_input = []
    gc.disable()
    t = 0
    for triangle_index in lipids_to_delete.keys():
        t = t + 1
        lipid_indices_to_delete = lipids_to_delete[triangle_index]
        some_input.append((triangle_index, lipid_indices_to_delete, molecules_by_triangle[triangle_index][1], t, params))
    gc.enable()
    
    tmp = multi_threading(some_input, params['number_of_processors'], remove_clashes2_multiprocessing, params, "REMARK            (step 2) ")
    
    # update molecules_by_triangle
    if params['use_disk_instead_of_memory'] != "TRUE": # so you need to reconstruct molecules_by_triangle
        for triangle_id, lipids_list in tmp.results:
            if len(molecules_by_triangle[triangle_id]) == 2: molecules_by_triangle[triangle_id] = (molecules_by_triangle[triangle_id][0], lipids_list)
            else: molecules_by_triangle[triangle_id] = (molecules_by_triangle[triangle_id][0], lipids_list, molecules_by_triangle[triangle_id][2])

def two_lipids_clash(mol1, mol2, cutoff, num_sub_partitions, params, very_large_distance_check=True):
    """Determine whether two lipid molecules clash

    Arguments:
    mol1 -- A Molecule object, the first lipid
    mol2 -- A Molecule object, the second lipid
    cutoff -- A float, how close the two lipids must be to constitute a "clash"
    num_sub_partitions -- An integer, the number of partitions into which the atoms of each lipid molecule are divided. Clashes are then determined pairwise on each partition, rather than comparing every atom of one lipid to every atom of the other. Important for large system to avoid memory problems, but it can usually just be set to 1.
    params -- A dictionary, the user-specified comand-line parameters
    very_large_distance_check -- An optional Boolean, to enable the option of eliminating steric clashes early by examining the distance between the first atoms of each lipid

    Returns:
    A Boolean, True if the two lipids clash, False otherwise.
    
    """
    
    # of the user provided a Molecule object, use just the coordinate numpy array
    if not type(mol1) is numpy.ndarray: mol1 = mol1.all_atoms_numpy
    if not type(mol2) is numpy.ndarray: mol2 = mol2.all_atoms_numpy
    
    # first, check if the bounding boxes of each lipid couldn't possiblely overlap
    margin = numpy.array([cutoff, cutoff, cutoff])
    mol1_min = numpy.min(mol1,0) - margin
    mol2_max = numpy.max(mol2,0) + margin
    
    if mol1_min[0] > mol2_max[0]: return False
    if mol1_min[1] > mol2_max[1]: return False
    if mol1_min[2] > mol2_max[2]: return False
    
    mol1_max = numpy.max(mol1,0) + margin
    mol2_min = numpy.min(mol2,0) - margin

    if mol2_min[0] > mol1_max[0]: return False
    if mol2_min[1] > mol1_max[1]: return False
    if mol2_min[2] > mol1_max[2]: return False
    
    # now check if the distance between their first atoms is so great, they are unlikely to overlap
    if very_large_distance_check == True and numpy.linalg.norm(mol1[0] - mol2[0]) > params['very_distant_lipids_cutoff']: return False
    
    # now reduce the points to the ones that really need to be compared. Basically, stripping the points
    # I think this just retains points in the overlapping regions of the bounding boxes
    reduced_boundary_max = numpy.min(numpy.vstack((mol1_max, mol2_max)),0)
    reduced_boundary_min = numpy.max(numpy.vstack((mol1_min, mol2_min)),0)
    
    mol1_to_remove = (mol1[:,0] < reduced_boundary_max[0])
    if not True in mol1_to_remove: return False
    mol1_to_remove = numpy.logical_and(mol1_to_remove, (mol1[:,0] > reduced_boundary_min[0]))
    if not True in mol1_to_remove: return False
    mol1_to_remove = numpy.logical_and(mol1_to_remove, (mol1[:,1] < reduced_boundary_max[1]))
    if not True in mol1_to_remove: return False
    mol1_to_remove = numpy.logical_and(mol1_to_remove, (mol1[:,1] > reduced_boundary_min[1]))
    if not True in mol1_to_remove: return False
    mol1_to_remove = numpy.logical_and(mol1_to_remove, (mol1[:,2] < reduced_boundary_max[2]))
    if not True in mol1_to_remove: return False
    mol1_to_remove = numpy.logical_and(mol1_to_remove, (mol1[:,2] > reduced_boundary_min[2]))
    if not True in mol1_to_remove: return False
    
    mol2_to_remove = (mol2[:,0] < reduced_boundary_max[0])
    if not True in mol2_to_remove: return False
    mol2_to_remove = numpy.logical_and(mol2_to_remove, (mol2[:,0] > reduced_boundary_min[0]))
    if not True in mol2_to_remove: return False
    mol2_to_remove = numpy.logical_and(mol2_to_remove, (mol2[:,1] < reduced_boundary_max[1]))
    if not True in mol2_to_remove: return False
    mol2_to_remove = numpy.logical_and(mol2_to_remove, (mol2[:,1] > reduced_boundary_min[1]))
    if not True in mol2_to_remove: return False
    mol2_to_remove = numpy.logical_and(mol2_to_remove, (mol2[:,2] < reduced_boundary_max[2]))
    if not True in mol2_to_remove: return False
    mol2_to_remove = numpy.logical_and(mol2_to_remove, (mol2[:,2] > reduced_boundary_min[2]))
    if not True in mol2_to_remove: return False
    
    mol1_reduced = mol1[numpy.nonzero(mol1_to_remove)[0]]
    mol2_reduced = mol2[numpy.nonzero(mol2_to_remove)[0]]
    
    # now do a pairwise distance comparison between the remaining atoms
    if num_sub_partitions == 1:
        if True in (cdist(mol1_reduced, mol2_reduced) < cutoff): return True
    else:
        a1s = numpy.array_split(mol1_reduced, num_sub_partitions) # in my benchmarks, 35 gave good results
        a2s = numpy.array_split(mol2_reduced, num_sub_partitions)
        
        for mol1_reduced in a1s:
            for mol2_reduced in a2s:
                if len(mol1_reduced) > 0 and len(mol2_reduced) > 0:
                    if True in (cdist(mol1_reduced, mol2_reduced) < cutoff):
                        return True
                
    return False

def indices_of_close_pts(points1, points2, cutoff, num_sub_partitions): # this is faster, potentially, and less memory intensive
    """Examine two sets of points and return the indices of the points that are close to each other

    Arguments:
    points1 -- A nx3 numpy array, a set of points to be considered
    points2 -- A nx3 numpy array, a second set of points to be considered
    cutoff -- A float, how close the two lipids must be to constitute a "clash"
    num_sub_partitions -- An integer, the number of partitions into which the atoms of each lipid molecule are divided. Clashes are then determined pairwise on each partition, rather than comparing every atom of one lipid to every atom of the other. Important for large system to avoid memory problems, but it can usually just be set to 1.

    Returns:
    A tuple, containing a numpy array with indices from the first molecule and a numpy array with indices from the second molecule
    
    """
    
    if num_sub_partitions == 1:
        dists = cdist(points1, points2) < cutoff # which ones clash
        return numpy.nonzero(dists) # these are indices that clash
    else:
        a1s = numpy.array_split(points1, num_sub_partitions) # in my benchmarks, 35 gave good results
        a2s = numpy.array_split(points2, num_sub_partitions)
        
        points1_indices = []
        points2_indices = []
        
        a1s_index = 0
        for points1 in a1s:
            a2s_index = 0
            for points2 in a2s:
                if len(points1) > 0 and len(points2) > 0:
                    dists = cdist(points1, points2) < cutoff
                    indices = numpy.nonzero(dists)
                    points1_indices.extend(indices[0] + a1s_index)
                    points2_indices.extend(indices[1] + a2s_index)
                a2s_index = a2s_index + len(points2)
            a1s_index = a1s_index + len(points1)
        
        points1_indices = numpy.array([points1_indices])
        points2_indices = numpy.array([points2_indices])
        
        return (results[:,0], results[:,1])

################## FILL HOLES IN THE BILAYER ##################

def fill_in_lipid_holes(molecules_by_triangle, params):
    """Position lipid molecules in the bilayer holes

    Arguments:
    molecules_by_triangle -- A list of tuples, where each tuple contains a Triangle object and a list of lipid molecules (Molecule objects) that belong to that triangle
    params -- A dictionary, the user-specified command-line parameters

    Returns:
    A list of tuples, where each tuple contains a Triangle object, a list of lipid molecules (Molecule objects) that belong to that triangle, and an integer representing the index in the molecules_by_triangle list

    """

    headgroup_locs = [] # this does need to be a list rather than presized numpy arrays because its hard to know a priori how big it needs to be
    
    index_to_try = 0
    while len(headgroup_locs) < 5: # I feel like there needs to be at least 5 to get any kind of representative sample. Because the first triangel could conceivably not have any headgroups in it or only have one, keep looking until you find one that does
        # I'm not going to search through all of them to find the max because that might require loading a lot of pickles
        headgroup_locs = []

        # the first step is to get the average and minimum distance between headgroups. 
        if params['use_disk_instead_of_memory'] == "TRUE": lipids = load_pickle(molecules_by_triangle[0][1], params)
        else: lipids = molecules_by_triangle[index_to_try][1] # Just look one of the triangles as representative
        
        gc.disable()
        for alipid in lipids: headgroup_locs.append(alipid.all_atoms_numpy[alipid.get_headgroup_index(params['lipid_headgroup_marker'])])
        gc.enable()
        
        index_to_try = index_to_try + 1
    
    headgroup_locs = numpy.vstack(headgroup_locs)
    headgroup_dists = scipy.spatial.distance.squareform(pdist(headgroup_locs))
    headgroup_min_dists = numpy.empty(len(headgroup_dists))

    for indx in range(len(headgroup_dists)):
        t = headgroup_dists[indx]
        headgroup_min_dists[indx] = numpy.min(t[numpy.nonzero(t)])
    
    average_dist_between_headgroups = int(numpy.round(numpy.average(headgroup_min_dists)))
    min_dist_between_headgroups = numpy.min(headgroup_min_dists)
    pt_step = max([1, int(average_dist_between_headgroups/3.0)]) # so grid points every third of the way between headgroups
    
    # now, determine which triangles are adjacent
    adjacent_triangles_map = {}
    gc.disable()
    for index1 in range(len(molecules_by_triangle)-1):
        triangle_pts1 = molecules_by_triangle[index1][0]
        for index2 in range(index1 + 1, len(molecules_by_triangle)):
            triangle_pts2 = molecules_by_triangle[index2][0]
            
            if triangle_pts1.near_other_triangle(triangle_pts2, params):
                
                try: adjacent_triangles_map[index1].append(index2)
                except:
                    adjacent_triangles_map[index1] = []
                    adjacent_triangles_map[index1].append(index2)
                    
                try: adjacent_triangles_map[index2].append(index1)
                except:
                    adjacent_triangles_map[index2] = []
                    adjacent_triangles_map[index2].append(index1)

    class lipid_inserts_multiprocessing(general_task):
        """A class for inserting lipid molecules into bilayer holes"""
        
        def value_func(self, item, results_queue): # so overwriting this function
            """Insert lipid molecules into bilayer holes

            Arguments:
            item -- A list or tuple, the input data required for the calculation
            results_queue -- A multiprocessing.Queue() object for storing the calculation output
            
            """

            molecules_by_triangle_index = item[0]
            triangle_pts = item[1]
            lipids = item[2]
            adjacent_lipids = item[3] # molecules of lipids in neighboring triangles, but NOT in this one (i.e., a triangle is not adjacent to itself)
            params = item[4]
            average_dist_between_headgroups = item[5]
            min_dist_between_headgroups = item[6]
            pt_step = item[7]
            
            self.print_star_if_appropriate(molecules_by_triangle_index)
            
            if params['use_disk_instead_of_memory'] == "TRUE": lipids = load_pickle(lipids, params)
            
            ########## GET THE PLANE GOING THROUGH THE TRIANGLE PONITS #########
            
            # now get the plane going between these three points
            
            # the order of the triangle points could potentially matter in the case of
            # right triangles. So we potentially need to make sure every order is considered,
            # though we can abort early if an acceptable solution is found.
            # basically, in the case of right triangles, the point opposite the hypotenuse
            # needs to be projected onto the hypotenuse. With other kinds of triangles,
            # it can really be any point projected onto the opposite side.
            
            combos = []
            combos.append((triangle_pts[0], triangle_pts[1], triangle_pts[2]))
            combos.append((triangle_pts[0], triangle_pts[2], triangle_pts[1]))
            combos.append((triangle_pts[1], triangle_pts[2], triangle_pts[0]))
            
            for combo in combos:
            
                pt1 = combo[0]
                pt2 = combo[1]
                pt3 = combo[2]

                # project pt3 onto the line segment pt1-pt2
                u = pt1 - pt2 
                v = pt1 - pt3
                u = u/numpy.linalg.norm(u)
                new_pt = pt1 - numpy.dot(u,v) * u # this is the projected point
                
                # make sure the project point isn't equal to one of the triangle verticies
                if not numpy.array_equal(pt3, new_pt) and not numpy.array_equal(pt1, new_pt): break

            vec1 = pt3 - new_pt
            vec2 = pt1 - new_pt

            vec1 = vec1/numpy.linalg.norm(vec1) # two perpenticular vectors in the plane
            vec2 = vec2/numpy.linalg.norm(vec2) # and a point in the plane
                    
            plane_normal = numpy.cross(vec1, vec2) # a normal to the plane
            plane_normal = plane_normal/numpy.linalg.norm(plane_normal)
            
            # good to get a scalar equation for the plane too: ax + by + cz + d = 0
            scalar_eq_a = plane_normal[0]
            scalar_eq_b = plane_normal[1]
            scalar_eq_c = plane_normal[2]
            scalar_eq_d = -numpy.dot(triangle_pts.center(), plane_normal)
            
            # now that the plane has been identified, find the average distance between the plane and lipid headgroups
            # also, start adding lipids that could clash with future inserted lipids into the neighborhood_lipids_that_could_clash list. All lipids in the margin and submargin of the central triangle will be added.
            lipid_head_indices = numpy.empty(len(lipids))
            for indx, lipid in enumerate(lipids):
                lipid_head_indices[indx] = lipid.get_headgroup_index(params['lipid_headgroup_marker'])
            
            all_lipid_heads_loc_in_central_triangle = numpy.empty((len(lipid_head_indices), 3))
            neighborhood_lipids_that_could_clash = [] 
            headgroup_locs_of_lipids_that_could_clash = []
            for t in range(len(lipid_head_indices)):
                all_lipid_heads_loc_in_central_triangle[t] = lipids[t].all_atoms_numpy[lipid_head_indices[t]]
                
                if lipids[t].in_triangle_margin == True or lipids[t].in_triangle_submargin == True:
                    neighborhood_lipids_that_could_clash.append(lipids[t])
                    headgroup_locs_of_lipids_that_could_clash.append(lipids[t].all_atoms_numpy[lipid_head_indices[t]])
            
            three_scalars = numpy.array([scalar_eq_a, scalar_eq_b, scalar_eq_c])
            dists2 = numpy.empty(len(all_lipid_heads_loc_in_central_triangle))
            for indx in range(len(all_lipid_heads_loc_in_central_triangle)):
                lipid_head_pt = all_lipid_heads_loc_in_central_triangle[indx]
                dist = numpy.fabs(numpy.dot(three_scalars, lipid_head_pt) + scalar_eq_d) / numpy.power(numpy.dot(three_scalars, three_scalars), 0.5)
                dists2[indx] = dist
            
            if len(dists2) == 0: # if there are no lipid headgroups in this triangle, so you can't proceed
                positioned_molecules = []
                if params['use_disk_instead_of_memory'] == "TRUE": self.results.append((molecules_by_triangle_index, save_pickle(positioned_molecules, params))) # here save the results for later compilation
                else: self.results.append((molecules_by_triangle_index, positioned_molecules)) # here save the results for later compilation
                return
            
            average_headgroup_dist_to_plane = numpy.mean(dists2)

            # Find the locations of the in-margin headgroups of all adjacent triangles
            gc.disable()
            for lipids2 in adjacent_lipids: # note that this does NOT include the lipids in the central triangle, which were identified above.
                if params['use_disk_instead_of_memory'] == "TRUE": lipids2 = load_pickle(lipids2, params)
                
                for alipid in lipids2:
                    
                    if alipid.in_triangle_margin == True: # so for neighboring triangles, we only care about the lipids that are in the margin, which might clash with future inserted lipids
                        neighborhood_lipids_that_could_clash.append(alipid)
                        headgroup_locs_of_lipids_that_could_clash.append(alipid.all_atoms_numpy[alipid.get_headgroup_index(params['lipid_headgroup_marker'])])
            gc.enable()
            
            # need to numpify headgroup_locs_of_lipids_that_could_clash
            headgroup_locs_of_lipids_that_could_clash = numpy.array(headgroup_locs_of_lipids_that_could_clash)
            
            # now flood the surface of both bilayers with points
            # first, generate a field of points
            s = numpy.arange(-triangle_pts.max_radius(), triangle_pts.max_radius(), pt_step, dtype=int)
            t = numpy.arange(-triangle_pts.max_radius(), triangle_pts.max_radius(), pt_step, dtype=int)
            pts = numpy.empty((len(s)*len(t),3))
            for s_index, s_val in enumerate(s): 
                for t_index, t_val in enumerate(t):
                    pt = s_val * vec1 + t_val * vec2
                    pts[s_index * len(t) + t_index][0] = pt[0]
                    pts[s_index * len(t) + t_index][1] = pt[1]
                    pts[s_index * len(t) + t_index][2] = pt[2]
            pts = numpy.array(pts) + triangle_pts.center()

            # check which of these points are within the central triangle
            indices_of_pts_in_triangle = triangle_pts.get_indices_of_points_within_triangle_boundaries(pts)
            pts_in_triangle = get_numpy_slice(pts,indices_of_pts_in_triangle)
            
            # now remove points that are too far in the interior. Fill only at triangle edges
            smaller_tri = triangle_pts.new_triangle_expanded_by_margin(-params['clashing_potential_margin'])
            indices_of_pts_in_triangle = smaller_tri.get_indices_of_points_within_triangle_boundaries(pts_in_triangle)
            pts_in_triangle = numpy.delete(pts_in_triangle, indices_of_pts_in_triangle,axis=0)
            
            # create points above and below each of these grid points
            local_pts_to_examine = numpy.empty((2*len(pts_in_triangle),3))
            for apt_index, apt in enumerate(pts_in_triangle):
                # now get the two points above and below the plane
                starting_pts = numpy.array([apt - plane_normal * average_headgroup_dist_to_plane, apt + plane_normal * average_headgroup_dist_to_plane])

                # place those two points into the local_pts_to_examine numpy array
                for starting_pt_index, starting_pt in enumerate(starting_pts):
                    local_pts_to_examine[apt_index * 2 + starting_pt_index][0] = starting_pt[0]
                    local_pts_to_examine[apt_index * 2 + starting_pt_index][1] = starting_pt[1]
                    local_pts_to_examine[apt_index * 2 + starting_pt_index][2] = starting_pt[2]

            # remove all pts that are too close to the headgroups
            indices_of_clashing_pts = indices_of_close_pts(headgroup_locs_of_lipids_that_could_clash, local_pts_to_examine, average_dist_between_headgroups, params['memory_optimization_factor'])[1]
            local_pts_to_examine = numpy.delete(local_pts_to_examine, indices_of_clashing_pts, 0)
            
            # remove all remaining pts that clash with other lipid atoms (headgroups first to reduce number of pair-wise distance comparisons)
            for lip in neighborhood_lipids_that_could_clash:
                indices_of_clashing_pts = indices_of_close_pts(lip.all_atoms_numpy, local_pts_to_examine, min_dist_between_headgroups, params['memory_optimization_factor'])[1]
                
                # indices_of_clashing_pts could be empty, so just try
                try: local_pts_to_examine = numpy.delete(local_pts_to_examine, indices_of_clashing_pts, 0)
                except: pass

            # now position lipids 
            positioned_molecules = [] # can't know size, so can't preallocate in numpy array
            positioned_molecules_headgroup_locs = [] # can't know size, so can't preallocate
            gc.disable() # because appending complex objects to a list
            
            for t in range(params['fill_hole_exhaustiveness']):
                indxs = range(len(local_pts_to_examine))
                random.shuffle(indxs) # so not examining points sequentially
                for headgroup_loc_index in indxs:
                    if headgroup_loc_index < len(local_pts_to_examine): # the point could have been deleted, in which case you should skip
                        
                        new_head_group_loc = local_pts_to_examine[headgroup_loc_index]
            
                        # determine the directionality of the lipid (i.e., points "up" or "down")
                        candidates_pts = numpy.array([new_head_group_loc - plane_normal, new_head_group_loc + plane_normal])
                        dists_to_center = cdist(candidates_pts, numpy.array([triangle_pts.center()]))
                        
                        if dists_to_center[0] < dists_to_center[1]: directionality = 1
                        else: directionality = -1
                        
                        # pick a random lipid
                        lipid = random.choice(lipids) # maybe needs to be a copy?
                        lipid_head_loc_index = lipid.get_headgroup_index(params['lipid_headgroup_marker'])
                        lipid_head_loc = lipid.all_atoms_numpy[lipid_head_loc_index]
                        lipid_center_loc = numpy.mean(lipid.all_atoms_numpy, 0)
                        lipid_length = numpy.linalg.norm(lipid_head_loc - lipid_center_loc)
                        
                        # you should be working with a copy of the lipid, not the original
                        lipid = lipid.copy_of()
                        lipid.in_triangle_margin = True
                        lipid.in_triangle_submargin = False
                        
                        # get new guide (static) template. This specifies where the lipid will ultimately be moved to
                        lipid_center_guidepoint = new_head_group_loc - directionality * lipid_length * plane_normal
                        guide_static_template = numpy.array([new_head_group_loc, lipid_center_guidepoint])
                    
                        # get new dynamic template. this is the starting location of the lipid before it's moved to the new location.
                        dynamic_template = numpy.array([lipid_head_loc, lipid_center_loc])
                        
                        # get origin template. This is a destination at the origin. You'll move it here for rotating before moving it to the new location
                        origin = numpy.array([0.0, 0.0, 0.0])
                        origin2 = origin - lipid_length * numpy.array([0.0, 0.0, 1.0])
                        guide_origin_template = numpy.array([origin, origin2])
            
                        # move lipid to origin.
                        transform_data = get_transformation_data(guide_origin_template, dynamic_template)
                        apply_transformation(lipid, transform_data)
                        
                        # now rotate about z axis
                        theta = random.random() * numpy.pi * 2.0
                        rot_max = numpy.array([
                            [numpy.cos(theta), -numpy.sin(theta), 0.0],
                            [numpy.sin(theta), numpy.cos(theta), 0.0],
                            [0.0, 0.0, 1.0]
                        ])
                        lipid.all_atoms_numpy = numpy.dot(lipid.all_atoms_numpy, rot_max)
                        
                        # now move to correct location in bilayer
                        center_dynamic_pdb, rot_quat, center_static_pdb = get_transformation_data(guide_static_template, guide_origin_template)
                        lipid.all_atoms_numpy = lipid.all_atoms_numpy - center_dynamic_pdb
                        lipid.rotate_mol_quat(rot_quat)
                        lipid.all_atoms_numpy = lipid.all_atoms_numpy + center_static_pdb
                        
                        # redefine the lead group location now that things have been moved    
                        lipid_head_loc = lipid.all_atoms_numpy[lipid_head_loc_index]
                        
                        # check to see if the positioned lipid clashes with other lipids
                        some_clash = False
                        first_pt_dists = cdist(headgroup_locs_of_lipids_that_could_clash, numpy.array([lipid_head_loc]))
                        first_pt_close_indices = numpy.nonzero(first_pt_dists < params['very_distant_lipids_cutoff'])[0]
                        for indx in first_pt_close_indices:
                            if two_lipids_clash(lipid, neighborhood_lipids_that_could_clash[indx], params['clash_cutoff'], 1, params, False) == True:
                                some_clash = True
                                break
                            
                        if some_clash == False: 
                            
                            if len(positioned_molecules_headgroup_locs) > 0:
                                
                                positioned_pt_dists = cdist(numpy.array(positioned_molecules_headgroup_locs), numpy.array([lipid_head_loc]))
                                positioned_pt_close_indices = numpy.nonzero(positioned_pt_dists < params['very_distant_lipids_cutoff'])[0]
                                for indx in positioned_pt_close_indices:
                                    if two_lipids_clash(lipid, positioned_molecules[indx], params['clash_cutoff'], 1, params, False) == True:
                                        some_clash = True
                                        break
                        
                        if some_clash == False: # so it doesn't clash. save it.
                            
                            positioned_molecules.append(lipid) # remember, a copy
                            positioned_molecules_headgroup_locs.append(lipid_head_loc)
                            
                            # now remove surface points from local_pts_to_examine that come close to the newly positioned lipid
                            dists = cdist(lipid.all_atoms_numpy, local_pts_to_examine) < min_dist_between_headgroups # which ones clash
                            indices_of_clashing_pts = numpy.nonzero(dists)[1] # these are indices that clash
                            local_pts_to_examine = numpy.delete(local_pts_to_examine, indices_of_clashing_pts, 0)
                            
            # now add all these positioned lipids to the molecules_by_triangle list
            if params['use_disk_instead_of_memory'] == "TRUE":
                self.results.append((molecules_by_triangle_index, save_pickle(positioned_molecules, params))) # here save the results for later compilation
            else: self.results.append((molecules_by_triangle_index, positioned_molecules)) # here save the results for later compilation
            gc.enable()
            
    # fill the lipid holes using multiple processors if possible
    some_input = []
    for molecules_by_triangle_index in range(len(molecules_by_triangle)):
        triangle_pts = molecules_by_triangle[molecules_by_triangle_index][0]
        lipids = molecules_by_triangle[molecules_by_triangle_index][1]
        adjacent_lipids = [molecules_by_triangle[index][1] for index in adjacent_triangles_map[molecules_by_triangle_index]]
        some_input.append((molecules_by_triangle_index, triangle_pts, lipids, adjacent_lipids, params, average_dist_between_headgroups, min_dist_between_headgroups, pt_step))
    
    gc.enable()

    positioned_lipids = multi_threading(some_input, params['number_of_processors'], lipid_inserts_multiprocessing, params, "REMARK ").results

    # now organize the positioned_lipids into the same organization as molecules_by_triangle for subsequent processing
    positioned_lipids_by_triangle = []
    gc.disable()
    for molecules_by_triangle_index, positioned_molecules in positioned_lipids:
        positioned_lipids_by_triangle.append((molecules_by_triangle[molecules_by_triangle_index][0], positioned_molecules, molecules_by_triangle_index))
    gc.enable()

    return positioned_lipids_by_triangle

################## HANDLE SPECIAL OUTPUT FILES ##################

def print_out_mesh_points(all_triangles, params):
    """Save the mesh points to a PDB file

    Arguments:
    all_triangles -- A list of Triangle objects, the tesselation/triangulation
    params -- A dictionary, the user-specified command-line parameters

    """

    def create_pdb_line(numpy_array, letter):
       
        """Create a string formatted according to the PDB standard from the atomic information contained in this atom class.
       
        Arguments:
        numpy_array -- A numpy array, containing the atom coordinates
        letter -- A string, which will serve as the atom name, residue name, chain, etc.
       
        Returns:
        A string, formatted according to the PDB standard.
       
        """
       
        output = "ATOM "
        output = output + "0".rjust(6) + letter.rjust(5) + "0".rjust(4) + letter.rjust(2) + "0".rjust(4)
        output = output + ("%.3f" % numpy_array[0]).rjust(12)
        output = output + ("%.3f" % numpy_array[1]).rjust(8)
        output = output + ("%.3f" % numpy_array[2]).rjust(8)
        output = output + letter.rjust(24) # + "   " + str(uniqueID) #This last part must be removed
        return output
    
    toprint = []
    point_already_shown = []
    for tile_triangle in all_triangles:
        key = str(tile_triangle[0][0]) + "_" + str(tile_triangle[0][1]) + "_" + str(tile_triangle[0][2])
        if not key in point_already_shown:
            point_already_shown.append(key)
            toprint.append(create_pdb_line(tile_triangle[0], "X"))
    
        key = str(tile_triangle[1][0]) + "_" + str(tile_triangle[1][1]) + "_" + str(tile_triangle[1][2])
        if not key in point_already_shown:
            point_already_shown.append(key)
            toprint.append(create_pdb_line(tile_triangle[1], "X"))
    
        key = str(tile_triangle[2][0]) + "_" + str(tile_triangle[2][1]) + "_" + str(tile_triangle[2][2])
        if not key in point_already_shown:
            point_already_shown.append(key)
            toprint.append(create_pdb_line(tile_triangle[2], "X"))
    if params['output_directory'] == "": print "\n".join(toprint)
    else:
        f = open(params['output_directory'] + 'grid_points.pdb', 'w')
        f.write("\n".join(toprint))
        f.close()

def print_out_triangle_tcl_file(all_triangles, params): # reimplement the above similarly. This code is repeated.
    """Save the tesselation/triangulation to a TCL file

    Arguments:
    all_triangles -- A list of Triangle objects, the tesselation/triangulation
    params -- A dictionary, the user-specified command-line parameters

    """

    # draw triangles
    f = open(params['output_directory'] + 'triangles.tcl','w')
    f.write("draw delete all\n")
    f.write("draw color red\n")
    
    for triangle in all_triangles:
        ia_pt = triangle[0]
        ib_pt = triangle[1]
        ic_pt = triangle[2]
        
        f.write("draw triangle {" + str(ia_pt[0]) + " " + str(ia_pt[1]) + " " + str(ia_pt[2]) + "} {" + str(ib_pt[0]) + " " + str(ib_pt[1]) + " " + str(ib_pt[2]) + "} {" + str(ic_pt[0]) + " " + str(ic_pt[1]) + " " + str(ic_pt[2]) + "}\n")
    
    f.close()

################## RUN THE PROGRAM ##################

def run_program(argv):
    starttime = time.time() # to keep track of execution time
    
    current_step = 0 # used for output filenames
    
    print "\nREMARK      LipidWrapper " + version + "\n"
    
    # check for Tkinter
    try:
        import Tkinter
        print "REMARK      The Tkinter python module is available. You may prefer to"
        print "REMARK      use the LipidWrapper graphical user interface"
        print "REMARK      (lipidwrapperGUI.py).\n"
        del Tkinter
    except: pass # GUI not available
    
    params = get_commandline_parameters(argv) # get the commandline parameters
    
    # if you're going to be storing the growing model on the disk, make the temporary directory
    if params['use_disk_instead_of_memory'] == "TRUE":
        if os.path.exists(params['memory_store_dir']): shutil.rmtree(params['memory_store_dir'])
        os.mkdir(params['memory_store_dir'])
    
    # load mesh points and generate triangulation
    print "REMARK      Loading/creating and triangulating the mesh..."
    all_triangles = load_mesh_points_and_triangulations(params) # get the triangulations
    
    # load in the user-specified planar bilayer model
    print "REMARK      Loading the original lipid-bilayer model (" + params['lipid_pdb_filename'] + ")..."
    lipid, min_headgroups, max_headgroups = load_lipid_model(params) # get the lipid molecule object, properly centered on x-y plane, as well as bounding-box coordinates
    
    # fill the tessellated triangles with bilayer sections
    print "REMARK      Position copies of the lipid bilayer on the trianguled mesh..."
    molecules_by_triangle = position_lipid_model_on_triangulated_tiles(params, lipid, all_triangles, min_headgroups, max_headgroups) # position the lipids on the triangles
    
    # save the bilayer sections if user requested
    if params['output_directory'] != "": 
        print "REMARK      Saving positioned lipid bilayers..."
        current_step = current_step + 1
        
        dir_pathname = params['output_directory'] + "step_" + str(current_step) + '.positioned_lipid_triangles' + "/"
        if not os.path.exists(dir_pathname): os.mkdir(dir_pathname)
    
        class save_positioned_lipids_multiprocessing(general_task):
            """A class for saving the lipid molecules associated with each triangle"""
            
            def value_func(self, item, results_queue): # so overwriting this function
                """Save lipid molecules associated with a triangle
    
                Arguments:
                item -- A list or tuple, the input data required for the calculation
                results_queue -- A multiprocessing.Queue() object for storing the calculation output
                
                """

                lipids = item[0]
                i = item[1]
                current_step = item[2]
                params = item[3]
                dir_pathname = item[4]
                
                self.print_star_if_appropriate(i)
        
                f = openfile(dir_pathname + 'step_' + str(current_step) + '.original_positioned_lipid_triangle.' + str(i) + ".pdb", 'w', params)
                
                if params['use_disk_instead_of_memory'] == "TRUE": lipids = load_pickle(lipids, params)
                
                for lipid in lipids:
                    for index in range(len(lipid.all_atoms_numpy)): f.write(lipid.create_pdb_line(index) + "\n")
                f.close()
                
        some_input = []
        i = 1
        gc.disable()
        for discard, lipids in molecules_by_triangle:
            i = i + 1
            some_input.append((lipids, i, current_step, params, dir_pathname))
        gc.enable()
            
        multi_threading(some_input, params['number_of_processors'], save_positioned_lipids_multiprocessing, params, "REMARK ")
    
    # delete clashing lipids if user requested
    if params['delete_clashing_lipids'] == "TRUE":
        print "REMARK      Deleting lipids that sterically clash..."
        remove_steric_clashes(molecules_by_triangle, params) # remove steric clashes between lipids of adjacent tiles
        
        # save work from this step if user requested
        if params['output_directory'] != "": 
            current_step = current_step + 1
            dir_pathname = params['output_directory'] + "step_" + str(current_step) + '.remove_lipids_with_clashes' + "/"
            if not os.path.exists(dir_pathname): os.mkdir(dir_pathname)
    
            # print out remaining lipids
            print "REMARK            Saving the lipids that were not deleted for reference..."
    
            class save_nondeleted_lipids_multiprocessing(general_task):
                """A class for saving the lipid molecules that were not deleted"""
            
                def value_func(self, item, results_queue): # so overwriting this function
                    """Save lipid molecules that were not deleted
        
                    Arguments:
                    item -- A list or tuple, the input data required for the calculation
                    results_queue -- A multiprocessing.Queue() object for storing the calculation output
                    
                    """

                    triangle_index = item[0]
                    dir_pathname = item[1]
                    current_step = item[2]
                    lipids = item[3]
                    
                    self.print_star_if_appropriate(triangle_index)
            
                    f = openfile(dir_pathname + 'step_' + str(current_step) + ".retained_lipids_no_clash." + str(triangle_index + 1) + ".pdb", 'w', params)
                    if params['use_disk_instead_of_memory'] == "TRUE": triangle_lipids = load_pickle(lipids, params)
                    else: triangle_lipids = lipids
                    
                    for lipid in triangle_lipids:
                        for index in range(len(lipid.all_atoms_numpy)):
                            f.write(lipid.create_pdb_line(index) + "\n")
                    f.close()
    
            some_input = []
            gc.disable()
            for triangle_index in range(len(molecules_by_triangle)): some_input.append((triangle_index, dir_pathname, current_step, molecules_by_triangle[triangle_index][1]))
            gc.enable()
            multi_threading(some_input, params['number_of_processors'], save_nondeleted_lipids_multiprocessing, params, "REMARK ")
            
        # fill holes in bilayer if user requested
        if params['fill_holes'] == 'TRUE':

            # fill the holes
            print "REMARK      Filling holes in the bilayer with additional lipid molecules..."
            positioned_lipids_by_triangle = fill_in_lipid_holes(molecules_by_triangle, params)
            
            # remove filling lipids that clash
            print "REMARK            Removing added lipids that clash..."
            remove_steric_clashes(positioned_lipids_by_triangle, params)
            
            # save work from this step if user requested
            if params['output_directory'] != "":
                print "REMARK            Saving the lipids that were added for reference..."
                current_step = current_step + 1
                dir_pathname = params['output_directory'] + "step_" + str(current_step) + '.lipid_holes_plugged' + "/"
                if not os.path.exists(dir_pathname): os.mkdir(dir_pathname)
                
                class save_plugging_lipids_multiprocessing(general_task):
                    """A class for saving the lipid molecules that were placed in lipid holes"""
                    
                    def value_func(self, item, results_queue): # so overwriting this function
                        """Save lipid molecules that were placed in lipid holes
            
                        Arguments:
                        item -- A list or tuple, the input data required for the calculation
                        results_queue -- A multiprocessing.Queue() object for storing the calculation output
                        
                        """

                        triangle_lipids = item[0]
                        index = item[1]
                        dir_pathname = item[2]
                        current_step = item[3]
                        params = item[4]
                        
                        self.print_star_if_appropriate(index)
                        
                        if params['use_disk_instead_of_memory'] == "TRUE": triangle_lipids = load_pickle(triangle_lipids, params)
                        
                        f = openfile(dir_pathname + 'step_' + str(current_step) + ".lipids_added_into_bilayer_holes." + str(index + 1) + ".pdb", 'w', params)
                        for lipid in triangle_lipids:
                            for i in range(len(lipid.all_atoms_numpy)): f.write(lipid.create_pdb_line(i) + "\n")
                        f.close()
    
                some_input = []
                gc.disable()
                for ignore_var, lipids, index in positioned_lipids_by_triangle:
                    some_input.append((lipids, index, dir_pathname, current_step, params))
                gc.enable()
                
                multi_threading(some_input, params['number_of_processors'], save_plugging_lipids_multiprocessing, params, "REMARK ")
    
            # now add the positioned ligands into the lipid list associated with the original triangle
            print "REMARK            Adding the hole-filling lipids to the original models..."
            if params['use_disk_instead_of_memory'] == "TRUE":
                
                class add_plugging_lipids_multiprocessing(general_task):
                    """A class for adding the lipid molecules used to plug lipid holes to the list of lipids associated with the relevant triangle"""
                    
                    def value_func(self, item, results_queue): # so overwriting this function
                        """Add lipid molecules used to plug holes to their associated triangles
            
                        Arguments:
                        item -- A list or tuple, the input data required for the calculation
                        results_queue -- A multiprocessing.Queue() object for storing the calculation output
                        
                        """

                        existing_lipids_pickle_id = item[0]
                        position_lipids_pickle_id = item[1]
                        params = item[3]
                        
                        self.print_star_if_appropriate(item[2])
                        
                        tmp = load_pickle(existing_lipids_pickle_id, params)
                        tmp2 = load_pickle(position_lipids_pickle_id, params)
                        tmp.extend(tmp2)
                        save_pickle(tmp, params, existing_lipids_pickle_id)
                
                some_input = []
                gc.disable()
                for var_not_needed, lipids, index in positioned_lipids_by_triangle: some_input.append((molecules_by_triangle[index][1], lipids, index, params))
                gc.enable()
                multi_threading(some_input, params['number_of_processors'], add_plugging_lipids_multiprocessing, params, "REMARK ")
    
            else: # your going to have to merge all the lists into the main one regardless, so I see no utility in using multiple processors
                for var_not_needed, lipids, index in positioned_lipids_by_triangle: molecules_by_triangle[index][1].extend(lipids)
    
            # save work from this step if user requested
            if params['output_directory'] != "":
                # print out all lipids
                print "REMARK            Saving the bilayers with holes plugged for reference..."
    
                class save_plugged_lipids_multiprocessing(general_task):
                    """A class for saving the lipid molecules used to plug lipid holes"""
    
                    def value_func(self, item, results_queue): # so overwriting this function
                        """Save lipid molecules used to plug lipid holes
            
                        Arguments:
                        item -- A list or tuple, the input data required for the calculation
                        results_queue -- A multiprocessing.Queue() object for storing the calculation output
                        
                        """

                        index = item[0]
                        dir_pathname = item[1]
                        current_step = item[2]
                        params = item[3]
                        lipids = item[4]
                        
                        self.print_star_if_appropriate(index)
                        
                        f = openfile(dir_pathname + 'step_' + str(current_step) + ".all_lipids_with_holes_plugged." + str(index + 1) + ".pdb", 'w', params)
                        
                        if params['use_disk_instead_of_memory'] == "TRUE": triangle_lipids = load_pickle(lipids, params)
                        else: triangle_lipids = lipids
                        
                        for lipid in triangle_lipids:
                            for i in range(len(lipid.all_atoms_numpy)): f.write(lipid.create_pdb_line(i) + "\n")
                        f.close()
    
                some_input = []
                gc.disable()
                for index in range(len(molecules_by_triangle)): some_input.append((index, dir_pathname, current_step, params, molecules_by_triangle[index][1]))
                gc.enable()
                multi_threading(some_input, params['number_of_processors'], save_plugged_lipids_multiprocessing, params, "REMARK ")
    
    # now print out all the final molecules
    print "REMARK      Printing out or saving all lipids to a single file..."
    if params['output_directory'] != "":
        current_step = current_step + 1
        dir_pathname = params['output_directory'] + "step_" + str(current_step) + '.final_lipid_triangles' + "/"
        if not os.path.exists(dir_pathname): os.mkdir(dir_pathname)
        f = openfile(params['output_directory'] + "step_" + str(current_step + 1) + ".full_bilayer.pdb", 'w', params)
    
        class save_final_lipids_multiprocessing(general_task):
            """A class for saving the final lipid models"""
            
            def value_func(self, item, results_queue): # so overwriting this function
                """Save the final lipid models
    
                Arguments:
                item -- A list or tuple, the input data required for the calculation
                results_queue -- A multiprocessing.Queue() object for storing the calculation output
                
                """

                params = item[0]
                lipids_lists = item[1]
                dir_pathname = item[2]
                current_step = item[3]
                triangle_index = item[4]
                
                self.print_star_if_appropriate(triangle_index)
    
                if params['use_disk_instead_of_memory'] == "TRUE": triangle_lipids = load_pickle(lipids_lists, params)
                else: triangle_lipids = lipids_lists
                
                g = openfile(dir_pathname + 'step_' + str(current_step) + '.final_lipid_triangle.' + str(triangle_index + 1) + ".pdb", 'w', params)
                
                for lipid in triangle_lipids:
                    for index in range(len(lipid.all_atoms_numpy)): g.write(lipid.create_pdb_line(index) + "\n")
                
                g.close()
        
        some_input = []
        gc.disable()
        for triangle_index in range(len(molecules_by_triangle)): some_input.append((params, molecules_by_triangle[triangle_index][1], dir_pathname, current_step, triangle_index))
        gc.enable()
        multi_threading(some_input, params['number_of_processors'], save_final_lipids_multiprocessing, params, "REMARK ")
    
    # print out single files
    atomindex = 0
    resindex = 0
    
    for triangle_index in range(len(molecules_by_triangle)):
    
        if params['use_disk_instead_of_memory'] == "TRUE": triangle_lipids = load_pickle(molecules_by_triangle[triangle_index][1], params)
        else: triangle_lipids = molecules_by_triangle[triangle_index][1]
        
        for lipid in triangle_lipids:
            resindex = resindex + 1
            for index in range(len(lipid.all_atoms_numpy)):
                atomindex = atomindex + 1
                if params['output_directory'] != "": f.write(lipid.create_pdb_line(index, atomindex, resindex) + "\n")
                else: print lipid.create_pdb_line(index, atomindex, resindex)
    
    if params['output_directory'] != "": f.close()
    
    # optional output files
    if params['show_grid_points'] == "TRUE":
        print "REMARK      Printing out or saving the grid points for reference..."
        print_out_mesh_points(all_triangles, params)
    if params['output_directory'] != "" or params['create_triangle_tcl_file'] == "TRUE":
        print "REMARK      Creating a VMD TCL file showing the triangulations..."
        print_out_triangle_tcl_file(all_triangles, params)
    
    # if the disk was used instead of memory, delete the temporary directory
    if params['use_disk_instead_of_memory'] == "TRUE": shutil.rmtree(params['memory_store_dir'])
    
    # tell the user how long it took for the program to execute
    print "REMARK      Execution time: " + str(time.time() - starttime) + " seconds"

if __name__=="__main__":
    run_program(sys.argv)
