#!/usr/bin/python
import sys, subprocess, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.Polypeptide import *
from Bio.PDB import *
#from Bio.PDB import PDBParser
#from Bio.PDB import PDBIO
#from Bio.PDB import Dice
from collections import OrderedDict
import math
import numpy

# This code reads in a given input PDB file, then calculates the coordinates of the Center-Of-Mass (COM) of the side chain (SC) atoms and the Center-Of-Mass (COM) of the entire Amino Acid (AA) atoms.
# It also outputs the atomic coordinates for N,C,CA,O, and CB atoms in the backbone and side chain.
# Then writes these coordinates in temporary files according to the required format of the input files to VORO++, then VORO++ will be prompted to calculate the Voronoi volume of individual amino acids of the requested PDB structures.
# Three coordinate sets will be used and fed to VORO++ sperately to represent amino acids in the PDB file.
#       1. the coordinates of CA atoms.
#       2. the coordinates of the Center-Of-Mass (COM) of the side chain atoms
#       3. the coordinates of the Center-Of-Mass (COM) of the entire Amino Acid atoms

# In the temporary files, on each row of the output, first appears the ID of the coordinate, then the X Y Z coordinates respectively.
# On the output, the results of VORO++ will written to the three filenames given by the user on the command line when running this code, each containing the VORO++ output for a set of coordinates.

# INPUT:  pdb files in ../structures/*  and the name of the output file.
# OUTPUT: Two files that will contain the output of Voro++ for the COM of Side Chain atoms of Amino Acids in the pdb file
# Amir Shahmoradi, 3:57 PM, Tuesday Aug 16 2014, Wilke Lab, iCMB, The University of Texas at Austin.

# This is a dictionary of residue volumes taken from the book "Protein Structure" Darby & Creighton. 
resvol_dict = { 'ALA':  67.0, 'ARG': 148.0, 'ASN':  96.0, 'ASP':  91.0, 'CYS':  86.0, \
                'GLN': 114.0, 'GLU': 109.0, 'GLY':  48.0, 'HIS': 118.0, 'ILE': 124.0, \
                'LEU': 124.0, 'LYS': 135.0, 'MET': 124.0, 'PHE': 135.0, 'PRO':  90.0, \
                'SER':  73.0, 'THR':  93.0, 'TRP': 163.0, 'TYR': 141.0, 'VAL': 105.0  }

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
res_dict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
             'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
             'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
             'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }


def main():
    if len( sys.argv ) != 9:
        print '''
 
Usage:'''
        print "     ", sys.argv[0], "<input PDB file>", "<Voronoi output file for SC coordinates>" , "<Voronoi output file for AA coordinates>" , "<Voronoi output file for AA coordinates>" , "<Voronoi output file for CA coordinates>" , "<Voronoi output file for CB coordinates>" , "<Voronoi output file for N coordinates>" , "<Voronoi output file for C coordinates>" , "<Voronoi output file for O coordinates>" , '\n'
        sys.exit('Program stopped.\n')
    else:
        pdb_in = sys.argv[1]  # path for the input PDB file to be read
        sum_SC = sys.argv[2]   # summary file containing all residue vornoi properties for SC coordinates
        sum_AA = sys.argv[3]   # summary file containing all residue vornoi properties for AA coordinates
        sum_CA = sys.argv[4]   # summary file containing all residue vornoi properties for CA coordinates
        sum_CB = sys.argv[5]   # summary file containing all residue vornoi properties for CB coordinates
        sum_N  = sys.argv[6]   # summary file containing all residue vornoi properties for N  coordinates
        sum_C  = sys.argv[7]   # summary file containing all residue vornoi properties for C  coordinates
        sum_O  = sys.argv[8]   # summary file containing all residue vornoi properties for O  coordinates

    p = PDBParser()
    pdb_name = pdb_in[-10:-4]
    pdb_chain = pdb_in[-5:-4]
    structure = p.get_structure(pdb_name,pdb_in)
    
    resnam   = []     # A list containing all Residue Names
    resnum   = []     # A list containing all Residue Numbers
    reschain = []     # A list containing the chain name in which the residues lie
    crdSC    = []     # Side Chain (SC) geometrical average coordinates
    crdAA    = []     # Amino Acid (AA) geometrical average coordinates
    crdCA    = []     # backbone atom CA coordinates
    crdCB    = []     # side chain atom CB coordinates
    crdN     = []     # backbone atom N coordinates
    crdC     = []     # backbone atom C coordinates
    crdO     = []     # backbone atom O coordinates
    sizeSC   = []   # A list containing the total number of atoms in each residue Side Chain (SC).
    sizeAA   = []   # A list containing the total number of atoms in each Amino Acid (AA).
    
    filenameSC = open('crdSC','w')
    filenameAA = open('crdAA','w')
    filenameCA = open('crdCA','w')
    filenameCB = open('crdCB','w')
    filenameN  = open('crdN' ,'w')
    filenameC  = open('crdC' ,'w')
    filenameO  = open('crdO' ,'w')
    pdb_size = 0
    for residue in structure.get_residues():
        pdb_size += 1
        resnam.append(residue.resname)
        resnum.append(residue.get_full_id()[3][1])
        reschain.append(residue.get_full_id()[2])
        rescrd_SC = []  # A list containing the coordinates of all side chain atoms of the current residue. Will be used to calculate the COM of the side chain.
        rescrd_AA = []  # A list containing the coordinates of all atoms of the current Amino Acid. Will be used to calculate the COM of the Amino Acid.
        noSC = True
        noCA = True
        noCB = True
        noN  = True
        noC  = True
        noO  = True

        for atom in structure.get_atoms():
            # atom.name is equivalent to atom.get_id()
            if atom.parent.id == residue.id and atom.name not in ['C','CA','O','N']:
                noSC = False
                rescrd_SC.append(atom.get_coord())
            
            if atom.parent.id == residue.id:
                rescrd_AA.append(atom.get_coord())

            if atom.parent.id == residue.id and atom.name == 'CA':
                noCA = False
                crdCA.append(atom.get_coord())
            elif atom.parent.id == residue.id and atom.name == 'CB':
                noCB = False
                crdCB.append(atom.get_coord())
            elif atom.parent.id == residue.id and atom.name == 'N':
                noN = False
                crdN.append(atom.get_coord())
            elif atom.parent.id == residue.id and atom.name == 'C':
                noC = False
                crdC.append(atom.get_coord())
            elif atom.parent.id == residue.id and atom.name == 'O':
                noO = False
                crdO.append(atom.get_coord())

        if noSC:
            print 'missing side chain in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            if resnam[-1] == 'GLY':
                crdSC.append(crdCA[-1])
                sizeSC.append(0)
            else:
                print 'FATAL: missing no side chain detected while the residue is NOT GLYCINE amino acid.'
                sys.exit()
        else:
            # Calculate side chain properties:
            sizeSC.append(len(rescrd_SC))
            crdSC.append(sum(rescrd_SC)/float(sizeSC[-1]))
            #######if sizeSC[-1] != len(resbf_SC):
            #######    print 'something is terribly wrong with the code!: sizeSC[-1] != len(resbf_SC)', sizeSC[-1], len(resbf_SC)
            #######    sys.exit()

        if noCA:
            print 'FATAL: missing CA backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            crdCA.append(['NA','NA','NA'])
            sys.exit()
        if noCB:
            #print 'missing CB backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            if resnam[-1] == 'GLY':
                crdCB.append(crdCA[-1])
            else:
                print 'FATAL: missing CB atom detected while the residue is NOT GLYCINE amino acid.'
                sys.exit()
        if noN:
            print 'missing N backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            crdN.append(crdCA[-1])
        if noC:
            print 'missing C backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            crdC.append(crdCA[-1])
        if noO:
            print 'missing O backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            crdO.append(crdCA[-1])
        
        # Now calculate the Amino Acid properties:
        sizeAA.append(len(rescrd_AA))
        crdAA.append(sum(rescrd_AA)/float(sizeAA[-1]))
        ######if sizeAA[-1] != len(resbf_AA):
        ######    print 'something is terribly wrong with the code!: sizeSC[-1] != len(resbf_SC)', sizeSC[-1], len(resbf_SC)
        ######    sys.exit()

        # Now write out (or append to) the ouput file
        if pdb_chain != reschain[-1]:
            print 'FATAL: residue chain is not A!'

        filenameSC.write( str(resnum[-1]) + '   ' + str(crdSC[-1][0]) + '   ' + str(crdSC[-1][1]) + '   ' + str(crdSC[-1][2]) + '\n' )
        filenameAA.write( str(resnum[-1]) + '   ' + str(crdAA[-1][0]) + '   ' + str(crdAA[-1][1]) + '   ' + str(crdAA[-1][2]) + '\n' )
        filenameCA.write( str(resnum[-1]) + '   ' + str(crdCA[-1][0]) + '   ' + str(crdCA[-1][1]) + '   ' + str(crdCA[-1][2]) + '\n' )
        filenameCB.write( str(resnum[-1]) + '   ' + str(crdCB[-1][0]) + '   ' + str(crdCB[-1][1]) + '   ' + str(crdCB[-1][2]) + '\n' )
        filenameN.write(  str(resnum[-1]) + '   ' + str(crdN[-1][0])  + '   ' + str(crdN[-1][1])  + '   ' + str(crdN[-1][2])  + '\n' )
        filenameC.write(  str(resnum[-1]) + '   ' + str(crdC[-1][0])  + '   ' + str(crdC[-1][1])  + '   ' + str(crdC[-1][2])  + '\n' )
        filenameO.write(  str(resnum[-1]) + '   ' + str(crdO[-1][0])  + '   ' + str(crdO[-1][1])  + '   ' + str(crdO[-1][2])  + '\n' )

    filenameSC.close()
    filenameAA.close()
    filenameCA.close()
    filenameCB.close()
    filenameN.close()
    filenameC.close()
    filenameO.close()
    
    # NOW CALL VORO++ TO CALCULATE VORONOI VOLUMES:
    # OVERALL, 6 VORO++ FILES WILL BE GENERATED. 2 FOR EACH SET OF COORDINATES: CA, SC, AA. EACH OF THE TWO CORRESPONDING TO A DIFFERENT DISTANCE OF THE WALLS FROM THE PDB INSIDE THE BOX.
    
    crdALL = [ numpy.array(crdSC) , numpy.array(crdAA) , numpy.array(crdCA) , numpy.array(crdCB) , numpy.array(crdN) , numpy.array(crdC) , numpy.array(crdO) ]
    filename  = [ filenameSC.name , filenameAA.name , filenameCA.name , filenameCB.name , filenameN.name , filenameC.name , filenameO.name ]
    summary_file = [ sum_SC , sum_AA , sum_CA , sum_CB , sum_N , sum_C , sum_O ]
    wall_distance = [ 30. , 60. ]      # minimum distance of the walls of the box from the pdb structure inside, in units of angstroms
    
    for i,crd in enumerate(crdALL):

        
        for j,distance in enumerate(wall_distance):
            
            # print wall_distance[j]
            # FIRST DETERMINE THE SMALLEST AND THE LARGEST X Y Z VALUES IN THE PDB FILE:
            xmin = min(crd[:,0]) - distance
            ymin = min(crd[:,1]) - distance
            zmin = min(crd[:,2]) - distance
            xmax = max(crd[:,0]) + distance
            ymax = max(crd[:,1]) + distance
            zmax = max(crd[:,2]) + distance

            # NOW USE VORO++ TO FIND THE RESIDUE VOLUMES
            # Reference for the meaning of the flags:  http://math.lbl.gov/voro++/doc/custom.html
            # %w --->  The number of vertices in the Voronoi cell.
            # %g --->  The number of edges of the Voronoi cell.   
            # %E --->  The total edge distance.
            # %s --->  The number of faces of the Voronoi cell.
            # %F --->  The total surface area of the Voronoi cell.
            # %v --->  The volume of the Voronoi cell.
            # %c --->  The centroid of the Voronoi cell, relative to the particle center (three float numbers).
            voro_flags = 'voro++ -o -c "%w %g %E %s %F %v %c" '
            commandline = voro_flags + str(xmin) + ' ' + str(xmax) + ' ' + str(ymin) + ' ' + str(ymax) + ' ' + str(zmin) + ' ' + str(zmax) + ' ' + filename[i]
            process = subprocess.Popen(commandline, shell = True, stdout = subprocess.PIPE)
            process.wait() # Wait until child process is finished.
            #print commandline

            # Now read and parse the output from VORO++
            voro_filename = open( filename[i] + '.vol','r')
            filecontent = voro_filename.readlines(); voro_filename.close()
            
            # data accuracy check point
            if len(filecontent) != len(resnum):
                print pdb_name, '---> while computing Voronoi cells for CA coordinates:' 
                print '   Something fishy was detected: The total number of Voronoi Cells are NOT the same as the total number of residues in the pdb file.'
                sys.exit()
            
            if j == 0:
                pdb_voro_data = []
                for record in filecontent:
                    #print ( record ), len(filecontent), len(resnum), pdb_size
                    residue_voro_data = record.split()[:-3]
                    eccentricity = math.sqrt( float(record.split()[-3])**2 + float(record.split()[-2])**2 + float(record.split()[-1])**2 )
                    residue_voro_data.append(str(eccentricity))
                    pdb_voro_data.append(residue_voro_data)
                    #print residue_voro_data
    
            else:
                 for k,record in enumerate(filecontent):
                     residue_extended_volume = float(record.split()[-4])
                     free_volume = residue_extended_volume - resvol_dict[resnam[k]]
                     if free_volume < 0.0:
                         print 'Negative Voronoi free volume detected!', pdb_name, resnam[k], str(resnum[k]), str(free_volume)
                     pdb_voro_data[k].append(str(free_volume))
                     volume_change_diff  = residue_extended_volume - float(pdb_voro_data[k][-3])
                     volume_change_ratio = residue_extended_volume / float(pdb_voro_data[k][-3])
                     pdb_voro_data[k].append(str(volume_change_diff))
                     pdb_voro_data[k].append(str(volume_change_ratio))

        # Now write out data in the output files. First check if the output file currently exists. If not, then create the output and add the file header.
        if os.path.isfile(summary_file[i]):
            output_file = open( summary_file[i] , 'a' )
        else:
            output_file = open( summary_file[i] , 'w' )
            if filename[i] == filenameSC.name:
                output_file.write( 'pdb' + '\t' + 'resnam' + '\t' + 'resnum' + '\t' + 'sizeSC' + '\t' + 'sizeAA' + '\t' + 'resvol' + '\t' \
                                 + 'VSCnvertices' + '\t' + 'VSCnedges' + '\t' + 'VSCedge_length_total' + '\t' + 'VSCnfaces' + '\t' \
                                 + 'VSCarea' + '\t' + 'VSCvolume' + '\t' + 'VSCeccentricity' + '\t' + 'VSCfree_volume' + '\t' \
                                 + 'VSCvolume_change_diff' + '\t' + 'VSCvolume_change_ratio' + '\n' )
            elif filename[i] == filenameAA.name:
                output_file.write( 'pdb' + '\t' + 'resnam' + '\t' + 'resnum' + '\t' + 'sizeSC' + '\t' + 'sizeAA' + '\t' + 'resvol' + '\t' \
                                 + 'VAAnvertices' + '\t' + 'VAAnedges' + '\t' + 'VAAedge_length_total' + '\t' + 'VAAnfaces' + '\t' \
                                 + 'VAAarea' + '\t' + 'VAAvolume' + '\t' + 'VAAeccentricity' + '\t' + 'VAAfree_volume' + '\t' \
                                 + 'VAAvolume_change_diff' + '\t' + 'VAAvolume_change_ratio' + '\n' )
            elif filename[i] == filenameCA.name:
                output_file.write( 'pdb' + '\t' + 'resnam' + '\t' + 'resnum' + '\t' + 'sizeSC' + '\t' + 'sizeAA' + '\t' + 'resvol' + '\t' \
                                 + 'VCAnvertices' + '\t' + 'VCAnedges' + '\t' + 'VCAedge_length_total' + '\t' + 'VCAnfaces' + '\t' \
                                 + 'VCAarea' + '\t' + 'VCAvolume' + '\t' + 'VCAeccentricity' + '\t' + 'VCAfree_volume' + '\t' \
                                 + 'VCAvolume_change_diff' + '\t' + 'VCAvolume_change_ratio' + '\n' )
            elif filename[i] == filenameCB.name:
                output_file.write( 'pdb' + '\t' + 'resnam' + '\t' + 'resnum' + '\t' + 'sizeSC' + '\t' + 'sizeAA' + '\t' + 'resvol' + '\t' \
                                 + 'VCBnvertices' + '\t' + 'VCBnedges' + '\t' + 'VCBedge_length_total' + '\t' + 'VCBnfaces' + '\t' \
                                 + 'VCBarea' + '\t' + 'VCBvolume' + '\t' + 'VCBccentricity' + '\t' + 'VCBfree_volume' + '\t' \
                                 + 'VCBvolume_change_diff' + '\t' + 'VCBvolume_change_ratio' + '\n' )
            elif filename[i] == filenameN.name:
                output_file.write( 'pdb' + '\t' + 'resnam' + '\t' + 'resnum' + '\t' + 'sizeSC' + '\t' + 'sizeAA' + '\t' + 'resvol' + '\t' \
                                 + 'VNnvertices' + '\t' + 'VNnedges' + '\t' + 'VNedge_length_total' + '\t' + 'VNnfaces' + '\t' \
                                 + 'VNarea' + '\t' + 'VNvolume' + '\t' + 'VNeccentricity' + '\t' + 'VNfree_volume' + '\t' \
                                 + 'VNvolume_change_diff' + '\t' + 'VNvolume_change_ratio' + '\n' )
            elif filename[i] == filenameC.name:
                output_file.write( 'pdb' + '\t' + 'resnam' + '\t' + 'resnum' + '\t' + 'sizeSC' + '\t' + 'sizeAA' + '\t' + 'resvol' + '\t' \
                                 + 'VCnvertices' + '\t' + 'VCnedges' + '\t' + 'VCedge_length_total' + '\t' + 'VCnfaces' + '\t' \
                                 + 'VCarea' + '\t' + 'VCvolume' + '\t' + 'VCeccentricity' + '\t' + 'VCfree_volume' + '\t' \
                                 + 'VCvolume_change_diff' + '\t' + 'VCvolume_change_ratio' + '\n' )
            elif filename[i] == filenameO.name:
                output_file.write( 'pdb' + '\t' + 'resnam' + '\t' + 'resnum' + '\t' + 'sizeSC' + '\t' + 'sizeAA' + '\t' + 'resvol' + '\t' \
                                 + 'VOnvertices' + '\t' + 'VOnedges' + '\t' + 'VOedge_length_total' + '\t' + 'VOnfaces' + '\t' \
                                 + 'VOarea' + '\t' + 'VOvolume' + '\t' + 'VOeccentricity' + '\t' + 'VOfree_volume' + '\t' \
                                 + 'VOvolume_change_diff' + '\t' + 'VOvolume_change_ratio' + '\n' )

        for j,data in enumerate(pdb_voro_data):
            output_file.write( pdb_name + '\t' + resnam[j] + '\t' + str(resnum[j]) + '\t' + str(sizeSC[j]) + '\t' + str(sizeAA[j]) + '\t' + str(resvol_dict[resnam[j]]) + '\t' \
                             + '\t'.join(data) + '\n' )
            #print sizeSC[j], sizeAA[j]
        
        output_file.close()
        
    # NOW REMOVE TEMPORARY FILES
    for item in filename:
        commandline = 'rm ' + item + ' ' + item + '.vol'
        print item, commandline
        process = subprocess.Popen(commandline, shell = True, stdout = subprocess.PIPE)
        process.wait() # Wait until child process is finished.
    
if __name__ == "__main__":
   main()
