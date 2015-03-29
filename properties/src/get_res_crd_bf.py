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

# This code reads in a given input PDB file, and the name of the atom representing the individual residues in proteins. Then outputs, in a given output file, the residue names and numbers and the coordinates of the representative atoms in all residues and their corresponding B factors. There will be 7 output files, corresponding to atom names: C, CA, O, N, CB, SC: the Center-Of-Mass (COM) of the side chain, and AA: the COM of the BackBone and side chain atoms of each residue.
# The B factor for SC and AA are calculated by averaging over all the corresponding atomic B factors.
# INPUT:  pdb files in ../structures/*  and the name of the 7 output files in a VERY SPECIFIC ORDER!.
# OUTPUT: 7 files containing all relevant residue properties described above, in chronological order of pdbs in the directory.
# ATTN: If the output files do not exist, the requested output files will be created. Otherwise, the output will be appended to the existing output files and no HEADER for the definitions of columns will be added to the outfile!
# ATTN: If a residue does not have the requested atom, then the properties of the corresponding CA atom for that residue will be instead reported.
# Amir Shahmoradi, 11:10 PM, Saturday March 28 2015, Wilke Lab, iCMB, The University of Texas at Austin.
    
def main():
    if len( sys.argv ) != 9:
        print '''
 
Usage:'''
        print "     ", sys.argv[0], "<input PDB file>", "<output file for SC coordinates>" , "<output file for AA coordinates>" , "<output file for AA coordinates>" , "<output file for CA coordinates>" , "<output file for CB coordinates>" , "<output file for N coordinates>" , "<output file for C coordinates>" , "<output file for O coordinates>" , '\n'
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


    # Now check if the output summary file alreadu exists. If so, then append data for the new pdb to the data already in the file.
    if os.path.isfile(sum_SC):
       sum_SC_file = open(sum_SC,'a')
    else:
       sum_SC_file = open(sum_SC,'w')
       sum_SC_file.write('pdb'   + '\t' + 'resnam'+ '\t' + 'resnum' + '\t' + 'x' + '\t' + 'y' + '\t' + 'z' + '\t' + 'bf' + '\n')
       
    if os.path.isfile(sum_AA):
       sum_AA_file = open(sum_AA,'a')
    else:
       sum_AA_file = open(sum_AA,'w')
       sum_AA_file.write('pdb'   + '\t' + 'resnam'+ '\t' + 'resnum' + '\t' + 'x' + '\t' + 'y' + '\t' + 'z' + '\t' + 'bf' + '\n')

    if os.path.isfile(sum_CA):
       sum_CA_file = open(sum_CA,'a')
    else:       
       sum_CA_file = open(sum_CA,'w')
       sum_CA_file.write('pdb'   + '\t' + 'resnam'+ '\t' + 'resnum' + '\t' + 'x' + '\t' + 'y' + '\t' + 'z' + '\t' + 'bf' + '\n')

    if os.path.isfile(sum_CB):
       sum_CB_file = open(sum_CB,'a')
    else:
       sum_CB_file = open(sum_CB,'w')
       sum_CB_file.write('pdb'   + '\t' + 'resnam'+ '\t' + 'resnum' + '\t' + 'x' + '\t' + 'y' + '\t' + 'z' + '\t' + 'bf' + '\n')
       
    if os.path.isfile(sum_N):
       sum_N_file = open(sum_N,'a')
    else:
       sum_N_file = open(sum_N,'w')
       sum_N_file.write('pdb'   + '\t' + 'resnam'+ '\t' + 'resnum' + '\t' + 'x' + '\t' + 'y' + '\t' + 'z' + '\t' + 'bf' + '\n')
       
    if os.path.isfile(sum_C):
       sum_C_file = open(sum_C,'a')
    else:
       sum_C_file = open(sum_C,'w')
       sum_C_file.write('pdb'   + '\t' + 'resnam'+ '\t' + 'resnum' + '\t' + 'x' + '\t' + 'y' + '\t' + 'z' + '\t' + 'bf' + '\n')
       
    if os.path.isfile(sum_O):
       sum_O_file = open(sum_O,'a')
    else:
       sum_O_file = open(sum_O,'w')
       sum_O_file.write('pdb'   + '\t' + 'resnam'+ '\t' + 'resnum' + '\t' + 'x' + '\t' + 'y' + '\t' + 'z' + '\t' + 'bf' + '\n')

    p = PDBParser()
    pdb_name = pdb_in[-10:-4]
    pdb_chain = pdb_in[-5:-4]
    structure = p.get_structure(pdb_name,pdb_in)
    
    resnam   = []     # A list containing all Residue Names
    resnum   = []     # A list containing all Residue Numbers
    reschain = []     # A list containing the chain name in which the residues lie
    crdN    = []
    crdCA   = []
    crdC    = []
    crdO    = []
    crdCB   = []
    crdAA   = []
    crdSC   = []
    bfN     = []
    bfCA    = []
    bfC     = []
    bfO     = []
    bfCB    = []
    bfAA    = []
    bfSC    = []
    sizeSC  = []   # A list containing the total number of atoms in each residue Side Chain (SC).
    sizeAA  = []   # A list containing the total number of atoms in each Amino Acid (AA).
    
    #Ncounter  = 0
    #CAcounter = 0
    #Ccounter  = 0
    #Ocounter  = 0
    #CBcounter = 0

    for residue in structure.get_residues():
        #print residue
        resnam.append(residue.resname)
        resnum.append(residue.get_full_id()[3][1])
        reschain.append(residue.get_full_id()[2])
        noN  = True
        noCA = True
        noC  = True
        noO  = True
        noCB = True
        noSC = True
        rescrd_SC = []  # A list containing the coordinates of all side chain atoms of the current residue. Will be used to calculate the COM of the side chain.
        rescrd_AA = []  # A list containing the coordinates of all atoms of the current Amino Acid. Will be used to calculate the COM of the Amino Acid.
        resbf_SC  = []  # A list containing the Bfactors of all side chain atoms of the current residue. Will be used to calculate the side chain average Bfactor.
        resbf_AA  = []  # A list containing the Bfactors of all Amino Acid atoms of the current Amino Acid. Will be used to calculate the Amino Acid average Bfactor.
        for atom in structure.get_atoms():
            # atom.name is equivalent to atom.get_id()
            if atom.parent.id == residue.id and atom.name == 'N':
                noN = False
                crdN.append(atom.get_coord())
                bfN.append(atom.get_bfactor())
            elif atom.parent.id == residue.id and atom.name == 'CA':
                noCA = False
                crdCA.append(atom.get_coord())
                bfCA.append(atom.get_bfactor())
            elif atom.parent.id == residue.id and atom.name == 'C':
                noC = False
                #Ccounter += 1
                #print Ccounter
                crdC.append(atom.get_coord())
                bfC.append(atom.get_bfactor())
            elif atom.parent.id == residue.id and atom.name == 'O':
                noO = False
                crdO.append(atom.get_coord())
                bfO.append(atom.get_bfactor())
            elif atom.parent.id == residue.id and atom.name == 'CB':
                noCB = False
                crdCB.append(atom.get_coord())
                bfCB.append(atom.get_bfactor())
            
            if atom.parent.id == residue.id and atom.name not in ['C','CA','O','N']:
                noSC = False
                rescrd_SC.append(atom.get_coord())
                resbf_SC.append(atom.get_bfactor())
            
            if atom.parent.id == residue.id:
                rescrd_AA.append(atom.get_coord())
                resbf_AA.append(atom.get_bfactor())

        if noN:
            print 'missing N backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            crdN.append(crdCA[-1])
            bfN.append(bfCA[-1])
        if noCA:
            print 'FATAL: missing CA backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            crdCA.append(['NA','NA','NA'])
            bfCA.append('NA')
            sys.exit()
        if noC:
            print 'missing C backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            crdC.append(crdCA[-1])
            bfC.append(bfCA[-1])
        if noO:
            print 'missing O backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            crdO.append(crdCA[-1])
            bfO.append(bfCA[-1])
        if noCB:
            #print 'missing CB backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            if resnam[-1] == 'GLY':
                crdCB.append(crdCA[-1])
                bfCB.append(bfCA[-1])
            else:
                print 'FATAL: missing CB atom detected while the residue is NOT GLYCINE amino acid.'
                sys.exit()
        if noSC:
            print 'missing side chain in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            if resnam[-1] == 'GLY':
                crdSC.append(crdCA[-1])
                bfSC.append(bfCA[-1])
                sizeSC.append(0)
            else:
                print 'FATAL: missing no side chain detected while the residue is NOT GLYCINE amino acid.'
                sys.exit()
        else:
            # Calculate side chain properties:
            sizeSC.append(len(rescrd_SC))
            crdSC.append(sum(rescrd_SC)/float(sizeSC[-1]))
            bfSC.append(sum(resbf_SC)/float(sizeSC[-1]))
            if sizeSC[-1] != len(resbf_SC):
                print 'something is terribly wrong with the code!: sizeSC[-1] != len(resbf_SC)', sizeSC[-1], len(resbf_SC)
                sys.exit()
        
        # Now calculate the Amino Acid properties:
        sizeAA.append(len(rescrd_AA))
        crdAA.append(sum(rescrd_AA)/float(sizeAA[-1]))
        bfAA.append(sum(resbf_AA)/float(sizeAA[-1]))
        if sizeAA[-1] != len(resbf_AA):
            print 'something is terribly wrong with the code!: sizeSC[-1] != len(resbf_SC)', sizeSC[-1], len(resbf_SC)
            sys.exit()

    # Now write out (or append to) the ouput file
    for i in range(len(resnam)):
        if pdb_chain != reschain[i]:
            print 'FATAL: residue chain is not A!'
        sum_N_file.write(pdb_name + '\t' + resnam[i] + '\t' + str(resnum[i]) + '\t' + str(crdN[i][0]) + '\t' + str(crdN[i][1]) + '\t' + str(crdN[i][2])  + '\t' + str(bfN[i]) + '\n')
        sum_C_file.write(pdb_name + '\t' + resnam[i] + '\t' + str(resnum[i]) + '\t' + str(crdC[i][0]) + '\t' + str(crdC[i][1]) + '\t' + str(crdC[i][2])  + '\t' + str(bfC[i]) + '\n')
        sum_CA_file.write(pdb_name + '\t' + resnam[i] + '\t' + str(resnum[i]) + '\t' + str(crdCA[i][0]) + '\t' + str(crdCA[i][1]) + '\t' + str(crdCA[i][2])  + '\t' + str(bfCA[i]) + '\n')
        sum_O_file.write(pdb_name + '\t' + resnam[i] + '\t' + str(resnum[i]) + '\t' + str(crdO[i][0]) + '\t' + str(crdO[i][1]) + '\t' + str(crdO[i][2])  + '\t' + str(bfO[i]) + '\n')
        sum_CB_file.write(pdb_name + '\t' + resnam[i] + '\t' + str(resnum[i]) + '\t' + str(crdCB[i][0]) + '\t' + str(crdCB[i][1]) + '\t' + str(crdCB[i][2])  + '\t' + str(bfCB[i]) + '\n')
        sum_SC_file.write(pdb_name + '\t' + resnam[i] + '\t' + str(resnum[i]) + '\t' + str(crdSC[i][0]) + '\t' + str(crdSC[i][1]) + '\t' + str(crdSC[i][2])  + '\t' + str(bfSC[i]) + '\n')
        sum_AA_file.write(pdb_name + '\t' + resnam[i] + '\t' + str(resnum[i]) + '\t' + str(crdAA[i][0]) + '\t' + str(crdAA[i][1]) + '\t' + str(crdAA[i][2])  + '\t' + str(bfAA[i]) + '\n')
    
if __name__ == "__main__":
   main()


    
#crd.append( numpy.array( [ float(record.split()[-3]) , float(record.split()[-2]) , float(record.split()[-1]) ] ) ) # stores CA atom triplet coordinates as individual elements of the list crd.
## Now calculate residue wcn :
#wcn = []   # A list containing all CA-atom WCN in the pdb file
#for i in range(len(crd)) :
#    sum_terms = 0.
#    for j in range(len(crd)) :
#        if i != j :
#            sum_terms += 1./( (crd[i][0]-crd[j][0])**2 + (crd[i][1]-crd[j][1])**2 + (crd[i][2]-crd[j][2])**2 )
#    wcn.append(sum_terms)
#input.close()
