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

# This code reads in a given input PDB file and outputs, in a given output file, the coordinates of backbone atoms: C, CA, O, N and the CB atom and the Center-Of-Mass (COM) of the side chains and the COM of the BackBone of each residue.
# Also on the output is the Bfactors for each of the corresponding atoms and the average Bfactor for the case Side-Chain (SC) COM and the entire Amino Acid (AA) COM (including backbone atoms).
#
# INPUT:  pdb files in ../structures/*  and the name of the output file.
# OUTPUT: A file containing all relevant residue properties described above, in chronological order of pdbs in the directory.
# Amir Shahmoradi, 12:10 PM, Tuesday Aug 12 2014, Wilke Lab, iCMB, The University of Texas at Austin.
    
def main():
    if len( sys.argv ) != 3:
        print '''
 
Usage:'''
        print "     ", sys.argv[0], "<input PDB file>", "<output summary file>", '\n'
        sys.exit('Program stopped.\n')
    else:
        pdb_in  = sys.argv[1]  # path for the input PDB file to be read
        sum_out = sys.argv[2]   # summary file containing all residue coordinate & Bfactor data


    # Now check if the output summary file alreadu exists. If so, then append data for the new pdb to the data already in the file.
    if os.path.isfile(sum_out):
       sum_out_file = open(sum_out,'a')
    else:
       sum_out_file = open(sum_out,'w')
       sum_out_file.write('pdb' + '\t' + 'chain' + '\t' + 'resnam' + '\t' + 'resnum' + '\t' + 'sizeSC' + '\t' + 'sizeAA' + '\t'  \
                          'wcnSC' + '\t' + 'SC_bf'  + '\t' + \
                          'wcnAA' + '\t' + 'AA_bf'  + '\t' + \
                          'wcnN'  + '\t' + 'N_bf'   + '\t' + \
                          'wcnCA' + '\t' + 'CA_bf'  + '\t' + \
                          'wcnC'  + '\t' + 'C_bf'   + '\t' + \
                          'wcnO'  + '\t' + 'O_bf'   + '\t' + \
                          'wcnCB' + '\t' + 'CB_bf'  + '\n' )
                          #'bfc' + '\t' + 'bfca' + '\t' + 'bfo' + '\t' + 'bfn' + '\t' + 'bfcb' + '\t' + 'bfaa' + '\t' + 'bfsc' + '\n')
                          #'N.x' + '\t' + 'N.y' + '\t' + 'N.z' + '\t' + \
                          #'CA.x' + '\t' + 'CA.y' + '\t' + 'CA.z' + '\t' + \
                          #'C.x' + '\t' + 'C.y' + '\t' + 'C.z' + '\t' + \
                          #'O.x' + '\t' + 'O.y' + '\t' + 'O.z' + '\t' + \
                          #'CB.x' + '\t' + 'CB.y' + '\t' + 'CB.z' + '\t' + \
                          #'SC.x' + '\t' + 'SC.y' + '\t' + 'SC.z' + '\t' + \
                          #'AA.x' + '\t' + 'AA.y' + '\t' + 'AA.z' + '\t' + \
                          #'bfc' + '\t' + 'bfca' + '\t' + 'bfo' + '\t' + 'bfn' + '\t' + 'bfcb' + '\t' + 'bfaa' + '\t' + 'bfsc' + '\n')
                          
       #sum_out_file.write('pdb' + '\t' + 'chain' + '\t' + 'resnam' + '\t' + 'resnum' + '\t' + 'wcnc' + '\t' + 'wcnca' + '\t' + 'wcno' + '\t' + 'wcnn' + '\t' + 'wcncb' + '\t' + 'wcnaa' + '\t' + 'wcnsc' + '\t' + 'bfc' + '\t' + 'bfca' + '\t' + 'bfo' + '\t' + 'bfn' + '\t' + 'bfcb' + '\t' + 'bfaa' + '\t' + 'bfsc' + '\t' + 'sc_size' + '\n')

    p = PDBParser()
    pdb_name = pdb_in[-10:-6]
    pdb_chain = pdb_in[-5:-4]
    structure = p.get_structure(pdb_name,pdb_in)
    
    resnam   = []     # A list containing all Residue Names
    resnum   = []     # A list containing all Residue Numbers
    reschain = []     # A list containing the chain name in which the residues lie
    N_crd    = []
    CA_crd   = []
    C_crd    = []
    O_crd    = []
    CB_crd   = []
    AA_crd   = []
    SC_crd   = []
    N_bf     = []
    CA_bf    = []
    C_bf     = []
    O_bf     = []
    CB_bf    = []
    AA_bf    = []
    SC_bf    = []
    SC_size  = []   # A list containing the total number of atoms in each residue Side Chain (SC).
    AA_size  = []   # A list containing the total number of atoms in each Amino Acid (AA).
    
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
            if atom.get_full_id()[3][1] == resnum[-1] and atom.name == 'N':
                noN = False
                N_crd.append(atom.get_coord())
                N_bf.append(atom.get_bfactor())
            elif atom.get_full_id()[3][1] == resnum[-1] and atom.name == 'CA':
                noCA = False
                CA_crd.append(atom.get_coord())
                CA_bf.append(atom.get_bfactor())
            elif atom.get_full_id()[3][1] == resnum[-1] and atom.name == 'C':
                noC = False
                #Ccounter += 1
                #print Ccounter
                C_crd.append(atom.get_coord())
                C_bf.append(atom.get_bfactor())
            elif atom.get_full_id()[3][1] == resnum[-1] and atom.name == 'O':
                noO = False
                O_crd.append(atom.get_coord())
                O_bf.append(atom.get_bfactor())
            elif atom.get_full_id()[3][1] == resnum[-1] and atom.name == 'CB':
                noCB = False
                CB_crd.append(atom.get_coord())
                CB_bf.append(atom.get_bfactor())
            
            if atom.get_full_id()[3][1] == resnum[-1] and atom.name not in ['C','CA','O','N']:
                noSC = False
                rescrd_SC.append(atom.get_coord())
                resbf_SC.append(atom.get_bfactor())
            
            if atom.get_full_id()[3][1] == resnum[-1]:
                rescrd_AA.append(atom.get_coord())
                resbf_AA.append(atom.get_bfactor())

        if noN:
            print 'missing N backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            N_crd.append(CA_crd[-1])
            N_bf.append(CA_bf[-1])
        if noCA:
            print 'FATAL: missing CA backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            CA_crd.append(['NA','NA','NA'])
            CA_bf.append('NA')
            sys.exit()
        if noC:
            print 'missing C backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            C_crd.append(CA_crd[-1])
            C_bf.append(CA_bf[-1])
        if noO:
            print 'missing O backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            O_crd.append(CA_crd[-1])
            O_bf.append(CA_bf[-1])
        if noCB:
            #print 'missing CB backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            if resnam[-1] == 'GLY':
                CB_crd.append(CA_crd[-1])
                CB_bf.append(CA_bf[-1])
            else:
                print 'FATAL: missing CB atom detected while the residue is NOT GLYCINE amino acid.'
                sys.exit()
        if noSC:
            print 'missing side chain in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_in[-10:-4]
            if resnam[-1] == 'GLY':
                SC_crd.append(CA_crd[-1])
                SC_bf.append(CA_bf[-1])
                SC_size.append(0)
            else:
                print 'FATAL: missing no side chain detected while the residue is NOT GLYCINE amino acid.'
                sys.exit()
        else:
            # Calculate side chain properties:
            SC_size.append(len(rescrd_SC))
            SC_crd.append(sum(rescrd_SC)/float(SC_size[-1]))
            SC_bf.append(sum(resbf_SC)/float(SC_size[-1]))
            if SC_size[-1] != len(resbf_SC):
                print 'something is terribly wrong with the code!: SC_size[-1] != len(resbf_SC)', SC_size[-1], len(resbf_SC)
                sys.exit()
        
        # Now calculate the Amino Acid properties:
        AA_size.append(len(rescrd_AA))
        AA_crd.append(sum(rescrd_AA)/float(AA_size[-1]))
        AA_bf.append(sum(resbf_AA)/float(AA_size[-1]))
        if AA_size[-1] != len(resbf_AA):
            print 'something is terribly wrong with the code!: SC_size[-1] != len(resbf_SC)', SC_size[-1], len(resbf_SC)
            sys.exit()

    # Now calcualte the Contact numbers for differnt sets of coordinates and output the results :
    
    wcnN     = []
    wcnCA    = []
    wcnC     = []
    wcnO     = []
    wcnCB    = []
    wcnSC    = []
    wcnAA    = []
    
    for i in range(len(resnam)):
        
        wcnNi = 0.      # WCN for atom N of the ith residue in the PDB file.
        for j in range(len(N_crd)) :
            if i != j :
                wcnNi += 1./( (N_crd[i][0]-N_crd[j][0])**2 + (N_crd[i][1]-N_crd[j][1])**2 + (N_crd[i][2]-N_crd[j][2])**2 )
        wcnN.append(wcnNi)
        
        wcnCAi = 0.      # WCN for atom CA of the ith residue in the PDB file.
        for j in range(len(CA_crd)) :
            if i != j :
                wcnCAi += 1./( (CA_crd[i][0]-CA_crd[j][0])**2 + (CA_crd[i][1]-CA_crd[j][1])**2 + (CA_crd[i][2]-CA_crd[j][2])**2 )
        wcnCA.append(wcnCAi)
        
        wcnCi = 0.      # WCN for atom C of the ith residue in the PDB file.
        for j in range(len(C_crd)) :
            if i != j :
                wcnCi += 1./( (C_crd[i][0]-C_crd[j][0])**2 + (C_crd[i][1]-C_crd[j][1])**2 + (C_crd[i][2]-C_crd[j][2])**2 )
        wcnC.append(wcnCi)
        
        wcnOi = 0.      # WCN for atom O of the ith residue in the PDB file.
        for j in range(len(O_crd)) :
            if i != j :
                wcnOi += 1./( (O_crd[i][0]-O_crd[j][0])**2 + (O_crd[i][1]-O_crd[j][1])**2 + (O_crd[i][2]-O_crd[j][2])**2 )
        wcnO.append(wcnOi)
        
        wcnCBi = 0.      # WCN for atom CB of the ith residue in the PDB file.
        for j in range(len(CB_crd)) :
            if i != j :
                wcnCBi += 1./( (CB_crd[i][0]-CB_crd[j][0])**2 + (CB_crd[i][1]-CB_crd[j][1])**2 + (CB_crd[i][2]-CB_crd[j][2])**2 )
        wcnCB.append(wcnCBi)
        
        wcnSCi = 0.      # WCN for atom SC of the ith residue in the PDB file.
        for j in range(len(SC_crd)) :
            if i != j :
                wcnSCi += 1./( (SC_crd[i][0]-SC_crd[j][0])**2 + (SC_crd[i][1]-SC_crd[j][1])**2 + (SC_crd[i][2]-SC_crd[j][2])**2 )
        wcnSC.append(wcnSCi)
        
        wcnAAi = 0.      # WCN for atom N of the ith residue in the PDB file.
        for j in range(len(AA_crd)) :
            if i != j :
                wcnAAi += 1./( (AA_crd[i][0]-AA_crd[j][0])**2 + (AA_crd[i][1]-AA_crd[j][1])**2 + (AA_crd[i][2]-AA_crd[j][2])**2 )
        wcnAA.append(wcnAAi)
        
        # Now write out (or append to) the ouput file
        if pdb_chain != reschain[i]:
            print 'FATAL: residue chain is not A!'
        sum_out_file.write(pdb_name + '\t' + pdb_chain + '\t' + resnam[i] + '\t' + str(resnum[i]) + '\t' + str(SC_size[i]) + '\t' + str(AA_size[i]) + '\t' + \
                           str(wcnSC[i])  + '\t' + str(SC_bf[i])  + '\t' + \
                           str(wcnAA[i])  + '\t' + str(AA_bf[i])  + '\t' + \
                           str( wcnN[i])  + '\t' + str( N_bf[i])  + '\t' + \
                           str(wcnCA[i])  + '\t' + str(CA_bf[i])  + '\t' + \
                           str( wcnC[i])  + '\t' + str( C_bf[i])  + '\t' + \
                           str( wcnO[i])  + '\t' + str( O_bf[i])  + '\t' + \
                           str(wcnCB[i])  + '\t' + str(CB_bf[i])  + '\n' )
                           #str(N_crd[i])  + '\t' + str(N_bf[i])  + '\t' + \
                           #str(CA_crd[i]) + '\t' + str(CA_bf[i]) + '\t' + \
                           #str(C_crd[i])  + '\t' + str(C_bf[i])  + '\t' + \
                           #str(O_crd[i])  + '\t' + str(O_bf[i])  + '\t' + \
                           #str(CB_crd[i]) + '\t' + str(CB_bf[i]) + '\t' + \
                           #str(SC_crd[i]) + '\t' + str(SC_bf[i]) + '\t' + \
                           #str(AA_crd[i]) + '\t' + str(AA_bf[i]) + '\t' + '\n')
    
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