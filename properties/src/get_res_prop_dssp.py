#!/usr/bin/python
import sys, subprocess, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.Polypeptide import *
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import Dice
from collections import OrderedDict
import math
import numpy

# I decided to break data collection into small pieces and hope for a correct sample at the end, instead of writting a comprehensive code that outputs all information in one call.
# So this code extracts residue properties only from dssp files.
# My assumption is that the ELJ entropy file contains the same residues as found in DSSP files and there is no inconsistency between the orders of the two.
# If there is, then I will have to change the code acordingly.
#
# INPUT:  dssp files in ../dssp_in/*  and the name of the output file.
# OUTPUT: A file containing all relevant DSSP residue properties including RSA (in addition to ASA), in chronological order of pdbs in the directory.
# Amir Shahmoradi, 11:08 AM, Monday Aug 4 2014, Wilke Lab, iCMB, The University of Texas at Austin.

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

#ASA normalization constants were taken from: M. Z. Tien, A. G. Meyer, D. K. Sydykova, S. J. Spielman, C. O. Wilke (2013). Maximum allowed solvent accessibilities of residues in proteins. PLOS ONE 8:e80635.
#This dictionary is not needed at the moment, but I am keeping it here for future use.
residue_max_acc = {'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, \
                   'C': 167.0, 'Q': 225.0, 'E': 223.0, 'G': 104.0, \
                   'H': 224.0, 'I': 197.0, 'L': 201.0, 'K': 236.0, \
                   'M': 224.0, 'F': 240.0, 'P': 159.0, 'S': 155.0, \
                   'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0}
    
def main():
    if len( sys.argv ) != 3:
        print '''
 
Usage:'''
        print "     ", sys.argv[0], "<input dssp file>", "<output summary file>", '\n'
        sys.exit('Program stopped.\n')
    else:
        dssp_in = sys.argv[1]  # path for the output dssp file to be read
        sum_out = sys.argv[2]   # summary file containing all residue DSSP-property data

    #pdb = get_res_props(pdb_path,dssp_in)

    # Now check if the output summary file alreadu exists. If so, then append data for the new pdb to the data already in the file.
    if os.path.isfile(sum_out):
       sum_out_file = open(sum_out,'a')
    else:
       sum_out_file = open(sum_out,'w')
       sum_out_file.write('pdb' + '\t' + 'resnam' + '\t' + 'resnum' + '\t' + 'asa' + '\t' + 'rsa' + '\t' + 'hbe_mean' + '\t' + 'rss' + '\n')
    
    ######################################################################################################################################
    # Now calculate the mean of all 4 hydrogen bonds energies and residue RSAs for each individual residue:
    
    input = open(dssp_in, 'r')
    fileContents = input.readlines()   # This is a list containing each line of the input file as an element.
    pdb_name = dssp_in[-11:-5]
    
    resnam = []     # A list containing all Residue Names
    resnum = []     # A list containing all Residue Numbers
    asa = []        # A list containing all normalized residue ASA values in the pdb file
    rsa = []        # A list containing all normalized residue ASA values in the pdb file
    crd = []        # A list containing all CA-atom WCN in the pdb file
    hbe = []        # A list containing average hydrogen bonds energies for each residue in the pdb file
    rss = []        # A list containing Residue Secondary Structure
    counter = 0
    for record in fileContents[25:len(fileContents)]:
        AA = record[13]
        if AA != ('!' or '*'):
            counter += 1
            if AA == ('X'):
                print 'potential Histidine variant HID, HIE or HIP found in pdb file ', dssp_in, 'residue number ', counter
                AA = 'H'    # correct Histidine variants code to the standard code
                #raw_input('Press <Enter> to continue ...')
            resnam.append(AA)
            resnum.append(record[6:10])
            asa.append(record[35:39])
            rsa.append(float(record[35:39])/residue_max_acc[AA])
            crd.append( numpy.array( [ float(record.split()[-3]) , float(record.split()[-2]) , float(record.split()[-1]) ] ) ) # stores CA atom triplet coordinates as individual elements of the list crd.
            hbe1 = float(record[38:51].split(',')[1])
            hbe2 = float(record[50:62].split(',')[1])
            hbe3 = float(record[61:73].split(',')[1])
            hbe4 = float(record[72:84].split(',')[1])
            hbe.append(numpy.mean([hbe1,hbe2,hbe3,hbe4]))  # This is the average of the four possible hydrogen bonds of a single residue. If a residue has no hydrogen bond, then all four are zero, so the mean is also zero.
            if record[16:17] == ' ' :
                rss.append('L')
            else :
                rss.append(record[16:17])

    # Now calculate residue wcn :  *****ATTN***** AUG 25 2014: wcnca is not needed to be claculated here anmymore. Python is too slow for this kind of calculation. All WCNs for different atoms are now calculated by a Fortran code separately.
    #wcn = []   # A list containing all CA-atom WCN in the pdb file
    #for i in range(len(crd)) :
    #    sum_terms = 0.
    #    for j in range(len(crd)) :
    #        if i != j :
    #            sum_terms += 1./( (crd[i][0]-crd[j][0])**2 + (crd[i][1]-crd[j][1])**2 + (crd[i][2]-crd[j][2])**2 )
    #    wcn.append(sum_terms)
    input.close()
    
    # Now write out (or append to) the ouput file
    for (i,j) in enumerate(rsa):
        sum_out_file.write(pdb_name + '\t' + resnam[i] + '\t' + resnum[i] + '\t' + asa[i] + '\t' + str(rsa[i]) + '\t' + str(hbe[i]) + '\t' + rss[i] + '\n')
    
    return

if __name__ == "__main__":
   main()