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

# Last modified by:  Amir Shahmoradi, Thursday 5:46 PM, July 31 2014, iCMB, UT Austin
# This Python (HPS stands for HydroPhobicity Scale) code takes in the path to a PDB file and outputs in the given path to the ouput file, 3 different HydroPhobicity Scales of all Amino Acids in the PDB file.
# kdh,wwh,hhh stand for the three residue hydrophobicity scales given in https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html
# The following are the three hydrophobicity scales: kdh, wwh, hhh corresponding to the three hydrophobicity scales in https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/hydrophob.html

kdh_dict = {'A':   1.8, 'R':  -4.5, 'N':  -3.5, 'D':  -3.5, \
            'C':   2.5, 'Q':  -3.5, 'E':  -3.5, 'G':  -0.4, \
            'H':  -3.2, 'I':   4.5, 'L':   3.8, 'K':  -3.9, \
            'M':   1.9, 'F':   2.8, 'P':  -1.6, 'S':  -0.8, \
            'T':  -0.7, 'W':  -0.9, 'Y':  -1.3, 'V':   4.2}

wwh_dict = {'A': -0.17, 'R': -0.81, 'N': -0.42, 'D': -1.23, \
            'C':  0.24, 'Q': -0.58, 'E': -2.02, 'G': -0.01, \
            'H': -0.96, 'I':  0.31, 'L':  0.56, 'K': -0.99, \
            'M':  0.23, 'F':  1.13, 'P': -0.45, 'S': -0.13, \
            'T': -0.14, 'W':  1.85, 'Y':  0.94, 'V': -0.07}

hhh_dict = {'A':  0.11, 'R':  2.58, 'N':  2.05, 'D':  3.49, \
            'C': -0.13, 'Q':  2.36, 'E':  2.68, 'G':  0.74, \
            'H':  2.06, 'I': -0.60, 'L': -0.55, 'K':  2.71, \
            'M': -0.10, 'F': -0.32, 'P':  2.23, 'S':  0.84, \
            'T':  0.52, 'W':  0.30, 'Y':  0.68, 'V': -0.31}

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
res_dict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
             'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
             'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
             'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

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
       sum_out_file.write( 'pdb' + '\t' + 'resnam' + '\t' + 'resnum' + '\t' + 'hpskd' + '\t' + 'hpsww' + '\t' + 'hpshh' + '\n' )

    p = PDBParser()
    pdb_name  = pdb_in[-10:-4]
    pdb_chain = pdb_in[-5:-4]
    structure = p.get_structure(pdb_name,pdb_in)
    

    for residue in structure.get_residues():
        resname = residue.resname
        resnum  = str(residue.id[1])
        AA1     = res_dict[resname]
        hpskd   = str(kdh_dict[AA1])
        hpsww   = str(wwh_dict[AA1])
        hpshh   = str(hhh_dict[AA1])
        sum_out_file.write(pdb_name + '\t' + resname + '\t' + resnum + '\t' + hpskd + '\t' + hpsww + '\t' + hpshh + '\n')
        if pdb_chain != residue.parent.id:
            print 'FATAL: residue chain is not A!'
            sys.exit()
        

if __name__ == "__main__":
   main()