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

# written by Amir Shahmoradi, Monday 7:41 PM, October 20 2014, Wilke Lab, iCMB, The University of Texas at Austin.
# I am writing this code to find the number of sequences in the alignments for each protein in the dataset.
# The goal is to see wether the number of sequence in the alignments has any significant effect on the observed sequence-structure correlations among the proteins.

pdb = OrderedDict([('name'        , 'NA'), ('naln'     , 'NA')])

def main():
    if len( sys.argv ) != 3:
        print '''
 
Usage:'''
        print "     ", sys.argv[0], "<input pdb file>", "<output dssp file>", "<output summary file>", '\n'
        sys.exit('Program stopped.\n')
    else:
        tree_path = sys.argv[1]  # path to the input pdb file
        sum_out = sys.argv[2]   # summary file containing the sequences in the alignments each pdb.

    input = open(tree_path, 'r')
    fileContent = input.readline()
    pdb['name'] = tree_path.split('/')[-1][0:6] # split element 4 of the filecontent list (space-delimited). the 0th element is number of residues in pdb
    pdb['naln'] = fileContent.split()[0] # read the number of sequence alingments

    # Now check if the output summary file alreadu exists. If so, then append data for the new pdb to the data already in the file.
    if os.path.isfile(sum_out):
       sum_out_file = open(sum_out,'a')
    else:
       sum_out_file = open(sum_out,'w')
       sum_out_file.write( 'pdb' + '\t' + 'naln' + '\n' )

    for key in pdb:
        sum_out_file.write(str(pdb[key]) + '\t')
    sum_out_file.write('\n')
    
    return

if __name__ == "__main__":
   main()
