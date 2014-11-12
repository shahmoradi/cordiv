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

# Amir Shahmoradi, 6:59 PM, Tuesday July 1, 2014, Wilke Lab, iCMB, The University of Texas at Austin.
# ATTN, AUG 25, 2014: There is a newer version of this code that writes out data in LONG TABLE format. Use the new code instead.
# I am writing this code to read run DSSP on a given input pdb file. The code then reads in DSSP output and extracts the useful information about the overall protein structure and stores it in a given output file by the user.

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

pdb = OrderedDict([('name'        , 'NA'), \
                   ('asa_pdb'     , 'NA'), \
                   ('nres'        , 'NA'), \
                   ('nchain'     , 'NA'), \
                   ('sum_nssb'    , 'NA'), \
                   ('mean_nssb'   , 'NA'), \
                   ('sum_nhbon'   , 'NA'), \
                   ('mean_nhbon'  , 'NA'), \
                   ('sum_nhbps'   , 'NA'), \
                   ('mean_nhbps'  , 'NA'), \
                   ('sum_nhbas'   , 'NA'), \
                   ('mean_nhbas'  , 'NA')  \
###                ('mhbe'        , 'NA'), \
###                ('vhbe'        , 'NA'), \
                  #('nah'         , 'NA'), \
                  #('nresah'      , 'NA'), \
                  #('nbs'         , 'NA'), \
###                ('mrsa'        , 'NA'), \
###                ('vrsa'        , 'NA'), \
###                ('mwcn'        , 'NA'), \
###                ('vwcn'        , 'NA')
                  ])
# Definition of variables:
  # name        : pdb name
  # nres        : Number of RESidues in pdb
  # nchain     : Number of CHAINS in pdb. Default is assumed to be 1 for the expected data set of mine.
  # nssb        : Number of SS Bridges in the pdb structure
  # asa         : total Accessible Surface Area of the protein
  # nhbon       : total Number of Hydrogen Bonds of the form O --> H-N(j), j = -5,5
  # nhbps       : total Number of Hydrogen Bonds in Parallel beta Sheets
  # nhbas       : total Number of Hydrogen Bonds in Antiparallel beta Sheets
  # mhbe        : Mean of the Hydrogen Bonds Energies in the protein
  # vhbe        : Variance of the Hydrogen Bonds Energies in the protein
  # nah         : total Number of Alpha-Helices in the protein
  # nresah      : total Number of RESidues in all Alpha-Helices in the protein
  # nbs         : total Number of Beta Sheets in the protein
  # mrsa        : Mean value of the Relative Solvent Accessibility of all residues in the protein
  # vrsa        : Variance of the Relative Solvent Accessibility of all residues in the protein
  # mwcn        : Mean value of the Weighted Contact Number of all CA atoms in the protein
  # vwcn        : Variance of the Weighted Contact Number of all CA atoms in the protein

#This calculates the RSA values for a PDB using DSSP. It returns a list of the Amino Acids and a list of their RSA values. 
def get_pdb_props(pdb_path,dssp_out):
    commandline = '$HOME/dssp' + ' -i ' + pdb_path + ' -o ' + dssp_out
    process = subprocess.Popen(commandline, shell = True, stdout = subprocess.PIPE)
    process.wait() # Wait until dssp is done processing the file and calculating the Solvent Acessiblility  values
    input = open(dssp_out, 'r')    
    pdb['name'] = pdb_path.split('/')[-1][-10:-4] # split element 4 of the filecontent list (space-delimited). the 0th element is number of residues in pdb
    fileContents = input.readlines()   # This is a list containing each line of the input file as an element.
    pdb['nres'] = fileContents[3].split()[0] # split element 4 of the filecontent list (space-delimited). the 0th element is number of residues in pdb
    pdb['nchain'] = fileContents[3].split()[1]
    if int(pdb['nchain']) != 1:
        print ""
        print " WARNING: There are more than ONE chain in the input pdb file! "
        print " input pdb file:   ", pdb_path
        print " number of chains: ", pdb['nchain']
        print ""
        #raw_input('Press <Enter> to continue ...')
    pdb['sum_nssb']        = fileContents[3].split()[2]                         
    pdb['mean_nssb']   = float(fileContents[3].split()[2])/float(pdb['nres'])
    pdb['asa_pdb']     = fileContents[4].split()[0]
    #pdb['asa_normed'] = str(round(float(pdb['asa'])/float(pdb['nres']),4)) 
    pdb['sum_nhbon']   = fileContents[5].split()[0]
    pdb['mean_nhbon']  = float(fileContents[5].split()[0])/float(pdb['nres'])
    pdb['sum_nhbps']   = fileContents[6].split()[0]
    pdb['mean_nhbps']  = float(fileContents[6].split()[0])/float(pdb['nres'])
    pdb['sum_nhbas']   = fileContents[7].split()[0]
    pdb['mean_nhbas']  = float(fileContents[7].split()[0])/float(pdb['nres'])
    # Now calculate the mean and average of all hydrogen bonds energies, residue RSAs, and WCN based on CA coordinates:
### hbe = []   # A list containing all hydrogen bonds energies in the pdb file
### rsa = []   # A list containing all normalized residue ASA values in the pdb file
### crd = []   # A list containing all CA-atom WCN in the pdb file
### counter = 0
### for record in fileContents[25:len(fileContents)]:
###     AA = record[13]
###     if AA != ('!' or '*'):
###         counter += 1
###         if AA == ('X'):
###             print 'potential Histidine variant HID, HIE or HIP found in pdb file ', pdb['name'], 'residue number ', counter
###             AA = 'H'    # correct Histidine variants code to the standard code
###             #raw_input('Press <Enter> to continue ...')
###         hbe1 = float(record[38:51].split(',')[1])
###         hbe2 = float(record[50:62].split(',')[1])
###         hbe3 = float(record[61:73].split(',')[1])
###         hbe4 = float(record[72:84].split(',')[1])
###         hbe.append(numpy.mean([hbe1,hbe2,hbe3,hbe4]))  # This is the average of the four possible hydrogen bonds of a single residue. If a residue has no hydrogen bond, then all four are zero, so the mean is also zero.
###         #print record[13],counter,hbe1,hbe2,hbe3,hbe4,hbe[-1]
###         rsa.append(float(record[35:39])/residue_max_acc[AA])
###         #print AA,rsa[-1]        
###         #crd = numpy.array( [ float(record.split()[-3]) , float(record.split()[-2]) , float(record.split()[-1]) ] )
###         crd.append( numpy.array( [ float(record.split()[-3]) , float(record.split()[-2]) , float(record.split()[-1]) ] ) ) # stores CA atom triplet coordinates as individual elements of the list crd.
###         #if counter > 1 : print crd[-1], crd[-1] + crd[-2]
### if int(pdb['nres']) != (counter or len(hbe) or len(rsa) or len(wcn)) :
###     print 'Something is fishy about DSSP output'
###     print 'The calculated number of residues in the pdb file do not via two different methods!'
###     print ( 'counter = %4d , nres = %4d , len(hbe) = , len(rsa) = , len(wcn) = ' %( counter, int(pdb['nres']), len(hbe) , len(rsa), len(wcn) ) )
###     sys.exit('Program stopped. Check PDB file and DSSP output file for inconsistencies')
###     pdb['mhbe'] = str(numpy.mean(hbe))
###     pdb['vhbe'] = str(numpy.var(hbe))
###     pdb['mrsa'] = str(numpy.mean(rsa))
###     pdb['vrsa'] = str(numpy.var(rsa))
###     # Now calculate residue wcn :
###     wcn = []   # A list containing all CA-atom WCN in the pdb file
###     for i in range(len(crd)) :
###         sum_terms = 0.
###         for j in range(len(crd)) :
###             if i != j :
###                 sum_terms += 1./( (crd[i][0]-crd[j][0])**2 + (crd[i][1]-crd[j][1])**2 + (crd[i][2]-crd[j][2])**2 )
###         wcn.append(sum_terms)
###         #print i,j,crd[i], wcn[-1]
###     pdb['mwcn'] = str(numpy.mean(wcn))
###     pdb['vwcn'] = str(numpy.var(wcn))
    #print len(rsa),len(hbe), len(wcn)
    input.close() #Close the file with the dssp output file
    #os.remove('pdbOutput.txt') #Deletes the dssp output file 
    return (pdb) #Return the RSA values and the SAList

def main():
    if len( sys.argv ) != 4:
        print '''
 
Usage:'''
        print "     ", sys.argv[0], "<input pdb file>", "<output dssp file>", "<output summary file>", '\n'
        sys.exit('Program stopped.\n')
    else:
        pdb_path = sys.argv[1]  # path to the input pdb file
        dssp_out = sys.argv[2]  # path for the output dssp file to be written
        sum_out = sys.argv[3]   # summary file containing pdb structure property data

    pdb = get_pdb_props(pdb_path,dssp_out)

    # Now check if the output summary file alreadu exists. If so, then append data for the new pdb to the data already in the file.
    if os.path.isfile(sum_out):
       sum_out_file = open(sum_out,'a')
    else:
       sum_out_file = open(sum_out,'w')
       #sum_out_file.write('name' + 'nres' + '\t' + 'nchain' + '\t' + 'nssb' + '\t' + 'asa' + '\t' + 'nhb' + '\t' + 'nhbp' + '\t' + 'nhba' + '\t' + 'nah' + '\t' + 'nresah' + '\t' + 'nbs' + '\t' + '\t' + '\n')
###    for key in pdb:
###        sum_out_file.write(key + '\t')
###    sum_out_file.write('\n')
       sum_out_file.write( 'pdb' + '\t' + 'pdb_asa' + '\t' + 'nres' + '\t' + 'nchain' + '\t' \
                         + 'sum.nssb' + '\t' + 'mean.nssb' + '\t' + 'sum.nhbon' + '\t' + 'mean.nhbon' + '\t' \
                         + 'sum.nhbps' + '\t' + 'mean.nhbps' + '\t' + 'sum.nhbas' + '\t' + 'mean.nhbas' + '\n' )

    for key in pdb:
        sum_out_file.write(str(pdb[key]) + '\t')
    sum_out_file.write('\n')
    #sum_out_file.write(pdb.name + '\t' + pdb.nres + '\t' + pdb.nchain + '\t' + pdb.nssb + '\n')
	#for i, aa in enumerate(AA_List):
	#	rsa_outputfile.write(str(i+1) + '\t' + str(aa) + '\t' + str(RSA[i]) + '\n')
    
    return

if __name__ == "__main__":
   main()
