import sys
import re
import ddG_var_helper
import os, subprocess
from Bio import PDB
from Bio import AlignIO
from Bio.PDB.Polypeptide import PPBuilder 
'''
This is a file that creates the runfiles, paramlist, and pdb lists needed to run FoldX for the proteins
'''

pdb_dict = {"1RD8_X": "HP", "2FP7_B" : "WNPB" , "2JLY_A": "DPH" , "2Z83_A" : "JEHN" ,
			"3GOL_A" : "HCP" , "3GSZ_A": "HCP" , "3I5K_A": "HCP" , "3LYF_A": "RVFVNP" , 
			"4AQF_B": "CCHFN" , "4GHA_A": "MRNABD" , "4IRY_A": "INP" }

#Description: This script creates the runfiles, paramlists and lists used in Position Scan mutagenesis. 
def create_residue_lists(pdb, chain):
	print pdb, chain
	res_lists = []
	aa_list = ""
	try:
		ppb = PPBuilder()	
		p=PDB.PDBParser(QUIET = True)
		s=p.get_structure('X', "../foldx_setup/repaired_pdbs/RepairPDB_" + pdb + "_" + chain + ".pdb")
		pp = ppb.build_peptides(s)[0]
		seq = str(pp.get_sequence())
		
		#Get the structures
		ref_struct = s[0]
		ref_chain = ref_struct[chain]
		ref_residues = []
		ref_res_nums = []
		for res in ref_chain:
			ref_residues.append( res.resname )
			ref_res_nums.append(res.id[1])				
	except KeyError:
		print "Something is wrong with the mapped residues"
		#continue
	seq_array = ddG_var_helper.convert_to_one_letter_code(ref_residues)
	
	count = 0
	i = 0
	while(i < len(seq_array)):	
		if (i < 101): #These statements split up the sequence into 100 aa chunks for performing a Position Scan
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			if(i == len(seq_array) - 1):
				aa_list = aa_list + ";"	
				res_lists.append(aa_list)
		elif(i == 101):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			aa_list = aa_list + ";"	
			res_lists.append(aa_list)
			aa_list = ""
		elif(i < 201 and i>101):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			if(i == len(seq_array) - 1):
				aa_list = aa_list + ";"	
				res_lists.append(aa_list)
		elif(i == 201):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			aa_list = aa_list + ";"	
			res_lists.append(aa_list)
			aa_list = ""
		elif(i < 301 and i>201):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			if(i == len(seq_array) - 1):
				aa_list = aa_list + ";"	
				res_lists.append(aa_list)
		elif(i == 301):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			aa_list = aa_list + ";"	
			res_lists.append(aa_list)
			aa_list = ""
		elif(i < 401 and i>301):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			if(i == len(seq_array) - 1):
				aa_list = aa_list + ";"	
				res_lists.append(aa_list)
		elif(i == 401):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			aa_list = aa_list + ";"	
			res_lists.append(aa_list)
			aa_list = ""
		elif(i < 501 and i>401):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			if(i == len(seq_array) - 1):
				aa_list = aa_list + ";"	
				res_lists.append(aa_list)			
		elif(i == 501):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			aa_list = aa_list + ";"	
			res_lists.append(aa_list)
			aa_list = ""
		if (i < 601 and i>501):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			if(i == len(seq_array) - 1):
				aa_list = aa_list + ";"	
				res_lists.append(aa_list)
		elif(i == 601):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			aa_list = aa_list + ";"	
			res_lists.append(aa_list)
			aa_list = ""
		if(i < 701 and i>601):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			if(i == len(seq_array) - 1):
				aa_list = aa_list + ";"	
				res_lists.append(aa_list)
		elif(i == 701):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			aa_list = aa_list + ";"	
			res_lists.append(aa_list)
			aa_list = ""
		if (i < 801 and i>701):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			if(i == len(seq_array) - 1):
				aa_list = aa_list + ";"	
				res_lists.append(aa_list)
		elif(i == 801):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			aa_list = aa_list + ";"	
			res_lists.append(aa_list)
			aa_list = ""
		if (i < 901 and i>801):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			if(i == len(seq_array) - 1):
				aa_list = aa_list + ";"	
				res_lists.append(aa_list)
		elif(i == 901):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			aa_list = aa_list + ";"	
			res_lists.append(aa_list)
			aa_list = ""
		elif(i < 1001 and i>901):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			if(i == len(seq_array) - 1):
				aa_list = aa_list + ";"	
				res_lists.append(aa_list)
		elif(i == 1001):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			aa_list = aa_list + ";"	
			res_lists.append(aa_list)
			aa_list = ""
		elif(i < 1101 and i>1001):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			if(i == len(seq_array) - 1):
				aa_list = aa_list + ";"	
				res_lists.append(aa_list)
		elif(i == 1101):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			aa_list = aa_list + ";"	
			res_lists.append(aa_list)
			aa_list = ""
		elif(i < 1201 and i>1101):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			if(i == len(seq_array) - 1):
				aa_list = aa_list + ";"	
				res_lists.append(aa_list)
		elif(i == 1201):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			aa_list = aa_list + ";"	
			res_lists.append(aa_list)
			aa_list = ""
		elif(i < 1301 and i>1201):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			if(i == len(seq_array) - 1):
				aa_list = aa_list + ";"	
				res_lists.append(aa_list)
		elif(i == 1301):
			aa_list = aa_list + ","
			res_num = ref_res_nums[i]
			aa_list = aa_list + seq_array[i] + chain + str(res_num) + "a"
			aa_list = aa_list + ";"	
			res_lists.append(aa_list)
		i = i + 1
	return res_lists

def create_runfiles(pdb, chain): #This function creates a runfile used for performing a Position Scan in FoldX
	runfiles = []
	residue_lists  = create_residue_lists(pdb, chain)
	for i in xrange(0, len(residue_lists)):
		filename = "../foldx_setup/runfiles/runfile_" + pdb + "_" + chain + "_" + str(i+1) + ".txt"
		runfile = "runfile_" + pdb + "_" + chain + "_" + str(i+1) + ".txt"
		runfiles.append(runfile)
		listname = "pdb_list_" + pdb + "_" + chain + ".txt"	
		out = open(filename,"w") 
		res_string = residue_lists[i]
		out.write("<TITLE>FOLDX_runscript;" + "\n")
		out.write("<JOBSTART>#;" + "\n")
		out.write("<PDBS>#;" + "\n")
		out.write("<BATCH>" + listname + ";" + "\n")
		out.write("<COMMANDS>FOLDX_commandfile;" + "\n")
		out.write("<PositionScan>#" + res_string + "\n")
		out.write("<END>#;" + "\n")
		out.write("<OPTIONS>FOLDX_optionfile;" + "\n")
		out.write("<Temperature>298;" + "\n")
		out.write("<R>#;" + "\n")
		out.write("<pH>7;" + "\n")
		out.write("<IonStrength>0.050;" + "\n")
		out.write("<water>-CRYSTAL;" + "\n")
		out.write("<metal>-CRYSTAL;" + "\n")
		out.write("<VdWDesign>2;" + "\n")
		out.write("<OutPDB>false;" + "\n")
		out.write("<pdb_hydrogens>false;" + "\n")
		out.write("<complex_with_DNA> true;" + "\n")
		out.write("<END>#;" + "\n")
		out.write("<JOBEND>#;" + "\n")
		out.write("<ENDFILE>#;")
		out.close()
	return runfiles	

#This function creates a paramlist used mutating the using launcher on a HPC cluster
def create_paramlist(pdb, chain, runfiles):	
	i = 1
	out = open("../foldx_setup/paramlist_mutate_pdbs", "a+") 
	for file in runfiles:
		pdb_file = "RepairPDB_" + pdb + "_" + chain + ".pdb"
		out.write("./foldx_mutate.sh " + pdb + "_" + chain + " " + pdb + "_" + chain + "_" + str(i) + " " + pdb_file + " " + file + "\n")
		i = i + 1 
	out.close()

#Main function
def main(argv = sys.argv):
	pdbs = []
	chains = []
	for key in pdb_dict.keys():
		p = re.split("_", key)
		pdbs.append(p[0])
		chains.append(p[1])
	
	subprocess.call("mkdir -p ../foldx_setup/pdb_lists", shell = True) #Open directories for the files
	subprocess.call("mkdir -p ../foldx_setup/runfiles", shell = True)
	for i in xrange(len(pdbs)):	 #For each pdb in the dataset
		filename = "../foldx_setup/pdb_lists/pdb_list_" + pdbs[i] + "_" + chains[i] + ".txt"
		out = open(filename, "w")
		out.write("RepairPDB_" + pdbs[i] + "_" + chains[i] + ".pdb")
		out.close()
		runfiles = create_runfiles(pdbs[i], chains[i]) #Creates the runfiles needed for foldx
		create_paramlist(pdbs[i], chains[i], runfiles) #Creates the paramlist 
	
if __name__ == "__main__":
	main(sys.argv)