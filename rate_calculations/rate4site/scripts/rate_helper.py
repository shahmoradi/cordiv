from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio import AlignIO
from Bio import PDB

'''
A helper file needed to run the scripts in the evolutionary_rates directory 
'''

def get_sequence_array(seq):
	"""
	Turns a string representing a string into a list of characters.
	
	Args:
		seq: A string representing a sequence (ex. nucleotide, protein)
		
	Returns:
		seq_array: A list of single letter elements representing the sequence (ex. each amino acid in a protein is an element in the list)
	
	"""

	seq_array = []
	seq_length = len(seq)
	count = 0
	while count < seq_length:
		aa = seq[count]
		seq_array.append(aa)
		count = count + 1
	return seq_array

def write_seq_to_file(seqs, headers, out_fasta):
	"""
	Writes a group of sequences to a file
	
	Args:
		seqs: A list of sequences
		headers: A list of descriptors (often fasta headers > included)
		out_fasta: The name of the file that the sequences will be written to
		
	Returns:
		out_fasta: A string that is the filename that was given to the function
	
	"""
	file = open(out_fasta, "w")
	j = 0
	for seq in seqs:
		file.write(headers[j]+ "\n")
		file.write(seq + "\n")
		j = j+1
	file.close()
	return out_fasta

def get_sequences(file):
	"""
	This function takes a file with a bunch of sequences in fasta format and returns a list of the sequences. 
	
	Args: 
		file: The name of the fasta file of the given sequences
		
	Returns:
		Returns TWO lists
		It returns: all_sequences, headers
		all_sequences is a list of every sequence in the fasta file
		headers is a list of every header for each sequence in the fasta file
		The headers and the sequences are in corresponding elements in the two lists 
		(so the first element in the header list is the header for the first sequence in the sequence list.)				
	"""
		
	all_sequences = []  
	natural_sequences = []
	designed_sequences = []
	headers = []

	file_data = open(file, "r")
	seq_data = file_data.readlines()
	file_data.close()
	
	#print seq_data
	#This block of code basically just removes all of the headers from the array of sequence information.
	#It creates all_sequences, which is a list of all the sequences from the natural alignment file.
	string  = ''
	finished_sequence = ''
	for sequence in seq_data:
		if sequence[0] == '>': #If it is a header append the last sequence that was processed
			headers.append(sequence.rstrip())
			if(string != ''):
				finished_sequence = string #Just in case the sequence was translated from Stockholm format
				all_sequences.append(finished_sequence) #Appends the sequence into the array full of aligned sequences
				string = '' #Empty the string that represents the sequence (you don't want to append the old with the new)
		else:       
			string = string + sequence.rstrip("\n") #Strip the new line that is at the end of each sequence in alignment
	all_sequences.append(string)
	num_sequences = len(all_sequences)
	return all_sequences, headers  
	
def convert_alignment_format(input_alignment, in_format, out_format, out_file, all_caps = False):
	"""
	Convert alignment format
	
	Args: 
		in_format: The name of the file containing the fasta alignment
		out_format
		out_file
		final_phylip_format: A boolean that determines where the final fasta that is returned should be have all characters (ex. amino acids) uppercased.
		If TRUE, the characters are all uppercased
		
	Returns:
		out_file: The name of the fasta file that is created

	"""
	#fileparts = re.split("_", input_alignment)
	#id = fileparts[0]
	#out_fasta = id + "_aligned.fasta"

	#This block looks uses biopython's Bio.AlignIO to convert the fasta file into phylip 
	input_handle = open(input_alignment, "rU")
	output_handle = open(out_file, "w")
	alignment = AlignIO.parse(input_handle, in_format)
	AlignIO.write(alignment, output_handle, out_format)
	
	output_handle.close()
	input_handle.close()
	
	#Read the from the created fasta file into a list
	file = open(out_file, "r")
	f_lines = file.readlines()
	file.close()

	file = open(out_file, "w")
	#Writes the alignment back to a fasta (either as is with the insertion states lowercased or not depending on the value of final_fasta_format)
	for line in f_lines:
		if(all_caps == True):
			upper_line = line.upper()
		else:
			upper_line = line
		file.write(upper_line)	
	file.close()
	
def write_seq_to_file(seqs, headers, out_fasta):
	"""
	Writes a group of sequences to a file
	
	Args:
		seqs: A list of sequences
		headers: A list of descriptors (often fasta headers > included)
		out_fasta: The name of the file that the sequences will be written to
		
	Returns:
		out_fasta: A string that is the filename that was given to the function
	
	"""
	file = open(out_fasta, "w")
	j = 0
	for seq in seqs:
		file.write(headers[j]+ "\n")
		file.write(seq + "\n")
		j = j+1
	file.close()
	return out_fasta

def get_pdb_info(data_file):
	"""
	Extracts the list of the pdbs and chains from the Huang et al dataset
	
	Args:
		data_files: The location of the Huang et al datafile (Table S1)
	"""

	pdbs = []
	chains = []
	pdb_file = open(data_file, "r") #Open the list of pdbs to be downloaded
	pdb_file_data = pdb_file.readlines()
	pdb_file_data.pop(0)
	
	for entry in pdb_file_data: #These lines read in the details from the Pfam Database File 
		pdbs.append(entry[1:5]) #Append the pdb names to a list
		chains.append(entry[8]) #Append the chain names to a list		
	return pdbs, chains