#!/opt/local/bin/python
#Descripion: This is a helper file that can be imported which some helper functions
#Most functions have an associated docstring. To print out the docstring for a function: do  my_func.__doc__
#	Ex. print back_var_helper.extract.__doc__

import re, os, math, string
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio import AlignIO
import subprocess
import numpy as np

_hydrogen=re.compile("[123 ]*H.*") 

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

class ChainSelector(object): 
	"""
	Adapted from the extract() function of Biopython 
	"""

	def __init__(self, chain_id, start, end, model_id=0): 
		self.chain_id=chain_id 
		self.start=start 
		self.end=end 
		self.model_id=0 

	def accept_model(self, model): 
		# model - only keep model 0 
		if model.get_id()==self.model_id: 
			return 1 
		return 0 

	def accept_chain(self, chain): 
		if chain.get_id()==self.chain_id: 
			return 1 
		return 0 
		
	def accept_residue(self, residue): 
		# residue - between start and end 
		hetatm_flag, resseq, icode=residue.get_id() 
		if hetatm_flag!=" ": 
			# skip HETATMS 
			return 0 
		if self.start<=resseq<=self.end: 
			return 1 
			return 0 

	def accept_atom(self, atom): 
		# atoms - get rid of hydrogens 
		name=atom.get_id() 
		if _hydrogen.match(name): 
			return 0 
		else: 
			return 1 

def extract(structure, chain_id, start, end, filename): 
	""" 
	Write out selected residues from a particular given pdb structure to filename. 
	
	Args: 
		structure: The four letter pdb code for the protein (this pdb must currently exist!)
		chain_id: The chain in the pdb you are trying to extract
		start: The beginning residue that you want to start your extraction from
		end:
		filename: The filename (ex. my_protein.pdb) of the pdb that you will extract the selected residues to  	
	Returns:
		A pdb file that contains only the selected residues from the original pdb structure given as an argument
		
	Example of usage
		p=PDBParser()
		s=p.get_structure('X', '1RUZ.pdb')
		extract(s, 'H', 40, 60, 'extracted.pdb')
 
	""" 
	sel=ChainSelector(chain_id, start, end) 
	io=PDBIO() 
	io.set_structure(structure) 
	io.save(filename, sel) 


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

def is_diverged(divergence_cutoff, seq, seqs_to_be_compared):
	"""
	This function takes a sequence and compares it to all the sequences in a list to determined if is diverged enough from the sequences in the list.


	Args:
		divergence_cutoff: A decimal (0 to 1.0) that represents how similar a given sequence can be to each sequence in the list being compared.
		seq: The sequence you want to compare to the list of sequences
		seqs_to_be_compared: A list of sequences that you are checking your sequence against for divergence. 

	Returns:
		A boolean (True or False). If True the seq being comparing is diverged by at least (1 - divergence_cutoff). 
		For example, if divergence_cutoff is 0.75, then the function will return True only if the sequence is at least 0.25 (or 25%) diverged
		from every sequence in the list seqs_to_compared.

	"""

	sequence_similar = True
	for sequence in seqs_to_be_compared: #This compares the seq to all of the sequences in the list AA by AA
		similarity = calculate_identity(seq, sequence, ignore_insert_states=True) 
		if(similarity>divergence_cutoff): #If sequence is too similar to the sequence compared return False immediately 
			sequence_similar = False
			return sequence_similar
		else:
			continue
	return sequence_similar 



def calculate_identity(seq1, seq2, ignore_insert_states=False):
	"""
	Calculates the sequence identity between two sequences.
	
	Args: 
		seq1: The first sequence (to be compared with the second)
		seq2: The second sequence (to be compared with the first)
		ignore_insert_states: A boolean controls whether insert states are to be ignored when calculating identity
	
	Returns:
		A decimal the sequence identity between two sequences	
	
	"""
	#print "IN FUNCTION"
	#print "seq1"
	#print seq1
	#print "seq2"
	#print seq2
	seq1_array = get_sequence_array(seq1) #Take each sequence and turn it into a list of characters
	seq2_array = get_sequence_array(seq2)
	similiar_count = 0
	aa_counts = 0
	smallest_length = len(min([seq1, seq2], key = len)) #Take the minimum in case the sequences are not of equal length
	for i in xrange(0, smallest_length): #For each character in the array
		aa1 = seq1[i] #Save the character for one site as aa1 for the first sequence and aa2 sort the second sequence
		aa2 = seq2[i]
		aa_equals = (aa1 == aa2)
		if ignore_insert_states and ( aa1.islower() or aa2.islower() ) :
			continue
		if aa1 == '-' or aa2 == '-':
			continue
		#print "aa_counts: " + str(aa_counts)
		aa_counts = aa_counts + 1
		if(aa1 == aa2):
			similiar_count = similiar_count + 1
	#print float(similiar_count), float(aa_counts)
	if(aa_counts == 0):
		identity = 0.0 #Could be the case that there are no amino acids found at corresponding positions
		return identity
	else:
		identity = float(similiar_count)/float(aa_counts)
		return identity


def align_sequences(seq_file, id, HMM, final_fasta_format = False):
	"""
	Aligns a file of sequences using a Hidden Markov (HMM) using HMMER (for example aligning pfam protein families)
	
	Args: 
		seq_file: The name of the file that will be aligned
		id: An identifier that describes the sequences to be aligned (for example the pfam protein family id). This is used for naming the output	
		HMM: The file that is used by HMMER to do the alignment of the sequences
		final_fasta_format: A boolean that determines where the final fasta that is returned should be have all characters (ex. amino acids) uppercased.
		If TRUE, the characters are all uppercased
		
	Returns:
		This function does not return anything but a fasta file is created with the aligned sequences.

	"""
	
	out = id + "_aligned.txt"
	out_fasta = id + "_aligned.fasta"
	align_string = "hmmalign -o " + out + " " + HMM + " " + seq_file #String representing a command line call to HMMER
	process = subprocess.call(align_string, shell=True) #Execute the command line call to HMMER to align the sequences
	print "Finished Aligning ", seq_file 
	
	#This block looks uses biopython's Bio.AlignIO to convert the stockholm file returned by HMMER into a fasta file
	input_handle = open(out, "rU")
	output_handle = open(out_fasta, "w")
	alignments = AlignIO.parse(input_handle, "stockholm")
	AlignIO.write(alignments, output_handle, "fasta")
	
	output_handle.close()
	input_handle.close()
	
	#Read the from the created fasta file into a list
	file = open(out_fasta, "r")
	fasta_lines = file.readlines()
	file.close()

	file = open(out_fasta, "w")
	#Writes the alignment back to a fasta (either as is with the insertion states lowercased or not depending on the value of final_fasta_format)
	for line in fasta_lines:
		if(final_fasta_format == True):
			upper_line = line.upper()
		else:
			upper_line = line
			#print line
		file.write(upper_line)	
	file.close()

def count_insert_states(seq):
	"""
	Counts the number of insert states in a sequence that was been aligned using a Hidden Markov Model (HMM)	
	
	Args:
		seq: A string representing a sequence whose insert states will be counted 
		
	Returns:
		The number of insert states in the sequence
	
	"""
	num_insert_states = 0
	seq_array = get_sequence_array(seq) #Take each sequence and turn it into a list of characters
	for aa in seq_array:
		if (aa.islower() and aa != "-"):
			num_insert_states = num_insert_states + 1
	return num_insert_states

def count_match_states(seq):
	"""
	Counts the number of match states in a sequence that was been aligned using a Hidden Markov Model (HMM)
	
	Args: 
		seq: A string representing a sequence whose match states will be counted
	
	Returns:
		The number of match states in the sequence
	
	"""
	num_match_states = 0
	seq_array = get_sequence_array(seq) #Take each sequence and turn it into a list of characters
	for aa in seq_array:
		if (aa.isupper() and aa != "-"):
			num_match_states = num_match_states + 1
	return num_match_states
	
def calculate_occupancy(fasta_alignment):
	"""
	Calculates the occupancy (the percent number of sequences that have a amino acid) of each column in an alignment

	Args: 
		fasta_alignment: The filename of a fasta file containing an alignment	
	
	Returns:
		A list called occupancy, where each element is the percent occupancy of each column
		For example, occupancy[0], would be the percent occupancy of the first column in the alignment
		If occupancy[0] = 0.75, the 75 percent of the sequences in the alignment have an amino acid at the first position
	
	"""

	seq_list = []
	occupancy = []
	[seqs, headers] = get_sequences(fasta_alignment)
	for s in seqs:
		seq_list.append(get_sequence_array(s))

	alignment_array = np.array(seq_list)
	rows,cols = alignment_array.shape

	for j in xrange(0, cols):
		site = list(alignment_array[:,j])
		#print site
		num_sites_occupied = rows - site.count("-")
		percent_occupancy = float(num_sites_occupied)/float(rows)
		occupancy.append(percent_occupancy)
	return occupancy
	
def is_occupied_enough(seq, alignment_occupancy, max_unoccupied):
	"""
	Determines whether a given sequence is occupied
	
	Args: 
		seq: A string representing the sequence to be assessed for occupancy
		alignment_occupancy: A list whose elements represent the occupancy of an alignments columns
		max_unoccupied: The max number of columns that can be occupied
		
	Returns:
		occupied: A boolean that determines whether a sequence has enough amino acids at highly occupied positions
		
	"""
	num_unoccupied = 0
	occupied = True
	seq_array = get_sequence_array(seq)
	for i in xrange(0, len(seq_array)):
		if(seq_array[i]== "-" and (alignment_occupancy[i]>= 0.6)):
			num_unoccupied = num_unoccupied + 1
	if (num_unoccupied > max_unoccupied):
		occupied = False

	return occupied
	
def convert_to_one_letter_code(residue_list):
	"""
	Takes a list of amino acids in letter-code and changes it to one-letter code.
	
	Args:
		residue_list: A list with an amino acid sequence represented in three-letter code
	
	Returns:
		one_letter_residue_list: A list whose elements are the amino acids in one-letter code
	
	"""
	one_letter_residue_list = []
	
	for res in residue_list:
		abbrev = resdict[res]
		one_letter_residue_list.append(abbrev)
	return one_letter_residue_list