#!/usr/bin/python
import cor_div_helper as ch
import numpy as np
import subprocess
import sys

'''
Last Updated: November 12, 2014
Description: This is a script that calculates entropy values for the viral sequence alignments
'''

aa_list = ['G', 'A', 'L', 'V', 'I', 'P', 'R', 'T', 'S', 'C', 'M', 'K', 'E', 'Q', 'D', 'N', 'W', 'Y', 'F', 'H'] #List of amino acids

def calculate_aa_frequencies(list_of_aligned_seqs): 
	'''
	This function calculates the amino acid frequencies for each site in an alignment
	Args: 
		list_of_aligned_seqs: A list of sequences that have been aligned
	
	Returns:
		all_freqs: A list of lists where each list are the amino acid frequencies for a site
	
	'''
	all_freqs = []
	aligned_seqs = []
	for seq in list_of_aligned_seqs:
		aligned_seqs.append(ch.get_sequence_array(seq)) #Turn the lists of sequences into a list of arrays (where each array represents a sequence 
	seq_alignment = np.array(aligned_seqs)
	rows, cols = seq_alignment.shape

	for i in xrange(0, cols): #For each column in the alignment numpy array (each site in the alignment)
		#print "site: ", i
		site_freqs = []
		site = seq_alignment[:,i].tolist()
		num_aas_at_site = len(site) - site.count('-') #Count the number of amino acids at the site
		for aa in aa_list: #For each of the 20 amino acid types
			aa_freq = float(site.count(aa))/float(num_aas_at_site) #Calculate the AA's frequences
			site_freqs.append(aa_freq) #Append it to the list of frequencies for that site
		all_freqs.append(site_freqs) #Append it to the list of lists for the frequencies
	return all_freqs

def calculate_entropy(data_array):
	'''
	Calculates the entropy of list of lists where each list represents a site
	
	Args:
		data_array: A 2D numpy array where the rows are sites and the columns represent a value for each of the 20 amino acids
			(The values could be its frequency at that site or another meaningful calculated quantity)
	Returns:
		entropy_values: The entropy value for each site	
	'''
	
	entropy_values = []
	num_residues,num_AA = data_array.shape  
	for i in xrange(0, num_residues): #For each site
		probs_values = filter(lambda a: a != 0.0, data_array[i] )
		prob_sum = sum(probs_values)
		#print prob_sum
		entropy_number = 0.0
		for j in xrange(0,len(probs_values)): #Calculate the entropy (Can look up this formula, just the native entropy)
			value = (float(probs_values[j])*np.log(float(probs_values[j])))
			entropy_number = entropy_number + value
			#print -entropy_number
		if (entropy_number == 0.0):
			entropy_values.append(0.0) #Append entropy value to entropy at sites. 
		else:
			entropy_values.append(-entropy_number) #Append entropy value to entropy at sites. 
	return entropy_values

def main(argv = sys.argv):
	process = subprocess.call("mkdir entropies", shell = "True")
	alignments = ["CCHFN", "DPH", "HCP", "HCP", "HCP", "HCP", "HP", "INP", "JEHN","JEHN", "MRNABD", "RVFVNP", "WNPB"]
	maps = ["CCHFN_4AQF_B", "DPH_2JLY_A", "HCP_3GOL_A", "HCP_3GSZ_A", "HCP_3I5K_A", "HP_1RD8_AB", "INP_4IRY_A", "JEHN_2JLY_A", "JEHN_2Z83_A", "MRNABD_4GHA_A", "RVFVNP_3LYF_A", "WNPB_2FP7_B" ]
	alignment_counter = 0
	for m in maps: #For each viral alignment
		freq_list = []
		alignment_file = "../alignments/" + alignments[alignment_counter] + ".fasta"
		map_file = "../../pdb_maps/" + m + ".map"
		entropy_file = open("entropies/" + m + "_entropies.csv", "w")
		entropy_file.write("alignment_pos,pdb_pos,entropy\n")
		aligned_seqs, headers = ch.get_sequences(alignment_file)
		aa_freqs = calculate_aa_frequencies(aligned_seqs)
		mapped_residues = np.genfromtxt(map_file, skip_header = 1, dtype = "string", usecols = (2)) #Get the mapping of the sites to the alignment
		for site in aa_freqs: #Calculate the frequencies of each amino acid of each site in the alignment
			freq_list.append(ch.get_sequence_array(site))
		freq_array = np.array(freq_list)
		entropies = calculate_entropy(freq_array) #calculate the entropy values
		entropy_counter = 0
		#This prints out the entropy values to a file
		for i in xrange(0, len(mapped_residues)): 
			if(mapped_residues[i] == "NA"): 
				continue #If the site in the alignment is NA meaning it is not in the pdb skip and do not write
			else: #else write the entropy values and associated alignment and pdb positions
				entropy_file.write(str(i) + "," + mapped_residues[i] + "," + str(entropies[entropy_counter]) + "\n")
			entropy_counter = entropy_counter + 1
		entropy_file.close()
		alignment_counter = alignment_counter + 1
if __name__ == "__main__":
	main(sys.argv)