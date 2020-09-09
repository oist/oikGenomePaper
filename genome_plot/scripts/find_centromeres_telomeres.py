#!/bin/python3

import binascii
import gzip
from tqdm import tqdm
import multiprocessing
import main_statistics
import os


"""
Simply checks if a file is GZIP compressed or not. Instead of hacking away with file
extensions, instead this checks the first two bytes of a file after opening. For all
GZIP files, the first two bytes of the file are "1f8b", which will not likely be the
case if it's not a compressed file.
"""
def is_gzip_file(filePath):
	with open(filePath, "rb") as test_file:
		return binascii.hexlify(test_file.read(2)) == b'1f8b'


"""
Simple FASTA reader. A FASTA parser can be written as a generator instead, but since
it passes through a few if/etc. logical separators, so this is easier to read.
"""


def read_fasta_file(fastaFile, minContigLength=0):
	fasta_dictionary = {}
	if is_gzip_file(fastaFile):
		num_lines = sum(1 for line in gzip.open(fastaFile, "r"))
		with gzip.open(fastaFile, "rt") as fasta_handler:
			tqdm_obj = tqdm(fasta_handler, total=num_lines, desc="Decompressing and reading FASTA file")
			for line in tqdm_obj:
				line = line.strip()
				if line.startswith(">"):
					header = line.split(">")[1]
					fasta_dictionary[header] = []
				else:
					fasta_dictionary[header].append(line)
	else:
		num_lines = sum(1 for line in open(fastaFile, "r"))
		with open(fastaFile, "r") as fasta_handler:
			tqdm_obj = tqdm(fasta_handler, total=num_lines, desc="Reading FASTA file")
			for line in tqdm_obj:
				line = line.strip()
				if line.startswith(">"):
					header = line.split(">")[1]
					fasta_dictionary[header] = []
				else:
					fasta_dictionary[header].append(line)
	fasta_dictionary = {header : "".join(fasta_dictionary[header]) for header in fasta_dictionary if len(fasta_dictionary[header]) > minContigLength}
	return(fasta_dictionary)


def find_repeat_matches(fasta_dict, repeat):
	results = {}
	repeat_re = re.compile(repeat)
	for fasta_key in fasta_dict:
		#results[fasta_key] = []
		results[fasta_key] = {repeat:list()}
		sequence = fasta_dict[fasta_key]
		for match in repeat_re.finditer(sequence):
			#results[fasta_key].append((match.start(), match.group()))
			results[fasta_key][repeat].append(match.start())
	return(results)

genome = read_fasta_file('../inputs/I69-5.masked.fa')
telomeres = find_repeat_matches(genome, 'TTAGGG')
#centromeres = find_repeat_matches(genome, '')

