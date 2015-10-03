#!/usr/bin/python
import sys, os, re, getopt, glob, numpy as np
import timeit

start = timeit.default_timer()

'''usage: ./FIND_GENOME_ACCESSION_NCBI.py GENOMES.list

List of all genome names for which ACCESSIONS are desired'''

FILE = sys.argv[1]
GENOME_NAMES = open(FILE, "r")
output = open("ACCESSION.out","w")

# READ GENOME NAMES
for NAME in GENOME_NAMES:

	NAME = NAME.strip()

	#Perform esearch
	os.system(' '.join([
		"esearch -db genome -query",
		"\""+NAME+"\"",
		"| efetch -format docsum | grep \"<Project_Accession>\" > FOO.txt"
	]))

	#Parse output
	with open('FOO.txt', 'r') as f:
		parse = f.readline()
	parse = parse.strip()
	parse = re.split("<|>",parse)

	if np.size(parse) > 1:
		ACCESSION = parse[2]

	else:
		ACCESSION = "Not Found"

	#Write
	output.write(NAME+"\t"+ACCESSION+"\n")


stop = timeit.default_timer()

os.system("rm FOO.txt")

print "This operation took " + str(stop - start) + " seconds."
