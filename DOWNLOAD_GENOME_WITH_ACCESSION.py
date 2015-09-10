#!/usr/bin/python
import sys, os, re, getopt, glob, numpy as np
import timeit
from Bio import SeqIO

start = timeit.default_timer()

GENOME_NAMES = open("INPUT.tsv", "r")
output = open("ACCOUNT of DOWNLOADED GENOMES and ACCESSIONS.out","w")
output.write("GENOME\tACCESSION(S)\n")

HOLD_DICT = {}
ACCESSION_BIN = []

# READ GENOME NAMES
for FILE in GENOME_NAMES:

	FILE = FILE.strip()
	FILE = FILE.split("\t")
	NAME = FILE[0]
	NAME = NAME.strip()

	ACCESSION = FILE[1]
	ACCESSION = ACCESSION.strip()

	#Perform esearch
	os.system(' '.join([
		"esearch -db genome -query",
		ACCESSION,
		"| elink -target nuccore | efetch -format docsum | xtract -Pattern DocumentSummary -element Title -element Extra > FOO.txt"
	]))

	#Parse output
	readme = open('FOO.txt', 'r')
	chooseme = []

	count = 1

	for line in readme:
		line = line.strip()
		ENTRY = line.split("\t")
		
		chooseme.append([str(count)+"    "+ENTRY[0], ENTRY[1]])
		count = count + 1

	if np.size(chooseme) > 4:
		HOLD_DICT[NAME] = chooseme

	else:
		FOO = chooseme[0][1]
		FOO = FOO.split("|")
		DOWNLOAD = FOO[3]

		print "DOWNLOADING: "+DOWNLOAD
		output.write(NAME+"\t"+DOWNLOAD+"\n")

		NAME = NAME.replace(" ", "_")

		os.system(' '.join([
			"esearch -db genome -query",
			DOWNLOAD,
			"| elink -target assembly | elink -target nuccore | efetch -format fasta >>",
			NAME+".genome.fa"
		]))

for key, value in HOLD_DICT.iteritems():
	NAME = key
	chooseme = value

	GO_ON = "N"
	print "\nChoose all genome files you wish to download under the organism name: "+NAME+"\n"
	for ENTRY in chooseme:
		print ENTRY[0]

	while GO_ON == "N":
		print "\nPlease input the number(s) delimited by comma. For example: 1,3,4,7,8\n" 
		NUMBERS = raw_input()
		NUMBERS = NUMBERS.split(",")
		
		print "Are these the numbers you selected: "+str(NUMBERS)+"\n"
		print "Y/N?"
		GO_ON = raw_input()

	for i in NUMBERS:
		FOO = chooseme[int(i)-1][1]
		FOO = FOO.split("|")
		DOWNLOAD = FOO[3]

		print "DOWNLOADING: "+DOWNLOAD
		output.write(NAME+"\t"+DOWNLOAD+"\n")

		NAME = NAME.replace(" ", "_")

		os.system(' '.join([
			"esearch -db genome -query",
			DOWNLOAD,
			"| elink -target assembly | elink -target nuccore | efetch -format fasta >>",
			NAME+".genome.fa"
		]))

output.close()

#De-replicate all genome files based on ID and Sequence
error = open("ACCESSION_TRACKER.out","w")

for FILE in GENOME_NAMES:
	FILE = FILE.strip()
	FILE = FILE.split("\t")
	NAME = FILE[0]
	NAME = NAME.strip()
	OUT = NAME.replace(" ", "_")
	OUT = OUT+".genome.fa"

	print "\n\n"+OUT

	Seq_LIST = []

        #Initialize Dict
        Unique_Dictionary = {}

        # Open output file
        genome_output = open(re.sub(".fa","",OUT)+".uniq.fa", "w")

        # Create Unique Dictionary entry for each sequence
        for record in SeqIO.parse(open(OUT, "rU"), "fasta") :
		error.write(NAME+"\t"+str(record.id)+"\tALL\n")

		# Remove identical sequences
		if not str(record.seq[:100]+record.seq[100:]) in Seq_LIST:
	                Unique_Dictionary[record.id] = str(record.seq)
			Seq_LIST.append(str(record.seq[:100]+record.seq[100:]))
			error.write(NAME+"\t"+str(record.id)+"\tKEPT\n")

	print len(Seq_LIST)

        # Print to File
        for key, value in Unique_Dictionary.iteritems():
                id = key
                seq = value

                genome_output.write(">"+id+"\n"+seq+"\n")

        genome_output.close()

	os.system("rm "+NAME+".genome.fa")


error.close()
os.system("rm FOO.txt")


stop = timeit.default_timer()
print "This operation took " + str(stop - start) + " seconds."
