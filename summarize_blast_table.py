#!/usr/bin/python
import sys, os, re, getopt, numpy as np
import timeit

start = timeit.default_timer()

Usage = """
REQUIRED ARGUMENTS:

                -i      blast output
		-s	blasted sequences (fasta "\\n"-stripped)
		-o	output file name

Utility: To take blast output from the below command and organize into table. Further, you must provide the fasta file used for the query and the script will calculate the total coverage of the blast hit.

blastn -query my_query.txt -db databasename.db -out ouput.txt -outfmt 7 -max_target_seqs 1



"""

if len(sys.argv)<2:
        print Usage
        sys.exit()

header = "query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, length of query, % overlap of full query"

# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"i:s:o:")

###############################
# o == option
# a == argument passed to the o
###############################
for o, a in myopts:
    if o == '-i':
        BLAST= a
    if o == '-s':
        SEQ= a
    if o == '-o':
        OUTPUT= a

#Open file to write
output = open(OUTPUT+".tsv", "w")

#Write Table Header
header = re.sub(",","\t", header)
output.write(header+"\n")

#Calculate length of Input sequence
SEQ_LENGTH = {}
infile = open(SEQ, "r")
for line in infile:
	if line.startswith(">"):
		line = re.sub(">","",line)
		line = line.strip("\n\r")
		SEQ_LENGTH[line] = len(next(infile))	

#Extract OTU List According to Desired Cut-off
for line in open(BLAST, "r"):
	if not line.startswith("#"):
		line = line.strip("\n\r")
		split = line.split()
		TOTAL_LENGTH = SEQ_LENGTH[split[0]]
		ALIGNMENT_LENGTH = int(split[3])
		PERC_LENGTH = float(ALIGNMENT_LENGTH) / TOTAL_LENGTH
		output.write(line+"\t"+str(TOTAL_LENGTH)+"\t"+str(round(PERC_LENGTH, 3)*100)+"\n")

output.close()

stop = timeit.default_timer()

print stop - start

