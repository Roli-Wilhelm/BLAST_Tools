#!/usr/bin/python
import sys, os, re, getopt, glob, numpy as np
import timeit, datetime
from Bio import SeqIO

start = timeit.default_timer()
input = sys.argv[1]
OUTPUT = sys.argv[2]

if OUTPUT:
        if os.path.exists('./' + OUTPUT):
                print "\n-- Output Folder Exists - Caution: Files May Be Over-written --"
        else:
                os.mkdir(OUTPUT)

# OPEN INPUT FILE
GENOME_NAMES = open(input, "r")
GENOME_NAMES_DICT = {}

# OPEN OUTPUT FILES
output = open(OUTPUT+"/ACCOUNT of DOWNLOADED GENOMES and ACCESSIONS.out","w")
output.write("GENOME\tACCESSION(S)\tSTATUS\n")

NCBI_list = open(OUTPUT+"/ALL_GENOMES_DISPLAYED_FOR_USER_INPUT.txt","w")
NCBI_list.write("Date of Processing: "+str(datetime.date.today())+"\n")

error = open(OUTPUT+"/ACCESSION_TRACKER.out","w")
error.write("GENOME\tACCESSION(S)\tSTATUS\n")

# INITIALIZE DICTIONARIES | LISTS
HOLD_DICT = {}
ACCESSION_BIN = []
MULTIPLE_GENOME_LIST = []
SINGLE_DOWNLOADS = []
SKIPPED_GENOMES = []

# READ GENOME NAMES
for FILE in GENOME_NAMES:

	FILE = FILE.strip()
	FILE = FILE.split("\t")
	NAME = FILE[0]
	NAME = NAME.strip()

	ACCESSION = FILE[1]
	ACCESSION = ACCESSION.strip()

	SUBBED_NAME = NAME.replace(" ", "_")
	GENOME_NAMES_DICT[SUBBED_NAME] = FILE

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
		# CATCH EMPTY HITS
		try:
	
			FOO = chooseme[0][1]
			FOO = FOO.split("|")
			DOWNLOAD = FOO[3]

			# Some Project Accessions Appear Unique, but Refer to the same Nucleotide Accession
			# Thus occassionally the same genome is downloaded.
			if DOWNLOAD in SINGLE_DOWNLOADS:
				print "SKIPPED DOWNLOADING: "+DOWNLOAD
				output.write(NAME+"\t"+DOWNLOAD+"\tSKIPPED\n")
				SKIPPED_GENOMES.append(SUBBED_NAME)

			else:
				print "DOWNLOADING: "+DOWNLOAD
				output.write(NAME+"\t"+DOWNLOAD+"\tDOWNLOADED\n")

				# KEEP TRACK OF ALL PREVIOUSLY DOWNLOADED
				SINGLE_DOWNLOADS.append(DOWNLOAD)
	
				NAME = NAME.replace(" ", "_")

				os.system(' '.join([
					"esearch -db genome -query",
					DOWNLOAD,
					"| elink -target assembly | elink -target nuccore | efetch -format fasta >>",
					OUTPUT+"/"+NAME+".genome.fa"
				]))

		# Some entries in NCBI have the actual data wharehoused elsewhere, like JGI etc.
		except IndexError:
			output.write(NAME+"\tSTORED_IN_EXTERNAL_DATABASE\n")


for key, value in HOLD_DICT.iteritems():
	NAME = key
	chooseme = value

	GO_ON = "N"

	# PRINT AVAILABLE GENOMES TO SCREEN
	Number_of_entries = int(np.size(chooseme))/2

	if Number_of_entries > 10:
		count = 0
		
		while count < Number_of_entries:
			for ENTRY in chooseme[count:count+15]:
				print ENTRY[0]
				NCBI_list.write(ENTRY[0]+"\n")

			if Number_of_entries-(count+15) > 0:
				print "There are "+str(Number_of_entries-(count+15))+" more entries to display.\nPress ENTER to list the next set."

			MORE = raw_input()
			count = count+15

	else:
		for ENTRY in chooseme:
			print ENTRY[0]
			NCBI_list.write(ENTRY[0]+"\n")

	# LET USER DECIDE TO DOWNLOAD THEM ALL or SELECT A SUBSET
	while GO_ON == "N":
		print "\nDo you wish to download all the listed genomes (1), or select individual files to download (2)?"

		try:
			FORK = int(raw_input())


			if FORK == 1:
				NUMBERS = range(0, Number_of_entries+1)
				GO_ON = "Y"

			if FORK == 2:
				print "\nChoose all genome files you wish to download under the organism name:\n"+NAME+"\n"
				print "\nPlease input the number(s) delimited by comma. For example: 1,3,4,7,8\n" 
				NUMBERS = raw_input()
				NUMBERS = NUMBERS.split(",")
		
				print "Are these the numbers you selected: "+str(NUMBERS)+"\n"
				print "Y/N?"
				GO_ON = raw_input()

		except ValueError:
			GO_ON == "N"

	for i in NUMBERS:
		FOO = chooseme[int(i)-1][1]
		FOO = FOO.split("|")
		DOWNLOAD = FOO[3]

		print "DOWNLOADING: "+DOWNLOAD

		if np.size(NUMBERS) > 1:
			output.write(NAME+"\t"+DOWNLOAD+"\tMULTIPLE_GENOMES\n")

			NAME = NAME.replace(" ", "_")

			os.system(' '.join([
				"esearch -db nuccore -query",
				DOWNLOAD,
				"| efetch -format fasta >>",
				OUTPUT+"/"+NAME+".multiple.genomes.fa"
			]))

			MULTIPLE_GENOME_LIST.append(NAME)

		else:
			output.write(NAME+"\t"+DOWNLOAD+"\tSINGLE\n")

			NAME = NAME.replace(" ", "_")

			os.system(' '.join([
				"esearch -db nuccore -query",
				DOWNLOAD,
				"| efetch -format fasta >>",
				OUTPUT+"/"+NAME+".genome.fa"
			]))
output.close()

###De-replicate all genome files based on ID and Sequence
for NAME in GENOME_NAMES_DICT.iterkeys():

	# LOCATE CORRECT FASTA FILE
	if NAME in MULTIPLE_GENOME_LIST:
		# MAKE SURE FILE EXISTS (SINCE SOME WERE SKIPPED)
		if NAME not in SKIPPED_GENOMES:
			OUT = NAME + ".multiple.genomes.fa"
			PROCEED = "GO"

		else:
			PROCEED = "SKIP"
	else:
		# MAKE SURE FILE EXISTS (SINCE SOME WERE SKIPPED)
		if NAME not in SKIPPED_GENOMES:
			OUT = NAME + ".genome.fa"
			PROCEED = "GO"

		else:
			PROCEED = "SKIP"

	if PROCEED != "SKIP":
		print "\n\n"+OUT

		Seq_LIST = []

	        #Initialize Dict
        	Unique_Dictionary = {}

	        # Open output file
        	genome_output = open(OUTPUT+"/"+re.sub(".fa","",OUT)+".uniq.fa","w")

	        # Create Unique Dictionary entry for each sequence
        	for record in SeqIO.parse(open(OUTPUT+"/"+OUT, "rU"), "fasta") :
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

		if NAME in MULTIPLE_GENOME_LIST:
			os.system("rm "+OUTPUT+"/"+NAME+".multiple.genomes.fa")

		else:
			os.system("rm "+OUTPUT+"/"+NAME+".genome.fa")
			
error.close()
os.system("rm FOO.txt")

stop = timeit.default_timer()
print "This operation took " + str(stop - start) + " seconds.\n"
