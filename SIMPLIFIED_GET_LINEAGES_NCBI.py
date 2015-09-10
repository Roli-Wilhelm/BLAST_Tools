#!/usr/bin/python
import sys, os, re, getopt, glob, subprocess, os.path, numpy as np, time
import timeit
from cogent.parse.ncbi_taxonomy import NcbiTaxonomyFromFiles

start = timeit.default_timer()

Usage = """
REQUIRED ARGUMENTS:
	-d	a file containing two columns:

		1) Taxonomy Name
		2) Taxonomy ID (NCBI)

	-r	the name of rank desired ("superkingdom","phylum","class","order","family","genus")

DEPENDENCIES:
		'names.dmp' & 'nodes.dmp'	Downloadable from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip

Usage:  ./GET_LINEAGES_NCBI.py -d FILE_NAMES.tsv

"""

if len(sys.argv)<2:
        print Usage
        sys.exit()

# Initialize all containers for input
FILE_LIST=''

# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"d:r:")

###############################
# o == option
# a == argument passed to the o
###############################
for o, a in myopts:
    if o == '-d':
	FILE_LIST= a
    if o == '-r':
	ranks= [a]

## Define function for pulling lineage info from NCBI nodes and names files
def get_lineage(node, my_ranks):
	ranks_lookup = dict([(r,idx) for idx,r in enumerate(my_ranks)])
	lineage = [None] * len(my_ranks)
	curr = node
	while curr.Parent is not None:
		if curr.Rank in ranks_lookup:
			lineage[ranks_lookup[curr.Rank]] = curr.Name
		curr = curr.Parent
	return lineage

## IMPORT DIR NAMES into DICTIONARY

tree = NcbiTaxonomyFromFiles(open('/home/roli/db/nodes.dmp'), open('/home/roli/db/names.dmp'))
root = tree.Root

output = open("LINEAGES.tsv", "w")

with open(FILE_LIST) as f:
        for line in f:
		line = line.strip("\r\n")
		line = line.split("\t")

		NAME = line[0]
		TaxID = int(line[1])

		node = tree.ById[TaxID]
		tax = get_lineage(node, ranks)
		tax = str(tax[0]).lower()

		tax = tax[:1].upper() + tax[1:]

		output.write(NAME+"\t"+tax+"\n")

stop = timeit.default_timer()

print "This operation took " + str(stop - start) + " seconds."


