#!/usr/bin/env python -u
import re, sys, os, getopt
import subprocess
import argparse
import tempfile as tmp
import allel
from operator import itemgetter
print('Libraries loaded...')



def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--vcf", help="input vcf")
	parser.add_argument("--ref", help="reference fasta")
	parser.add_argument("--gff", help="gff3 file corresponding to reference fasta")
	parser.add_argument("--gene", help="gene name to query and return sequences for")
	parser.add_argument("--bed", help= "BED file of regions to extract and translate. Strand must be indicated in the 4th column!!! If specified, takes priority over --gff and --gene"
	parser.add_argument("--positions", help="Comma separated list of variant positions to swap in while processing reference fasta. REQUIRED")
	parser.add_argument('--prefix', help="Prefix for output file")

	print("Parsing arguments...")

	# Parse arguments
	args = parser.parse_args()
	vcf = args.vcf
	if '.gz' in vcf:
		comp = True
	ref = args.ref
	gff = args.gff
	gene = args.gene
	positions = args.positions

	## Process region info ##
	if gene:
		print("Using input gff file to grab sequences")
	else:
		print("Using gff and gene ID '%s' to grab sequences" % gene)
		tmpgff = tempfile.NamedTemporaryFile(delete=False)
		bash_cmd = ('''grep "\<{gene}\>" {gff} > {ougff}\n''').format(gene=gene,gff=gff,outgff=tmpgff.name)
		gff = tmpgff.name
	with open(gff,'r') as data:
		for idx, line in enumerate(data):
			line = line.strip('\n')
			line = line.split('\t')
			if idx ==1:
				chrom, start_pos, end_pos = [line[i] for i in [0,3,4]]
			else:
				chrom, blank, end_pos = [line[i] for i in [0,3,4]]
	data.close()

	## Process fasta ##
	tmpfasta = tempfile.NamedTemporaryFile(delete=False)
	gffread_cmd = ('''module load gffread/2.2.1-fasrc01\n'''
					''' gffread -w {outfasta} -g {ref} {gff}''').format(ref=ref,gff=gff,outfasta=tmpfasta.name)

	## Process vcf ##
	variants = allel.read_vcf(vcf,region='%s:%s-%s' % (chrom,start_pos,end_pos))
	vcf_positions = list(variants['variants/POS'])
	alts = list(variants['variants/ALT'])
	replacements = dict(zip(vcf_positions, alts))

	## Process positions string
	focal_vars = positions.split(',')

	## Create alternate fastas ##
	with open(tmpfasta.name,'r') as data:
		pos = start_pos
		old_sequence = ''
		new_sequence = ''
		for line_count, l_ref in enumerate(data):
			if '>' in l_ref: 
				continue
			else:
				for ix, b_ref in enumerate(l_ref.strip()):
					pos = pos + ix
					if (pos in replacements.keys()) and (pos in focal_vars):
						# Currently setup to only grab first alternate allele
						new_ref = replacements[pos][0]
					else:
						new_ref = b_ref
					old_sequence = old_sequence + b_ref
					new_sequence = new_sequence + b_ref
		sequences = {'old':Seq(old_sequence),'new':Seq(new_sequence)}
	data.close()

	## Translate to amino acid sequence ##
	# First, reverse complement if necessary # 
	if strand == '-':
		for ID, seq in sequences.items():
			sequences[ID] = sequences[ID].reverse_complement()
	aminos = {}
	for ID, seq in sequences.items():
		aminos[ID] = seq.translate()

	## Write out nucleotide and amino acid fastas
	for ID, seq in sequences.items():


	

if __name__ == "__main__":
	main()