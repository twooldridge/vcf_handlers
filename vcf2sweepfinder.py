import re, sys, os, sets, getopt
from itertools import izip        
import gzip
import re
import sys
import subprocess
import os
import argparse
import tempfile as tmp

def process_vcf(vcf_path):
    file = gzip.open(vcf_path, 'r')
    var = dict()
    for l in file:
        # If 'l' is an info line, skip
        if not re.match("#", l):
                    d = re.split('\t', l)
                    pos = d[1]
                    ref = [d[3]]
                    alt = re.split(",", d[4])
                    alleles = ref + alt

                    # If polyallelic, skip to next site now. Only interested in biallelic sites
                    if len(alleles) > 2:
                        continue
                    
                    # Handle (skip) indels
                    indel = False
                    for allele in alleles:
                        # Check length of allele string
                        if len(allele) > 1:
                            #Added this additional if statement because with only 1 sample GATK cannot distinguish the possibility of another allele. For our purposes we just set that to the reference allele though.
                            if allele != "<NON_REF>":
                                indel = True
                    
                    # Now get genotype info for snps
                    if not indel:
                            allele_counts = dict()
                            for i in range(len(alleles)):
                                    allele_counts[alleles[i]] = len(re.findall('\s%s[\/|\:]' % i, l)) + len(re.findall('\/%s' % i, l))
                            AN = sum(allele_counts.values())
                    else:
                        continue
                    
                    # Check if there are at least some genotyped samples and alt alleles. 
                    # If the input files have been prepared  & filtered correctly, 
                    # this should always be the case
                    AC = list()
                    AC.append(allele_counts[ref[0]])
                    AC.append(allele_counts[alt[0]])        
                    var[pos] = [alleles,AC,AN]
    return(var)

def get_chromosome(genome, chrom):
    (fd,outfile) = tmp.mkstemp()
    subprocess.call('''samtools faidx {genome} {chrom} > {outfile}'''.format(genome=genome,outfile=outfile,chrom=chrom),shell=True)
    out_f = open(outfile, 'r+')
    locus_name = out_f.next()
    return out_f,outfile

def get_state(f_ref,processed_vcf):
    states = dict()
    seq_length = 0
    for line_count, (l_ref) in enumerate(f_ref):
        if not re.match(">", l_ref):
            l_ref = l_ref.strip()
            line_length = len(l_ref)
            for ix, (b_ref) in enumerate(l_ref):
                # get the current genome position, need to iterate through reference and get each base value.	
                # need to add 1 because python indexes at 0; genomes are at index 1
                pos = seq_length + ix + 1
                pos = str(pos)
                if pos in processed_vcf:
                    states[pos] = b_ref
            seq_length = seq_length + line_length
    return states

    
def polarize_variants(states,var,mono=False):
    polarized_counts = dict()
    for pos in var.keys():
        AN = var[pos][2]
        vcf_alleles = var[pos][0]
        vcf_counts = var[pos][1]
        vcf_ref = vcf_alleles[0]
        vcf_alt = vcf_alleles[1]
        fasta_ref = states[pos]
        fasta_ref = fasta_ref.upper()
	##
        if fasta_ref == 'N':
            ## Ancestral state unknown
            AC = vcf_counts[1]
            fold=1
        else: 
            if fasta_ref == vcf_ref:
                ## No need to repolarize
                AC = vcf_counts[1]
                fold=0
            else:
                ## Repolarization necessary. This statement should work with both normal SNPs and monomorphic sites
                AC = vcf_counts[0]
                fold=0 
        ## The SweepFinder2 manual recommends omitting sites where AC=0, given that the site is polarized. 
        ## If the ancestral state is unknown, AC=0 could actually be a fixed substitution, so such sites are
        ## left in the output file. 
	if ((AC == 0 and fold == 0) and mono == False):
            continue
        else:
	    polarized_counts[pos] = [AC,AN,fold]
    
    return polarized_counts

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chr", help="chromosome for which to run analysis")
    parser.add_argument("--ref", help="reference genome with which to polarize snps")
    parser.add_argument("--vcf", help="input vcf")
    parser.add_argument("--fold", help="whether reference genome should be considered:\n\t'ancestral' (--fold 0)\n\t'reference' (--fold 1)\nIn the later case, the polarization step is skipped entirely")
    parser.add_argument("--keep_monomorphic",action='store_true',default=False,help="output sites where x (alt allele count) is 0")
    parser.add_argument("--out", help="path to output file")
    
    # Parse arguments
    args = parser.parse_args()
    chrom = args.chr
    ref = args.ref
    vcf = args.vcf
    fold = args.fold
    keep_mono = args.keep_monomorphic
    outfile = args.out
    
    # Run
    vcf_dict = process_vcf(vcf)
    if fold == '0':
        f_ref,tempfile = get_chromosome(ref,chrom)
        states = get_state(f_ref,vcf_dict)
        polarized_variants = polarize_variants(states,vcf_dict,mono=keep_mono)
        sweepfinder_format = [[int(k),v[0],v[1],v[2]] for k,v in polarized_variants.iteritems()]
    else:
        sweepfinder_format = [[int(k),v[1][1],v[2],1] for k,v in vcf_dict.iteritems()]
    
    os.remove(tempfile)

    sweepfinder_format = sorted(sweepfinder_format, key = lambda x: int(x[0]))
    out = open(outfile, 'w')
    out.write('position\tx\tn\tfolded\n')
    for line in sweepfinder_format:
        pos,x,n,fold = line[:4]
        out.write('%s\t%s\t%s\t%s\n' % (pos,x,n,fold))
    
if __name__ == "__main__":
    main()

