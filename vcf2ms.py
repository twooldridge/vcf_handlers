#!/usr/bin/env/python

import re, sys, os, getopt
import gzip
import sys
import subprocess
import os
import argparse
import tempfile as tmp

def process_vcf(vcf_path):
    file = gzip.open(vcf_path, 'rt')
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
                    genotypes = d[9:]
                    genotypes = [re.sub(r'/|\|','',(i.split(':')[0])) for i in genotypes]
                    genotypes = "".join(genotypes)

            else:
                continue
      
            var[pos] = ["".join(ref),"".join(alt),genotypes]
    
    return(var)

def get_chromosome(genome, chrom):
    (fd,filename) = tmp.mkstemp()
    subprocess.call('samtools faidx %s %s > %s' % (genome, chrom, filename), shell=True)
    out_f = open(filename, 'r+')
    locus_name = next(out_f)
    return out_f,filename,locus_name

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

def polarize_variants(states,var):
    polarized_variants = {}
    substitutions = {"0": "1", "1": "0"}
    for pos,items in var.items():
        # Get all the info
        vcf_ref = items[0]
        vcf_alt = items[1]
        genotypes = items[2]
        fasta_ref = states[pos]
        fasta_ref = fasta_ref.upper()
        # Do we need to polarize
        if vcf_ref == fasta_ref:
            polarized_genotypes = genotypes
        elif vcf_alt == fasta_ref:
            # Alternate allele is same as ancestral, re-polarize
            polarized_genotypes = replace(genotypes,substitutions)
        else:
            # Ancestral allele doesn't match ref. or alt in vcf. Should only happen when ancestral is 'N'
            polarized_genotypes = None
        polarized_variants[pos] = polarized_genotypes
    return polarized_variants

def replace(string, substitutions):
    substrings = sorted(substitutions, key=len, reverse=True)
    regex = re.compile('|'.join(map(re.escape, substrings)))
    return regex.sub(lambda match: substitutions[match.group(0)], string)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chr", help="chromosome for which to run analysis")
    parser.add_argument("--ref", help="reference genome with which to polarize snps")
    parser.add_argument("--vcf", help="input vcf")
    parser.add_argument("--out", help="path to output file")
    
    # Parse arguments
    args = parser.parse_args()
    chrom = args.chr
    ref = args.ref
    vcf = args.vcf
    outfile = args.out
    
    # Run
    vcf_dict = process_vcf(vcf)
    f_ref,tempfile,fasta_header = get_chromosome(ref,chrom)
    states = get_state(f_ref,vcf_dict)
    polarized_variants = polarize_variants(states,vcf_dict)
    os.remove(tempfile)
    
    # Write
    out = open(outfile, 'w')
    out.write(r'\\')
    out.write('\n')
    out.write(fasta_header)
    for pos,genotypes in polarized_variants.items():
        if genotypes is None:
            continue
        out.write('%s\n' % (genotypes))


if __name__ == "__main__":
 main()
