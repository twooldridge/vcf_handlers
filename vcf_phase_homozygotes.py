#!/usr/bin/env python
import re, sys, os, getopt
import operator
import argparse
import gzip
import numpy as np
from collections import Counter

flatten = lambda l: [item for sublist in l for item in sublist]

def process_vcf(vcf_path, out_path, unphased_het_to_missing = False):
    if '.vcf.gz' in vcf_path:
        gzipped = True
        file = gzip.open(vcf_path, 'rb')
    else:
        gzipped = False
        file = open(vcf_path, 'r')
    outfile = open(out_path,'w')
    #var_dict = dict()
    for l in file:
        if gzipped:
            l = l.strip().decode('ascii')
        else:
            l = l.strip()
        # If 'l' is an info line, skip
        if re.match("#", l):
            outfile.write(l + '\n')
        else:
            d = l.strip()
            d = re.split('\t', d)
            SNPinfo = d[:9]
            genotypes = d[9:]
            #genotypes = [re.split(r'[:]+',x)[0] for x in genotypes]
            corrected_genotypes = []
            for genotype_string in genotypes:
                comps = re.split(r'[:]+',genotype_string)
                GEN = comps[0]
                STATS = comps[1:]
                if '|' in GEN:
                    corrected_genotype_string = [GEN] + STATS
                else:
                    if ((GEN == '0/1') or (GEN == '1/0')):
                        if (unphased_het_to_missing):
                            GEN = './.'
                            corrected_genotype_string = [GEN] + STATS
                        else:
                            print(l)
                            sys.exit('Found unphased heterozygote, this genotype cannot be corrected. If you would like to set this to a missing (.|.) genotype instead, set --unphased_het_to_missing True; exiting')
                    else:
                        if GEN == './.':
                            corrected_genotype_string = [GEN] + STATS
                        else:
                            GEN = re.sub('/','|',GEN)
                            corrected_genotype_string = [GEN] + STATS
                ##
                corrected_genotype_string = ':'.join(corrected_genotype_string)
                corrected_genotypes.append(corrected_genotype_string)
            # Now recombine fields and write
            newRECORD = SNPinfo + corrected_genotypes
            newRECORD = '\t'.join(newRECORD)
            outfile.write(newRECORD + '\n')
    file.close()
    outfile.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", help="uncompressed (.vcf) or compressed (.vcf.gz) vcf", metavar='STRING')
    parser.add_argument("--unphased_het_to_missing", help="set unphased heterozygous genotypes to '.|.' Default = False ", default = False, metavar='True/False')
    parser.add_argument("--output", help="output file path (uncompressed vcf)", metavar='STRING')

    # Parse arguments
    args = parser.parse_args()

    # Run
    process_vcf(vcf_path = args.vcf, out_path = args.output, unphased_het_to_missing = args.unphased_het_to_missing)

if __name__ == "__main__":
    main()
                       