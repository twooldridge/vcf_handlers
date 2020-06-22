import re, sys, os, getopt
import operator
import argparse
import gzip
import numpy as np
from collections import Counter

flatten = lambda l: [item for sublist in l for item in sublist]

def process_vcf(vcf_path,out_path):
    file = gzip.open(vcf_path, 'rb')
    outfile = open(out_path,'w')
    #var_dict = dict()
    for l in file:
        l = l.strip().decode('ascii')
        # If 'l' is an info line, skip
        if re.match("#", l):
            outfile.write(l)
        else:
            d = l.strip()
            d = re.split('\t', d)
            chrom, pos, ref, alt = operator.itemgetter(0,1,3,4)(d)
            INFO = d[7]
            genotypes = d[9:]
            genotypes = [re.split(r'[:]+',x)[0] for x in genotypes]
            genotypes = [re.split(r'[|/]+',x) for x in genotypes]
            genotypes = [y for x in genotypes for y in x]
            genotypes = [x for x in genotypes if x not in ['|','/']]
            counts = Counter(genotypes)
            newAN = counts['0'] + counts['1']
            newAC = counts['1']
            newAF = counts['1']/newAN
            # create new info field
            if alt == '.':
                newINFO = re.sub('^AN[^;]*','AN=%s' % newAN,INFO)
            else:
                newINFO = re.sub('^AC[^;]*','AC=%s' % newAC,INFO)
                newINFO = re.sub(';AF[^;]*',';AF=%s' % newAF,newINFO)
                newINFO = re.sub(';AN[^;]*',';AN=%s' % newAN,newINFO)
            # Put vcf line back together
            newRECORD = d
            newRECORD[7] = newINFO
            outfile.write('\t'.join(newRECORD) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", help="gzippeed vcf")
    parser.add_argument("--out", help="output path")

    # Parse arguments
    args = parser.parse_args()

    # Run
    process_vcf(args.vcf,args.out)

if __name__ == "__main__":
    main()
