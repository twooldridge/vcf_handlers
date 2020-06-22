#!/usr/bin/env python
import re, sys, os, getopt
import gzip
import sys
import subprocess
import os
import argparse
import tempfile as tmp

def process_vcf(vcf_path,sample_map,fraction_missing=None,maf=None):
    file = gzip.open(vcf_path, 'rt')
    pops = set(sample_map.values())
    snp_pop_counts = []
    pop_counts = dict.fromkeys(pops)
    header = list(pop_counts.keys())
    #snp_pop_counts.append('chr\tpos\tref\talt\t' + '\t'.join(header))
    num_alleles = len(sample_map.keys())*2
    for l in file:
        # Get sample info line as list
        if re.match("#CHROM",l):
            s = l.strip('\n')
            s = s.split('\t')
            vcf_sample_order = s[9:]
            pop_order = {}
            for idx,sample in enumerate(vcf_sample_order):
                if sample in sample_map.keys():
                    pop_order[idx] = sample_map[sample]
                else:
                    continue
            if len(pop_order) != len(sample_map):
                sys.exit("Error obtaining pop order from vcf, exiting")
            #pop_order = [sample_map[sample] for sample in sample_order]
            #print(sample_order)
            #print(pop_order)
            print("Sample order obtained...")
            continue
        # If 'l' is an info line, skip
        if not re.match("#", l):
            d = re.split('\t', l)
            pos = d[1]
            ref = [d[3]]
            alt = re.split(",", d[4])
            alleles = ref + alt
            # If polyallelic, skip to next site now. Only interested in biallelic sites
            #if len(alleles) > 2:
                #continue

            # Handle (skip) indels
            #indel = False
            #for allele in alleles:
                # Check length of allele string
                #if len(allele) > 1:
                    #Added this additional if statement because with only 1 sample GATK cannot distinguish the possibility of another allele. For our purposes we just set that to the reference allele though.
                    #if allele != "<NON_REF>":
                    #    indel = True

            # Now get genotype info for snps
            genotypes = d[9:]
            # Subest for the samples we're interested in
            genotypes = [genotypes[i] for i in pop_order.keys()]
            genotypes = [(i.split(':')[0]) for i in genotypes]
            genotypes = [re.split('/|\|',i) for i in genotypes]
            basic_genotypes = flatten(genotypes)
            # Do some filtering
            if fraction_missing is not None:
                num_missing = basic_genotypes.count('.')
                if (num_missing/num_alleles) > float(fraction_missing):
                    continue
            if maf is not None:
                num_minor = basic_genotypes.count('1')
                num_genotyped = num_minor + basic_genotypes.count('0')
            if float(num_genotyped) == 0: continue    
            if (num_minor/float(num_genotyped)) < float(maf):
                    continue
            #genotypes = [[int(float(j)) for j in i] for i in genotypes]
            pop_counts = dict(zip(pops,[[[0],[0] * (len(alleles)-1)] for i in range(len(pops))]))
            pop_list = list(pop_order.values())
            for i,gen in enumerate(genotypes):
                if '.' in gen:
                    continue
                current_pop = pop_list[i]
                #print(current_pop)
                #print(sample_order[i])
                #if len(alleles) > 2:
                num_ref = [0]
                num_alt = [0] * (len(alleles)-1)
                for allele in gen:
                    if allele == '0':
                        num_ref[0] = num_ref[0] + 1
                    else:
                        num_alt[int(allele)-1] = num_alt[int(allele)-1] + 1
                    #if allele == '1':
                    #    num_alt = num_alt + 1
                #print(d[0])
                #print(pos)
                #print(alleles)
                #print(pop_counts[current_pop][0])
                #print(num_ref)
                pop_counts[current_pop][0] = [a+b for a,b in zip(pop_counts[current_pop][0],num_ref)]
                pop_counts[current_pop][1] = [a+b for a,b in zip(pop_counts[current_pop][1],num_alt)]           
            
            if len(snp_pop_counts) == 0:
                snp_pop_counts.append('chr\tpos\tref\talt\t' + '\t'.join(pop_counts.keys()))
            pop_counts = list(pop_counts.values())
            pop_counts = [flatten(l) for l in pop_counts]
            pop_counts = [','.join(map(str, i)) for i in pop_counts]
            pop_counts = '\t'.join(pop_counts)
            snp_pop_counts.append('\t'.join([d[0],pos,ref[0],','.join(alt),pop_counts]))
    
    return(snp_pop_counts)

flatten = lambda l: [item for sublist in l for item in sublist]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gzvcf", help="input gzvcf;unzipped not supported yet, although should be an easy fix")
    parser.add_argument("--map", help="popmap file -- this is a 2-column, tab separated file with sample names in the first column, and the corresponding population in the second")
    parser.add_argument("--fraction_missing",help="fraction of genotypes that can be missing before a site is filtered out. 0 = no missing genotypes allowed")
    parser.add_argument("--maf",help="filter variants by minor allele frequency. Default = 0.025")
    parser.add_argument("--out", help="path to output file")
    

    # Parse arguments
    args = parser.parse_args()
    vcf = args.gzvcf
    map_file = args.map
    missing_cutoff = args.fraction_missing
    maf = args.maf
    if maf is None:
        maf_val = 0.025
    else:
        maf_val = maf
    outfile = args.out

    # Run
    samplemap = get_samplemap(map_file)
    print("Processing vcf...")
    treemix_lines = process_vcf(vcf,samplemap,fraction_missing = missing_cutoff,maf = maf_val)

    # Write
    out = open(outfile, 'w')
    for line in treemix_lines:
        writeable = line
        out.write('%s\n' % (writeable))

def get_samplemap(popmap_path):
    sample_dict = {}
    with open(popmap_path,'r') as data:
        for line in data:
            line = line.strip('\n')
            line = line.split('\t')
            sample = line[0]
            pop = line[1]
            sample_dict[sample] = pop
    data.close()
    return(sample_dict)


if __name__ == "__main__":
 main()
