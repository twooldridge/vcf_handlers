#!/usr/bin/env python3

import sys,os,re,glob,shutil,pickle,subprocess,time,itertools,random
import pandas as pd
import numpy as np
import argparse
import gzip

def process_vcf(vcf_path,sample_map,max_missing=None):
    #file = open(unzipped_vcf_path, 'rt')
    file = gzip.open(vcf_path, 'rt')
    pops = set(sample_map.values())
    snp_pop_counts = {}
    pop_counts = dict.fromkeys(pops)
    header = list(pop_counts.keys())
    num_alleles = len(sample_map.keys())*2
    for l in file:
        # Get sample info line as list
        if re.match("#CHROM",l):
            s = l.strip('\n')
            s = s.split('\t')
            sample_order = s[9:]
            pop_order = [sample_map[sample] for sample in sample_order]
            continue
            
        # If 'l' is an info line, skip
        if not re.match("#", l):
            d = re.split('\t', l)
            chrom = d[0]
            pos = d[1]
            ref = [d[3]]
            alt = re.split(",", d[4])
            alleles = ref + alt

            # If polyallelic, skip to next site now. Only interested in biallelic sites
            if len(alleles) > 2:
                continue

            # Handle (skip) indels
            indel = False
            #for allele in alleles:
            #    # Check length of allele string
            #    if len(allele) > 1:
            #        #Added this additional if statement because with only 1 sample GATK cannot distinguish the possibility of another allele. For our purposes we just set that to the reference allele though.
            #        if allele != "<NON_REF>":
            #            indel = True
#
            # Now get genotype info for snps
            if not indel:
                    genotypes = d[9:]
                    genotypes = [(i.split(':')[0]) for i in genotypes]
                    genotypes = [re.split('/|\|',i) for i in genotypes]
                    if max_missing is not None:
                        basic_genotypes = flatten(genotypes)
                        num_missing = basic_genotypes.count('.')
                        if (num_missing/num_alleles) > max_missing:
                            continue
                    #genotypes = [[int(float(j)) for j in i] for i in genotypes]
                    pop_counts = dict(zip(pops,[[0,0] for i in range(len(pops))]))
                    for i,gen in enumerate(genotypes):
                        if '.' in gen:
                            continue
                        current_pop = pop_order[i]
                        num_ref = 0
                        num_alt = 0
                        for allele in gen:
                            if allele == '0':
                                num_ref = num_ref + 1
                            if allele == '1':
                                num_alt = num_alt + 1
                        pop_counts[current_pop][0] = pop_counts[current_pop][0] + num_ref
                        pop_counts[current_pop][1] = pop_counts[current_pop][1] + num_alt           
                    
                    snp_pop_counts[pos] = pop_counts
                    #pop_counts = list(pop_counts.values())
                    #pop_counts = [','.join(map(str, i)) for i in pop_counts]
                    #pop_counts = ' '.join(pop_counts)
                    #snp_pop_counts.append(pop_counts)
            else:
                continue
    
    return(snp_pop_counts)

flatten = lambda l: [item for sublist in l for item in sublist]

def get_freq_from_count(snp_counts_dict):
    freq_dict = {}
    for snp,counts_dict in snp_counts_dict.items():
        d = {}
        for pop,counts in counts_dict.items():
            if(counts[1]+ counts[0]) == 0:
                continue
            freq = counts[1]/(counts[0] + counts[1])
            d[pop] = freq
        if len(d.keys()) < len(counts_dict.keys()):
            continue
        freq_dict[snp] = d
    return(freq_dict)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", help="input vcf")
    parser.add_argument("--map", help="popmap file -- this is a 2-column, tab separated file with sample names in the first column, and the corresponding population in the second")
    parser.add_argument("--fraction_missing",help="fraction of genotypes that can be missing before a site is filtered out. 0 = no missing genotypes allowed")
    parser.add_argument("--maf",help="filter variants by minor allele frequency. Default = 0.025")
    parser.add_argument("--out", help="path to output file")
    parser.add_argument("--pos", help="output single-column file of relative site positions", default="None")


    # Parse arguments
    args = parser.parse_args()
    vcf = args.vcf
    map_file = args.map
    missing_cutoff = args.fraction_missing
    maf = args.maf
    if maf is None:
        maf_val = 0.025
    else:
        maf_val = maf
    outfile = args.out
  
    # Run
    sample_map = {}
    with open(map_file,'r') as data:
        for line in data:
            line = line.strip()
            sample,pop = line.split('\t')
            sample_map[sample] = pop
    data.close()

    a = process_vcf(vcf,sample_map,max_missing=0.25)
    b = get_freq_from_count(a)
    c = {}
    for pos,freqs in b.items():
        if max(freqs.values()) < float(maf_val):
            continue
        elif min(freqs.values()) == 1.0: #Reove fixed monomorphic
            continue
        else:
            c[pos] = freqs
    ## 
    if args.pos is not None:
        positions = list(c.keys())
        start = float(positions[0])
        end = float(positions[-1])
        size = end - start 
        #rel_positions = [(float(x) - start)/size for x in positions]
        with open('%s.positions' % outfile,'w') as outpos:
            for i in positions:
                outpos.write(str(i) + '\n')    

    ## Format for writing out a POP(rows) X SNP(cols) dataframe
    lst = []
    pop_order = ['SW','BK','BW','NUB']
    for pos,freqs in c.items():
        if max(freqs.values()) < float(maf_val):
            continue
        neworder = [freqs[pop] for pop in pop_order]
        lst.append(neworder)

    df = pd.DataFrame(lst)
    df = df.transpose()

    ## Write out data
    df.to_csv(outfile,sep="\t",index=False,header=False)


if __name__ == "__main__":
    main()
