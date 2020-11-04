#!/usr/bin/env python
import re, sys, os, getopt
import gzip
import sys
import subprocess
import os
import argparse
import tempfile as tmp

def process_vcf(vcf_path, sample_map, allelecounts, fraction_missing=None, maf=None):
    
    ## Grab the vcf
    if '.vcf.gz' in vcf_path:
        opener = gzip.open
    elif '.vcf' in vcf_path:
        opener = open
    else:
        print('file type not recognized')
    file = opener(vcf_path, 'rt')
    
    ## Setup data structures that will record counts, etc.
    pops = set(sample_map.values())
    snp_pop_counts = []
    pop_counts = dict.fromkeys(pops)
    header = list(pop_counts.keys())
    num_alleles = len(sample_map.keys())*2
    
    ## Begin processing file line by line
    last_reported = 0
    for l in file:
        
        ##  If you hit the header line, grab the sample names and store the order in which they appear in the vcf
        if re.match("#CHROM",l):
            s = l.strip('\n')
            s = s.split('\t')
            vcf_sample_order = s[9:]
            pop_order = {}
            ## This is setup so that if a sample is in the vcf but not specified 
            ## in the supplied pop map (not of interest), then it's fine
            for idx,sample in enumerate(vcf_sample_order):
                if sample in sample_map.keys():
                    pop_order[idx] = sample_map[sample]
                else:
                    continue
            ## Sanity check
            if len(pop_order) != len(sample_map):
                sys.exit("Error obtaining pop order from vcf, exiting")
            print("Sample order obtained...")
            continue
        
        ## Skip info lines, start processing when you hit lines without '#'
        if not re.match("#", l):
            d = re.split('\t', l)
            pos = d[1]
            ## Progress report every ~1e6 BP
            if (((int(pos) % 1e6) < 50000) and ((int(pos) - last_reported) >= 1e6)):
                print('Processing %s:%s' %(d[0],pos))
                last_reported = int(pos)
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

            
            ## Now get genotype info for snps
            genotypes = d[9:]
            ## Grabbing the genotypes only for the samples we're interested in
            genotypes = [genotypes[i] for i in pop_order.keys()]
            ## This gets rid of the depth & QC info
            genotypes = [(i.split(':')[0]) for i in genotypes]
            ## This splits up diploid genotypes so that I can count the # of alleles. Will work for haploid as well
            genotypes = [re.split('/|\|',i) for i in genotypes]
            ## Formatting
            basic_genotypes = flatten(genotypes)
            

            ## Now we do some filtering
            ## Missing data
            if fraction_missing is not None:
                num_missing = basic_genotypes.count('.')
                if (num_missing/num_alleles) > float(fraction_missing):
                    continue
            ## Minor allele frequency
            if maf is not None:
                num_minor = basic_genotypes.count('1')
                num_genotyped = num_minor + basic_genotypes.count('0')
            if float(num_genotyped) == 0: continue    
            if (num_minor/float(num_genotyped)) < float(maf):
                    continue
            

            ## This is where the work actually happens
            ## Dictionary to store allele counts, can handle multiallelic sites
            pop_counts = dict(zip(pops,[[[0],[0] * (len(alleles)-1)] for i in range(len(pops))]))
            ## List of pop names that is the same length as the # of samples being processed
            pop_list = list(pop_order.values())
            for i,gen in enumerate(genotypes):
                ## Skip missing
                if '.' in gen:
                    continue
                current_pop = pop_list[i]
                ## Initialize 
                num_ref = [0]
                num_alt = [0] * (len(alleles)-1)
                ## The lines below do the actualy counting
                for allele in gen:
                    if allele == '0':
                        num_ref[0] = num_ref[0] + 1
                    else:
                        num_alt[int(allele)-1] = num_alt[int(allele)-1] + 1
                ## Store results
                pop_counts[current_pop][0] = [a+b for a,b in zip(pop_counts[current_pop][0],num_ref)]
                pop_counts[current_pop][1] = [a+b for a,b in zip(pop_counts[current_pop][1],num_alt)]           
            
            ## BTW, write header, if this hasn't been done already
            if len(snp_pop_counts) == 0:
            	if allelecounts:
                	snp_pop_counts.append('chr\tpos\tref\talt\t' + '\t'.join([i + '_2n' for i in pop_counts.keys()]) + '\t' + '\t'.join([i + '_raw_alts' for i in pop_counts.keys()]) + '\t' + '\t'.join([i + '_freqs' for i in pop_counts.keys()]))
            	else:
            		snp_pop_counts.append('chr\tpos\tref\talt\t' + '\t'.join([i + '_2n' for i in pop_counts.keys()]) + '\t' + '\t'.join([i + '_freqs' for i in pop_counts.keys()]))


            ## Now parse the count data, calculate frequencies, etc. 
            ## The below lines produce a nested list structure like [[1,1], [0,1]] where each nested bracket is a population
            pop_counts = list(pop_counts.values())
            pop_counts = [flatten(l) for l in pop_counts]
            
            ## Frequencies, accounting for how many alleles are actually called in the population
            altfreqs = []
            for popdata in pop_counts:
            	ff = []
            	for altallele in popdata[1:]:
            		ff.append(altallele/sum(popdata))
            	altfreqs.append(ff)
            altfreqs = [','.join(map(str, i)) for i in altfreqs]
            altfreq_string = '\t'.join(altfreqs)
            ## Get 2n value (How many alleles were sequenced?)
            pop_2n = [sum(popdata) for popdata in pop_counts]
            pop_2n = '\t'.join([str(i) for i in pop_2n])

            ## Save data to list that will be written to file
            if allelecounts:
                ## If you also want the raw allele counts, we do the formatting and storing here
                pop_counts = [d[1:] for d in pop_counts]
                pop_counts = [','.join(map(str, i)) for i in pop_counts]
                pop_string = '\t'.join(pop_counts)
            	snp_pop_counts.append('\t'.join([d[0], pos, ref[0],  ','.join(alt), pop_2n,pop_string,altfreq_string]))
            else:
                ## This line is if you just want frequencies
            	snp_pop_counts.append('\t'.join([d[0], pos, ref[0],  ','.join(alt), pop_2n,altfreq_string]))
    
    return(snp_pop_counts)

flatten = lambda l: [item for sublist in l for item in sublist]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", help="input, searches or .vcf or .vcf.gz suffix")
    parser.add_argument("--map", help="popmap file -- this is a 2-column, tab separated file with sample names in the first column, and the corresponding population in the second")
    parser.add_argument("--fraction_missing",help="fraction of genotypes that can be missing before a site is filtered out. 0 = no missing genotypes allowed, 1 = all missing genotypes allowed. Default is None")
    parser.add_argument("--maf",help="filter variants by minor allele frequency. Default = 0.025")
    #parser.add_argument('--freq', help ="output population frequencies in addition to raw allele counts", default=True, action='store_true')
    parser.add_argument('--counts', help ="output raw allele counts in addition to frequences", default=False, action='store_true')
    parser.add_argument("--out", help="path to output file")
    

    # Parse arguments
    args = parser.parse_args()
    vcf = args.vcf
    map_file = args.map
    missing_cutoff = args.fraction_missing
    maf = args.maf
    counts = args.counts
    if maf is None:
        maf_val = 0.025
    else:
        maf_val = maf
    outfile = args.out

    ## Run
    samplemap = get_samplemap(map_file)
    print("Processing vcf...")
    ## This is the bulk of the work
    count_lines = process_vcf(vcf, samplemap, counts, fraction_missing = missing_cutoff, maf = maf_val)

    ## Write file
    out = open(outfile, 'w')
    for line in count_lines:
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
