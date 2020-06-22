import re, sys, os, getopt
import operator
import argparse
import gzip
import numpy as np

def process_haps(hap_file_path):
    hap_dict = {}
    with open(hap_file_path,'r') as data:
        for line in data:
            line = line.strip()
            bits = line.split(' ')
            chrom, pos, ref, alt = operator.itemgetter(0,2,3,4)(bits)
            genotypes = bits[5:]
            if chrom not in hap_dict.keys(): hap_dict[chrom] = {}
            hap_dict[chrom][pos] = [ref,alt,genotypes]
    return(hap_dict)

def switch_phase(genotype,p=0):
    if np.random.binomial(1,p,1) == 1:
        genotype = genotype[::-1]
    return(genotype)

def is_het(genotype):
    het = False
    if (('1' in genotype) and ('0' in genotype)):
        het = True
    return(het)

def check_phase(genotype):
    if '|' in genotype:
        ver = genotype
    else:
        contents = re.split(r'[/]+',genotype)
        if '0' in contents and '1' in contents:
            ver = None
        else:
            ver = genotype
    return(ver)

def process_vcf(vcf_path,error=None):
    file = gzip.open(vcf_path, 'rb')
    var_dict = dict()
    for l in file:
        l = l.strip().decode('ascii')
        # If 'l' is an info line, skip
        if not re.match("#", l):
            d = l.strip()
            d = re.split('\t', d)
            chrom, pos, ref, alt = operator.itemgetter(0,1,3,4)(d)
            if chrom not in var_dict.keys(): var_dict[chrom] = {}
            # If polyallelic, skip to next site now. Only interested in biallelic sites
            if len(alt) > 2:
                continue
                    
            # Handle (skip) indels
            indel = False
            if len(alt) > 1:
                #Added this additional if statement because with only 1 sample GATK cannot distinguish the possibility of another allele. For our purposes we just set that to the reference allele though.
                if allele != "<NON_REF>":
                    indel = True
                    continue
                    
            # Now get genotype info for snps
            if error is not None:
                genotypes = [switch_phase(geno,p=error) if is_het(geno) else geno for geno in d[9:]]
            else:
                genotypes = d[9:]
            genotypes = [re.split(r'[:]+',x)[0] for x in genotypes]
            pre_gen = genotypes
            genotypes = [check_phase(x) for x in genotypes]
            print(pos)
            if None in genotypes:
                continue
            else:
                genotypes = [re.split(r'[|/]+',x) for x in genotypes]
                genotypes = [y for x in genotypes for y in x]
                if '.' in genotypes:
                    continue
                genotypes = [x for x in genotypes if x not in ['|','/']]        
                # Store values
                var_dict[chrom][pos] = [ref,alt,genotypes]
    
    return(var_dict)


def get_state(f_ref,geno_dict):
    ##
    states = dict()
    out_data = []
    seq_length = 0
    f_ref = open(f_ref,'r')
    ##
    for line_count, (l_ref) in enumerate(f_ref):
        l_ref = l_ref.strip()
        if re.match(">", l_ref):
            chrom = re.sub(">|:.*","",l_ref)
            print(chrom)
            seq_length = 0
            states[chrom] = {}
            continue
        if chrom not in geno_dict.keys():
            continue
        if not re.match(">", l_ref):
            line_length = len(l_ref)
            for ix, (b_ref) in enumerate(l_ref):
                # get the current genome position, need to iterate through reference and get each base value. 
                # need to add 1 because python indexes at 0; genomes are at index 1
                pos = seq_length + ix + 1
                pos = str(pos)
                if pos in geno_dict[chrom].keys():
                    if b_ref == geno_dict[chrom][pos][0]:
                        true_ref,true_alt = geno_dict[chrom][pos][0:2]
                    elif b_ref == geno_dict[chrom][pos][1]:
                        true_alt,true_ref = geno_dict[chrom][pos][0:2]
                    else:
                        true_ref,true_alt = ['N','N']
                        continue
                    genotypes = [x for x in geno_dict[chrom][pos][2]]
                    #genotypes = [{'0':true_ref,'1':true_alt}[x] for x in geno_dict[chrom][pos][2]]
                    
                    ## 
                    haps_line = [chrom, '%s_%s' % (chrom,pos),pos,true_ref,true_alt,' '.join(genotypes)]
                    out_data.append(haps_line)
            
            seq_length = seq_length + line_length
    
    return(out_data)




def write_data(out_data,out_prefix):
    out_haps = open('%s.haps' % out_prefix,'w')
    for snp in out_data:
        print(snp)
        out_haps.write(' '.join(snp) + '\n') 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", help="phased vcf file")
    parser.add_argument("--ref", help="reference genome with which to polarize snps")
    parser.add_argument("--out", help="output path. Files will be written to <out>.inp and <out>.thap")
    parser.add_argument("--error", help="Introduce phasing error with frequency <p>. Only compatible with vcf option. Default is no error", default=0)

    # Parse arguments
    args = parser.parse_args()

    # Run
    if args.vcf is not None:
        print("Processing vcf")
        geno_dict = process_vcf(args.vcf, error = float(args.error))
    recoded = get_state(args.ref,geno_dict)
    print("Recoded haplotypes, writing out data")
    write_data(recoded,args.out)


if __name__ == "__main__":
    main()


