import re, sys, os, getopt
import gzip
import sys
import subprocess
import os
import argparse
import tempfile as tmp

def process_vcf(vcf_path,sample_map,fraction_missing=None):
    file = gzip.open(vcf_path, 'rt')
    pops = set(sample_map.values())
    snp_pop_counts = {}
    pop_counts = dict.fromkeys(pops)
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
                    genotypes = [(i.split(':')[0]) for i in genotypes]
                    genotypes = [re.split('/|\|',i) for i in genotypes]
                    if fraction_missing is not None:
                        basic_genotypes = flatten(genotypes)
                        num_missing = basic_genotypes.count('.')
                        if (num_missing/num_alleles) > fraction_missing:
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
                    
                    
                    ref_counts = [v[0] for v in pop_counts.values()]
                    alt_counts = [v[1] for v in pop_counts.values()]
                    #pop_counts = list(pop_counts.values())
                    #pop_counts = [','.join(map(str, i)) for i in pop_counts]
                    #pop_counts = ' '.join(pop_counts)
                    snp_pop_counts[pos] = [ref,alt,ref_counts,alt_counts]
            else:
                continue
    
    return(snp_pop_counts)


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
    for pos,items in var.items():
        # Get all the info
        vcf_ref = items[0][0]
        vcf_alt = items[1][0]
        ref_counts = items[2]
        alt_counts = items[3]
        fasta_ref = states[pos]
        fasta_ref = fasta_ref.upper()
        # Do we need to polarize
        if vcf_ref == fasta_ref:
            new_ref_counts = ref_counts
            new_alt_counts = alt_counts
        elif vcf_alt == fasta_ref:
            # Alternate allele is same as ancestral, re-polarize
            new_ref_counts = alt_counts
            new_alt_counts = ref_counts
        else:
            # Ancestral allele doesn't match ref. or alt in vcf. Should only happen when ancestral is 'N'
            new_ref_counts = None
            new_alt_counts = None
        polarized_variants[pos] = [vcf_ref,fasta_ref,new_ref_counts,new_alt_counts]
    return polarized_variants

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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", help="input vcf")
    parser.add_argument("--map", help="popmap file -- this is a 2-column, tab separated file with sample names in the first column, and the corresponding population in the second")
    parser.add_argument("--anc", help="fasta file for ancestral reference sequence")
    parser.add_argument("--chrom", help="chromosome to operate on")
    parser.add_argument("--fraction_missing",help="fraction of genotypes that can be missing before a site is filtered out. 0 = no missing genotypes allowed")
    parser.add_argument("--out", help="path to output file")
    

    # Parse arguments
    args = parser.parse_args()
    vcf = args.vcf
    map_file = args.map
    anc_fasta = args.anc
    chrom = args.chrom
    fraction_missing = args.fraction_missing
    outfile = args.out

    # Run
    samplemap = get_samplemap(map_file)
    vcf_dict = process_vcf(vcf,samplemap,fraction_missing = fraction_missing)
    f_ref,tempfile,fasta_header = get_chromosome(anc_fasta,chrom)
    states = get_state(f_ref,vcf_dict)
    polarized_variants = polarize_variants(states,vcf_dict)
    os.remove(tempfile)


    # Write
    out = open(outfile, 'w')
    pops = ' '.join(set(samplemap.values()))
    out.write('Ref Anc Allele1 %s Allele2 %s Position\n' % (pops,pops))
    for pos,vt in polarized_variants.items():
        formatted_line = ('''-{ref_allele}- -{anc_allele}- {ref_allele} {ref_counts} {anc_allele} {alt_counts} {pos}''').format(ref_allele=vt[0],anc_allele=vt[1],ref_counts=' '.join(str(v) for v in vt[2]),alt_counts=' '.join(str(v) for v in vt[3]),pos=pos)
        out.write('%s\n' % (formatted_line))


if __name__ == "__main__":
 main()
