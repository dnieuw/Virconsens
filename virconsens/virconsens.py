#!/usr/bin/env python3
import pysam
import argparse
from collections import Counter
import re
import multiprocessing
import itertools

parser = argparse.ArgumentParser(description='Virconsens')

parser.add_argument('-b',
                    '--bam',
                    help="BAM file from which to create a consensus",
                    type=str,
                    required = True)

parser.add_argument('-o', 
                    '--out',
                    help='Output path for consensus fasta',
                    type=str,
                   required = True)

parser.add_argument('-n', 
                    '--outname',
                    help='Name to be given to the output consensus sequence',
                    type=str,
                   required = True)

parser.add_argument('-r', 
                    '--reference',
                    help='Reference genome fasta file',
                    type=str,
                   required = True)

parser.add_argument('-vf', 
                    '--variantfile',
                    help='Output path for variant tsv file',
                    type=str,
                   required = False)

parser.add_argument('-c', 
                    '--cores',
                    help='Number of cores to use for processing',
                    default=1,
                    type=int,
                   required = False)

parser.add_argument('-d', 
                    '--mindepth',
                    help='Minimal depth at which to not consider any alternative alleles',
                    default=30,
                    type=int,
                   required = False)

parser.add_argument('-af', 
                    '--minAF',
                    help='Minimal allele frequency to output',
                    default=0.1,
                    type=float,
                   required = False)

def process_batch(start, stop, bamfile, reference):
    bamfile = pysam.AlignmentFile(bamfile, "rb")
    ref_name = bamfile.references[0]
    refseq = pysam.Fastafile(reference).fetch(ref_name)
    
    pileup = bamfile.pileup(contig=ref_name, start=start, stop=stop, ignore_orphans=False, min_mapping_quality=0, min_base_quality=0, truncate=True)
    
    return([parse_column(p.reference_pos, p.get_query_sequences(add_indels=True),p.get_num_aligned(), refseq) for p in pileup])

def parse_column(ref_pos, allele_list, num_aln, refseq):

    COMB_allele = Counter()
    
    #Initialize the reference allele to 0 in case it is not present in any of the reads
    COMB_allele[(refseq[ref_pos].upper(),refseq[ref_pos].upper())] = 0
    
    def add_allele(ref, alt):
        if (ref.islower() | alt.islower()):
            COMB_allele[(ref.upper(),alt.upper())] += 1
        else:
            COMB_allele[(ref,alt)] += 1
    
    insert_finder = re.compile("(.*)\+\d+(.*)")

    for var in allele_list:
        #Ignore positions that represent deletions
        if '*' in var:
            continue

        # - means next nucleotide is a deletion
        if '-' in var:
            #Determine number of deletions based on the number in the pilup string
            n_del = int(''.join(filter(str.isdigit, var)))
            #Get the nucleotides that were deleted (add 1 to select the reference position plus the deleted nucleotidesy)
            if var.islower():
                var=refseq[ref_pos:(ref_pos+n_del+1)].lower()
            else:
                var=refseq[ref_pos:(ref_pos+n_del+1)]
            add_allele(var, refseq[ref_pos])
        # + means next nucleotide is a insertion
        elif '+' in var:
            var = ''.join(insert_finder.match(var).groups())
            add_allele(refseq[ref_pos], var)
        else:
            add_allele(refseq[ref_pos], var)

    #Sort alleles by counts, highest count is major alternative allele
    major = sorted(COMB_allele, key=COMB_allele.get, reverse=True)[0]
    ref_seq = major[0]
    alt_seq = major[1]
    alt_count = COMB_allele[major]
    alt_AF = alt_count/num_aln

    return([str(num_aln), ref_pos, ref_seq, alt_seq, alt_count, alt_AF])

def main():
    args = parser.parse_args()

    bamfile = pysam.AlignmentFile(args.bam, "rb")
    ref_name = bamfile.references[0]
    refseq = pysam.Fastafile(args.reference).fetch(ref_name)
    
    #Make an array of start-stop intervals to parallelize processing
    genome_length = len(refseq)
    split = int(genome_length/args.cores)
    batch = [[i*split+1,(i+1)*split+1, args.bam, args.reference] for i in range(args.cores)]

    #Adjust the last "stop" to be the genome length
    batch[-1][1] = genome_length

    with multiprocessing.Pool(processes=args.cores) as p:
        resultlist = p.starmap(process_batch, iter(batch))

    variant_dict = {}
    if args.variantfile:
        outfile = open(args.variantfile,"w")
        print("ref_pos", "num_aln", "REF", "ALT", "ALT_count", "ALT_AF", sep='\t', file=outfile)
        
    for result in itertools.chain.from_iterable(resultlist):
        num_aln, ref_pos, ref_seq, alt_seq, alt_count, alt_AF = result
        variant_dict[ref_pos] = result
        if args.variantfile:
            print(ref_pos, num_aln, ref_seq, alt_seq, alt_count, alt_AF, sep='\t', file=outfile)
    
    if args.variantfile:
        outfile.close()

    consensus = []
    pos = 0
    while pos < genome_length:
        if pos in variant_dict:
            num_aln, ref_pos, ref_seq, alt_seq, alt_count, alt_AF = variant_dict[pos]

            if (alt_AF < args.minAF or int(num_aln) < args.mindepth):
                consensus.append("N")
                pos += 1
                continue
            #Ignore indels of 1 or 2 nt
            elif abs(len(ref_seq)-len(alt_seq)) in [1,2]:
                consensus.append(ref_seq)
            else:
                consensus.append(alt_seq)

            #If ref is bigger than alt, we have a deletion and have to increment position by the deletion size to skip the following positions
            if len(ref_seq) > len(alt_seq):
                pos += (len(ref_seq)-len(alt_seq))
        else:
            consensus.append("N")
        pos += 1

    consensus = ''.join(consensus)

    with open(args.out, 'w') as out_consensus:
        print(''.join(['>',args.outname]), file=out_consensus)
        print(consensus, file=out_consensus)

if __name__ == '__main__':
    main()