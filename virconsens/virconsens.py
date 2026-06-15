#!/usr/bin/env -S python3 -u
import pysam
import argparse
from collections import Counter
import re
import multiprocessing
import itertools
import plotly


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

parser.add_argument('--freqfile',
                    help='Output path for per-position base/indel frequency CSV '
                         '(position,A_freq,C_freq,G_freq,T_freq,depth,del_freq,ins_freq; position is 1-based)',
                    type=str,
                    required = False)

parser.add_argument('-p',
                    '--coverageplot',
                    help='Output path for coverage plot (html file)',
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

parser.add_argument('-k', 
                    '--keepindels',
                    help='Keep 1 and 2 nt indels',
                    action='store_true',
                   required = False)

parser.add_argument('--maxdepth',
                    help='Maximum depth to consider at any position',
                    default=8000,
                    type=int,
                    required = False)


def process_batch(start, stop, bamfile, maxdepth, reference):
    bamfile = pysam.AlignmentFile(bamfile, "rb")
    ref_name = bamfile.references[0]
    refseq = pysam.Fastafile(reference).fetch(ref_name)
    
    pileup = bamfile.pileup(
        contig=ref_name,
        start=start,
        stop=stop,
        ignore_orphans=False,
        min_mapping_quality=0,
        min_base_quality=0,
        truncate=True,
        max_depth=maxdepth
    )
    
    pileup_dict = {}
    for p in pileup:
        pileup_dict[p.reference_pos] = (p.get_query_sequences(add_indels=True), p.get_num_aligned())
    
    variant_rows = []
    freq_rows = []
    for pos in range(start, stop):
        if pos in pileup_dict:
            allele_list, num_aln = pileup_dict[pos]
        else:
            allele_list, num_aln = [], 0
        variant_rows.append(parse_column(pos, allele_list, num_aln, refseq))
        freq_rows.append(parse_column_freq(pos, allele_list, num_aln, refseq))

    return (variant_rows, freq_rows)


def parse_column(ref_pos, allele_list, num_aln, refseq):

    COMB_allele = Counter()
    
    # Initialize the reference allele to 0 in case it is not present in any of the reads
    COMB_allele[(refseq[ref_pos].upper(),refseq[ref_pos].upper())] = 0

    def add_allele(ref, alt):
        if (ref.islower() | alt.islower()):
            COMB_allele[(ref.upper(),alt.upper())] += 1
        else:
            COMB_allele[(ref,alt)] += 1
    
    insert_finder = re.compile(r"(.*)\+\d+(.*)")

    for var in allele_list:
        # Ignore positions that represent deletions
        if '*' in var:
            continue

        # - means next nucleotide is a deletion
        if '-' in var:
            # Determine number of deletions based on the number in the pilup string
            n_del = int(''.join(filter(str.isdigit, var)))
            # Get the nucleotides that were deleted (add 1 to select the reference position plus the deleted nucleotidesy)
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

    # Sort alleles by counts, highest count is major alternative allele
    major = sorted(COMB_allele, key=COMB_allele.get, reverse=True)[0]
    ref_seq = major[0]
    alt_seq = major[1]
    alt_count = COMB_allele[major]
    if num_aln == 0:
        alt_AF = 0.0
    else:
        alt_AF = alt_count / num_aln

    return([str(num_aln), ref_pos, ref_seq, alt_seq, alt_count, alt_AF])


def parse_column_freq(ref_pos, allele_list, num_aln, refseq):
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    del_inside = 0   # '*' counts only
    ins_count = 0

    ref_anchor = refseq[ref_pos].upper()

    for var in allele_list:
        s = var.upper()

        # insertion
        if '+' in s:
            ins_count += 1
            continue

        # deletion
        if '*' in s:
            del_inside += 1
            continue

        # base
        if s and s[0] in base_counts:
            base_counts[s[0]] += 1
        else:
            if ref_anchor in base_counts:
                base_counts[ref_anchor] += 1

    informative = sum(base_counts.values()) + del_inside + ins_count
    if informative == 0:
        return [ref_pos + 1, 0.0, 0.0, 0.0, 0.0, num_aln, 0.0, 0.0]

    A_freq = base_counts['A'] / informative
    C_freq = base_counts['C'] / informative
    G_freq = base_counts['G'] / informative
    T_freq = base_counts['T'] / informative
    del_freq = del_inside / informative
    ins_freq = ins_count / informative

    return [ref_pos + 1, A_freq, C_freq, G_freq, T_freq, num_aln, del_freq, ins_freq]


def main():
    args = parser.parse_args()

    bamfile = pysam.AlignmentFile(args.bam, "rb")
    ref_name = bamfile.references[0]
    refseq = pysam.Fastafile(args.reference).fetch(ref_name)
    
    #Make an array of start-stop intervals to parallelize processing
    genome_length = len(refseq)
    split = int(genome_length/args.cores)
    batch = [[i*split+1,(i+1)*split+1, args.bam, args.maxdepth, args.reference] for i in range(args.cores)]

    #Adjust the last "stop" to be the genome length
    batch[-1][1] = genome_length

    with multiprocessing.Pool(processes=args.cores) as p:
        resultlist = p.starmap(process_batch, iter(batch))

    variant_results = list(itertools.chain.from_iterable(r[0] for r in resultlist))
    freq_results = list(itertools.chain.from_iterable(r[1] for r in resultlist))

    # Build variant_dict for consensus
    variant_dict = {}
    for result in variant_results:
        num_aln, ref_pos, ref_seq, alt_seq, alt_count, alt_AF = result
        variant_dict[ref_pos] = result

    # Variant file
    if args.variantfile:
        with open(args.variantfile, "w") as outfile:
            print("POS", "num_aln", "REF", "ALT", "ALT_count", "ALT_AF", sep='\t', file=outfile)
            for result in variant_results:
                num_aln, ref_pos, ref_seq, alt_seq, alt_count, alt_AF = result
                if int(num_aln) != 0:
                    print(ref_pos + 1, num_aln, ref_seq, alt_seq, alt_count, alt_AF, sep='\t', file=outfile)

    # Frequency file
    if args.freqfile:
        with open(args.freqfile, "w") as ffile:
            print("position", "A_freq", "C_freq", "G_freq", "T_freq", "depth", "del_freq", "ins_freq",
                  sep=',', file=ffile)
            for pos1, A_f, C_f, G_f, T_f, depth, del_f, ins_f in freq_results:
                print(pos1, A_f, C_f, G_f, T_f, depth, del_f, ins_f, sep=',', file=ffile)

    # Coverage plot, when hovering over a position, show the depth and the frequencies of A,C,G,T,del,ins
    if args.coverageplot:
        positions = [row[0] for row in freq_results]
        depths = [row[5] for row in freq_results]
        A_freqs = [row[1] for row in freq_results]
        C_freqs = [row[2] for row in freq_results]
        G_freqs = [row[3] for row in freq_results]
        T_freqs = [row[4] for row in freq_results]
        del_freqs = [row[6] for row in freq_results]
        ins_freqs = [row[7] for row in freq_results]

        fig = plotly.graph_objects.Figure(data=plotly.graph_objects.Scatter(x=positions, y=depths, mode='lines',
            hovertemplate=
                'Position: %{x}<br>'+
                'Depth: %{y}<br>'+
                'A_freq: %{customdata[0]:.2f}<br>'+
                'C_freq: %{customdata[1]:.2f}<br>'+
                'G_freq: %{customdata[2]:.2f}<br>'+
                'T_freq: %{customdata[3]:.2f}<br>'+
                'Del_freq: %{customdata[4]:.2f}<br>'+
                'Ins_freq: %{customdata[5]:.2f}<extra></extra>',
            customdata=list(zip(A_freqs, C_freqs, G_freqs, T_freqs, del_freqs, ins_freqs))
        ))

        fig.update_layout(title='Coverage Plot', xaxis_title='Position', yaxis_title='Depth')
        fig.write_html(args.coverageplot)
    
    # Create consensus sequence
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
            elif abs(len(ref_seq)-len(alt_seq)) in [1,2] and not args.keepindels:
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