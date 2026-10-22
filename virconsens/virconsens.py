#!/usr/bin/env -S python3 -u
import pysam
import argparse
from collections import Counter
import re
import multiprocessing
import itertools
import plotly
import plotly.graph_objects as go


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
                    help='Minimum allele frequency required for inclusion in the consensus',
                    default=0.1,
                    type=float,
                   required = False)

parser.add_argument('-k', 
                    '--keepindels',
                    help='Keep 1 and 2 nt indels',
                    action='store_true',
                   required = False)

parser.add_argument('-a',
                    '--ambiguous',
                    help='Use IUPAC ambiguity codes',
                    action='store_true')

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


IUPAC_CODE_MAP = {
    frozenset({'A'}): 'A',
    frozenset({'C'}): 'C',
    frozenset({'G'}): 'G',
    frozenset({'T'}): 'T',
    frozenset({'A', 'G'}): 'R',
    frozenset({'C', 'T'}): 'Y',
    frozenset({'G', 'C'}): 'S',
    frozenset({'A', 'T'}): 'W',
    frozenset({'G', 'T'}): 'K',
    frozenset({'A', 'C'}): 'M',
    frozenset({'A', 'C', 'G'}): 'V',
    frozenset({'A', 'C', 'T'}): 'H',
    frozenset({'A', 'G', 'T'}): 'D',
    frozenset({'C', 'G', 'T'}): 'B',
    frozenset({'A', 'C', 'G', 'T'}): 'N',
}


def bases_to_iupac(bases: set) -> str:
    return IUPAC_CODE_MAP.get(frozenset(bases), 'N')


def choose_iupac_base(freq_row, threshold: float):
    """Select an IUPAC code based on base frequencies above a threshold.
    
    Determines which bases have frequencies >= the specified threshold
    and returns the corresponding IUPAC ambiguity code.
    
    Args:
        freq_row (list): Frequency row [pos, A_freq, C_freq, G_freq, T_freq, depth, del_freq, ins_freq].
        threshold (float): Frequency threshold (0.0-1.0) for including a base.
    
    Returns:
        str or None: IUPAC code if bases meet threshold, None if no bases qualify.
    """
    base_freqs = {
        'A': freq_row[1],
        'C': freq_row[2],
        'G': freq_row[3],
        'T': freq_row[4],
    }
    bases = {base for base, freq in base_freqs.items() if freq >= threshold}
    if not bases:
        return None
    if len(bases) == 1:
        return next(iter(bases))
    return bases_to_iupac(bases)


def allele_is_indel(ref_seq, alt_seq):
    """Determine if an allele represents an insertion or deletion.
    
    Args:
        ref_seq (str): Reference sequence allele.
        alt_seq (str): Alternate sequence allele.
    
    Returns:
        bool: True if the alleles have different lengths (insertion or deletion).
    """
    return len(ref_seq) != len(alt_seq)


def select_consensus_allele(variant_row, freq_row, ambiguous, minAF, mindepth, keepindels):
    """Select the consensus base for a position based on multiple criteria.
    
    Decision logic (in order):
    1. Returns 'N' if depth < mindepth
    2. If ambiguous is True and allele is SNV: returns IUPAC code based on minAF threshold
    3. Returns 'N' if alt_AF < minAF
    4. For small indels (1-2 nt) when not keeping indels: returns reference to avoid frameshift
    5. Otherwise returns the alternate allele
    
    Args:
        variant_row (list): Variant info [num_aln, pos, ref_seq, alt_seq, alt_count, alt_AF].
        freq_row (list): Frequency info [pos, A_freq, C_freq, G_freq, T_freq, depth, del_freq, ins_freq].
        ambiguous (bool): Whether to use IUPAC ambiguity codes.
        minAF (float): Minimum allele frequency threshold.
        mindepth (int): Minimum depth threshold.
        keepindels (bool): Whether to keep small indels in consensus.
    
    Returns:
        str: Consensus base/IUPAC code for the position.
    """
    num_aln = int(variant_row[0])
    ref_seq = variant_row[2]
    alt_seq = variant_row[3]
    alt_AF = variant_row[5]

    if num_aln < mindepth:
        return 'N'

    if ambiguous is True and not allele_is_indel(ref_seq, alt_seq):
        iupac = choose_iupac_base(freq_row, minAF)
        return iupac if iupac is not None else 'N'

    if alt_AF < minAF:
        return 'N'

    if allele_is_indel(ref_seq, alt_seq) and abs(len(ref_seq) - len(alt_seq)) in [1, 2] and not keepindels:
        return ref_seq

    return alt_seq


def parse_column_freq(ref_pos, allele_list, num_aln, refseq):
    """Calculate per-position base and indel frequencies from aligned reads.
    
    Counts occurrences of A, C, G, T, deletions ('*'), and insertions ('+')
    at a genomic position, then normalizes to frequencies. All counts are
    relative to informative bases (excluding completely deleted reads).
    
    Args:
        ref_pos (int): 0-based reference position.
        allele_list (list): List of allele strings from pileup.
        num_aln (int): Total number of aligned reads at this position.
        refseq (str): Reference sequence.
    
    Returns:
        list: [position_1based, A_freq, C_freq, G_freq, T_freq, depth, del_freq, ins_freq]
              Note: position is converted to 1-based for output.
    """
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    del_inside = 0   # Count of '*' (bases deleted in reads)
    ins_count = 0    # Count of '+' (insertions in reads)

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


def create_coverage_plot(freq_results, output_path, minAF, mindepth):
    # Coverage plot, when hovering over a position, show the depth and
    # the frequencies of A,C,G,T,del,ins.
    # Positions containing insertions/deletions above minAF are marked.

    positions = [row[0] for row in freq_results]
    depths = [row[5] for row in freq_results]
    A_freqs = [row[1] for row in freq_results]
    C_freqs = [row[2] for row in freq_results]
    G_freqs = [row[3] for row in freq_results]
    T_freqs = [row[4] for row in freq_results]
    del_freqs = [row[6] for row in freq_results]
    ins_freqs = [row[7] for row in freq_results]

    fig = go.Figure()

    # Main coverage trace
    fig.add_trace(
        go.Scatter(x=positions, y=depths, mode='lines', name='Coverage', line=dict(color='royalblue', width=2),
            customdata=list(zip(A_freqs, C_freqs, G_freqs, T_freqs, del_freqs, ins_freqs)),
            hovertemplate=
                'Position: %{x}<br>' +
                'Depth: %{y}<br>' +
                'A_freq: %{customdata[0]:.2f}<br>'+
                'C_freq: %{customdata[1]:.2f}<br>'+
                'G_freq: %{customdata[2]:.2f}<br>'+
                'T_freq: %{customdata[3]:.2f}<br>'+
                'Del_freq: %{customdata[4]:.2f}<br>'+
                'Ins_freq: %{customdata[5]:.2f}<extra></extra>'
        )
    )

    # Collect insertion positions
    ins_positions = []
    ins_depths = []
    ins_hover = []

    # Collect deletion positions
    del_positions = []
    del_depths = []
    del_hover = []

    for row in freq_results:

        pos = row[0]

        A_freq = row[1]
        C_freq = row[2]
        G_freq = row[3]
        T_freq = row[4]

        depth = row[5]

        del_freq = row[6]
        ins_freq = row[7]

        if ins_freq >= minAF and depth > mindepth:
            ins_positions.append(pos)
            ins_depths.append(depth)

            ins_hover.append([
                A_freq,
                C_freq,
                G_freq,
                T_freq,
                del_freq,
                ins_freq
            ])

        if del_freq >= minAF and depth > mindepth:
            del_positions.append(pos)
            del_depths.append(depth)

            del_hover.append([
                A_freq,
                C_freq,
                G_freq,
                T_freq,
                del_freq,
                ins_freq
            ])

    # Insertion markers
    if len(ins_positions) > 0:

        fig.add_trace(
            go.Scatter(
                x=ins_positions,
                y=ins_depths,
                mode='markers',
                name='Insertion (AF ≥ threshold)',
                marker=dict(
                    symbol='hexagram',
                    color='#FFC107',
                    size=12,
                    
                ),
                customdata=ins_hover,
                hovertemplate=
                    '<b>INSERTION</b><br>' +
                    'Position: %{x}<br>' +
                    'Depth: %{y}<br>' +
                    'A_freq: %{customdata[0]:.2f}<br>'+
                    'C_freq: %{customdata[1]:.2f}<br>'+
                    'G_freq: %{customdata[2]:.2f}<br>'+
                    'T_freq: %{customdata[3]:.2f}<br>'+
                    'Del_freq: %{customdata[4]:.2f}<br>'+
                    'Ins_freq: %{customdata[5]:.2f}<extra></extra>'
            )
        )

    # Deletion markers
    if len(del_positions) > 0:

        fig.add_trace(
            go.Scatter(
                x=del_positions,
                y=del_depths,
                mode='markers',
                name='Deletion (AF ≥ threshold)',
                marker=dict(
                    symbol='diamond',
                    color='#DC267F',
                    size=12,
                    
                ),
                customdata=del_hover,
                hovertemplate=
                    '<b>DELETION</b><br>' +
                    'Position: %{x}<br>' +
                    'Depth: %{y}<br>' +
                    'A_freq: %{customdata[0]:.2f}<br>'+
                    'C_freq: %{customdata[1]:.2f}<br>'+
                    'G_freq: %{customdata[2]:.2f}<br>'+
                    'T_freq: %{customdata[3]:.2f}<br>'+
                    'Del_freq: %{customdata[4]:.2f}<br>'+
                    'Ins_freq: %{customdata[5]:.2f}<extra></extra>'
            )
        )

    fig.update_layout(
        title='Coverage Plot with Indel Markers',
        xaxis_title='Position',
        yaxis_title='Depth',
        hovermode='closest',
        template='plotly_white',
        legend=dict(
            orientation='h',
            yanchor='bottom',
            y=1.02,
            xanchor='left',
            x=0
        )
    )

    fig.write_html(output_path)

def main():
    args = parser.parse_args()

    bamfile = pysam.AlignmentFile(args.bam, "rb")
    ref_name = bamfile.references[0]
    refseq = pysam.Fastafile(args.reference).fetch(ref_name)
    
    #Make an array of start-stop intervals to parallelize processing
    genome_length = len(refseq)
    split = int(genome_length/args.cores)
    batch = [[i*split, (i+1)*split, args.bam, args.maxdepth, args.reference] for i in range(args.cores)]

    #Adjust the last "stop" to be the genome length
    batch[-1][1] = genome_length

    with multiprocessing.Pool(processes=args.cores) as p:
        resultlist = p.starmap(process_batch, iter(batch))

    variant_results = list(itertools.chain.from_iterable(r[0] for r in resultlist))
    freq_results = list(itertools.chain.from_iterable(r[1] for r in resultlist))

    # Build variant and frequency dictionaries for consensus
    variant_dict = {}
    freq_dict = {}
    for result in variant_results:
        num_aln, ref_pos, ref_seq, alt_seq, alt_count, alt_AF = result
        variant_dict[ref_pos] = result
    for row in freq_results:
        freq_dict[row[0] - 1] = row

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

    # Create coverage plot
    if args.coverageplot:
        create_coverage_plot(freq_results, args.coverageplot, args.minAF, args.mindepth)
    
    # Build consensus sequence by walking through genome positions
    consensus = []
    pos = 0
    while pos < genome_length:
        if pos in variant_dict:
            variant_row = variant_dict[pos]
            freq_row = freq_dict.get(pos)
            if freq_row is None:
                consensus.append('N')
                pos += 1
                continue
            
            # Use frequency data and criteria to select consensus base
            consensus_base = select_consensus_allele(
                variant_row,
                freq_row,
                args.ambiguous,
                args.minAF,
                args.mindepth,
                args.keepindels
            )

            if consensus_base == 'N':
                consensus.append('N')
                pos += 1
                continue

            # Small indels are skipped by returning reference sequence
            if allele_is_indel(variant_row[2], variant_row[3]) and abs(len(variant_row[2]) - len(variant_row[3])) in [1, 2] and not args.keepindels:
                consensus.append(variant_row[2])
            else:
                consensus.append(consensus_base)

            # Skip positions covered by deletions in the consensus allele
            if len(variant_row[2]) > len(variant_row[3]):
                pos += (len(variant_row[2]) - len(variant_row[3]))
        else:
            # Position with no variant calls defaults to 'N'
            consensus.append("N")
        pos += 1

    consensus = ''.join(consensus)

    # Write consensus to FASTA format
    with open(args.out, 'w') as out_consensus:
        print(''.join(['>',args.outname]), file=out_consensus)
        print(consensus, file=out_consensus)


if __name__ == '__main__':
    main()
