#!/usr/bin/env python3
"""
My MinION data contains a fairly high number of multimeric reads where PCR
products appear to be physically joined end-to-end. This creates trouble for
standard demulitplexing, but tossing multimers is a waste since the multimeric
reads contain usable information. As long as barcodes for constituent reads are
intact, multimers can be decomposed and salvaged.

Now with multiprocessing!

Requires BioPython and regex (not re!)

CL Weinrich 2026
"""

import argparse
import regex
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count
from os import makedirs
from os.path import join
from Bio import SeqIO
from Bio.Seq import Seq
from itertools import product

# Globals used by worker processes
_top_barcode_to_ID = {}
_bot_barcode_to_ID = {}
_top_F_seqs = set()
_top_R_seqs = set()
_bot_F_seqs = set()
_bot_R_seqs = set()

# Sets for degenerate base conversions used by find_barcodes()
# Not sure why you would ever have degenerate barcodes...
IUPAC_dict = {
    frozenset(['G']): 'G',
    frozenset(['A']): 'A',
    frozenset(['T']): 'T',
    frozenset(['C']): 'C',
    frozenset(['G', 'A']): 'R',
    frozenset(['T', 'C']): 'Y',
    frozenset(['A', 'C']): 'M',
    frozenset(['G', 'T']): 'K',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['A', 'C', 'T']): 'H',
    frozenset(['G', 'T', 'C']): 'B',
    frozenset(['G', 'C', 'A']): 'V',
    frozenset(['G', 'A', 'T']): 'D',
    frozenset(['G', 'A', 'T', 'C']): 'N'
}
IUPAC_reverse = {v: k for k,v in IUPAC_dict.items()}
IUPAC_dict[frozenset(['N'])] = 'N'


def read_barcodes(file, order, out_dir, suffix, delimiter='\t'):
    """
    Read primer and barcode combinations from text file
    
    :param file: Text file containing primers, barcodes, and ID values
    :param delimiter: Character delimiting columns in primer/barcode file
    :param order: Comma separated string of integers designating the column index of 
                    ID, F.barcode, F.primer, R.barcode, and R.primer 
                    in primer/barcode file
    """

    # Get info from the correct column and place into dict
    # bot = RC of top strand
    top_barcode_to_ID, bot_barcode_to_ID = {}, {}
    top_F_seqs, top_R_seqs, bot_F_seqs, bot_R_seqs = set(), set(), set(), set()
    order = [(int(i) - 1) for i in str(order).split(',')]

    with open(file, 'r') as f:
        for row in f:
            row_list = row.strip('\n').split(delimiter)
            # TODO: utilize primers to remove false positive barcode hits?
            ID, top_F_barcode, top_F_primer, bot_R_barcode, bot_R_primer = [row_list[o] for o in order]

            # Create empty file for each sample (prevents snakemake workflow breaking if there are no seqs for some samples)
            file_path = join(out_dir, f'{ID}{suffix}.fastq') if out_dir else f'{ID}{suffix}.fastq'
            with open(file_path, 'w') as f:
                pass

            # R seqs assumed to be given in standard 5'-3' bottom strand orientation 
            top_F_seq = top_F_barcode
            bot_R_seq = bot_R_barcode

            # RC 
            bot_F_seq = str(Seq(top_F_seq).reverse_complement())
            top_R_seq = str(Seq(bot_R_seq).reverse_complement())

            # barcode combinations to ID
            if (top_F_seq, top_R_seq) in top_barcode_to_ID:
                raise KeyError(f"Duplicate barcode combination detected. Please double check the barcode file: {top_F_seq, top_R_seq}.")
            else:
                # Top strand orientation
                top_F_seqs.add(top_F_seq)
                top_R_seqs.add(top_R_seq)
                top_barcode_to_ID[(top_F_seq, top_R_seq)] = ID

                # Bottom strand orientation
                bot_F_seqs.add(bot_F_seq)
                bot_R_seqs.add(bot_R_seq)
                bot_barcode_to_ID[(bot_R_seq, bot_F_seq)] = ID # order of seqs swaps on bottom strand

    return top_barcode_to_ID, bot_barcode_to_ID, top_F_seqs, top_R_seqs, bot_F_seqs, bot_R_seqs


def find_barcode(seq, barcode, max_edit_dist, max_deletions, max_insertions, max_substitutions):
    """
    Search for pattern in sequence, allowing degenerate positions and mismatches.
    """
    # Establish constraints
    constraints = []
    if max_substitutions:
        constraints.append(f's<={max_substitutions}')
    if max_insertions:
        constraints.append(f'i<={max_insertions}')
    if max_deletions:
        constraints.append(f'd<={max_deletions}')
    constraints.append(f'e<={max_edit_dist}')
    
    # Convert nucleotides to a pattern that account for degeneracy
    new_pattern = ''
    for position in barcode:
        new_pattern += f'[{"".join(set([*IUPAC_reverse[position], position]))}]'
    pattern = f'({new_pattern}){{{",".join(constraints)}}}'
    matches = list(regex.finditer(pattern, seq))

    return matches


def map_barcodes(search_seqs, record_seq):
    """
    Return a dictionary of matches for search_seq within record_seq
    """
    hit_map = {}
    for search_seq in search_seqs:
        matches = find_barcode(record_seq, search_seq, max_edit_dist=_max_edit_dist, max_deletions=_max_deletions, max_insertions=_max_insertions, max_substitutions=_max_substitutions)
        if matches:
            hit_map[search_seq] = matches

    return hit_map


def find_monomers(F_hit_map, R_hit_map, barcode_dict, strand):
    """
    Find valid monomers defined as: 
        [F][EL][R] 
        Where:
            [F] = Forward barcode sequence
            [EL] = Sequence of expected length
            [R] = Reverse primer sequence
        And (F, R) is a valid barcode combination found in the supplied barcode_dict
    """
    monomers = []

    if strand == 'top':
        seq_combos = list(product(F_hit_map, R_hit_map))
    elif strand == 'bot':
        seq_combos = list(product(R_hit_map, F_hit_map))

    for combo in seq_combos:
        if combo in barcode_dict.keys():

            # get index positions for combos
            if _trim:
                if strand == 'top':
                    F, R = combo[0], combo[1]
                    F_hits = [h.end() for h in F_hit_map[F]]
                    R_hits = [h.start() for h in R_hit_map[R]]
                elif strand == 'bot':
                    F, R = combo[1], combo[0]
                    F_hits = [h.start() for h in F_hit_map[F]]
                    R_hits = [h.end() for h in R_hit_map[R]]
            else:
                if strand == 'top':
                    F, R = combo[0], combo[1]
                    F_hits = [h.start()for h in F_hit_map[F]]
                    R_hits = [h.end() for h in R_hit_map[R]]
                elif strand == 'bot':
                    F, R = combo[1], combo[0]
                    F_hits = [h.end() for h in F_hit_map[F]]
                    R_hits = [h.start()for h in R_hit_map[R]]

            # check that seq of expected length is between indices
            for F, R in list(product(F_hits, R_hits)):
                if strand == 'top':
                    length = R - F
                    if _expected_length[0] <= length <= _expected_length[1]:
                        monomers.append((barcode_dict[combo], (F, R), 'top'))
                elif strand == 'bot':
                    length = F - R
                    if _expected_length[0] <= length <= _expected_length[1]:
                        monomers.append((barcode_dict[combo], (R, F), 'bot'))

    return monomers


def init_worker(top_F_seqs_, top_R_seqs_, bot_F_seqs_, bot_R_seqs_, top_barcode_to_ID_, bot_barcode_to_ID_, trim_, expected_length_, max_edit_dist_, max_deletions_, max_insertions_, max_substitutions_):
    global _top_F_seqs, _top_R_seqs, _bot_F_seqs, _bot_R_seqs
    global _top_barcode_to_ID, _bot_barcode_to_ID
    global _trim, _expected_length, _max_edit_dist, _max_deletions, _max_insertions, _max_substitutions

    _top_F_seqs = top_F_seqs_
    _top_R_seqs = top_R_seqs_
    _bot_F_seqs = bot_F_seqs_
    _bot_R_seqs = bot_R_seqs_
    _top_barcode_to_ID = top_barcode_to_ID_
    _bot_barcode_to_ID = bot_barcode_to_ID_
    _trim = trim_
    _expected_length = expected_length_
    _max_edit_dist = max_edit_dist_
    _max_deletions = max_deletions_
    _max_insertions = max_insertions_
    _max_substitutions = max_substitutions_


def record_to_tuple(record):
    """Read FASTQ record into tuple"""
    return (
        record.id,
        str(record.seq),
        record.letter_annotations.get('phred_quality', []),
        record.description,
    )


def tuple_to_record(record_tuple):
    """Convert tuple back to FASTQ record"""
    record_id, seq_str, quality, description = record_tuple
    record = SeqIO.SeqRecord(Seq(seq_str), id=record_id, description=description)
    if quality:
        record.letter_annotations['phred_quality'] = quality
    return record


def process_record(record_tuple):
    """Extract and identify valid monomers in the record tuple"""
    record = tuple_to_record(record_tuple)
    record_seq = str(record.seq)

    top_F_hit_map = map_barcodes(_top_F_seqs, record_seq)
    top_R_hit_map = map_barcodes(_top_R_seqs, record_seq)
    bot_F_hit_map = map_barcodes(_bot_F_seqs, record_seq)
    bot_R_hit_map = map_barcodes(_bot_R_seqs, record_seq)

    top_monomers = find_monomers(top_F_hit_map, top_R_hit_map, _top_barcode_to_ID, 'top')
    bot_monomers = find_monomers(bot_F_hit_map, bot_R_hit_map, _bot_barcode_to_ID, 'bot')
    monomers = top_monomers + bot_monomers

    outputs = []
    for j, monomer in enumerate(monomers):
        ID, index_tup, strand = monomer
        new_record = record[index_tup[0]:index_tup[1]]
        new_record.id = f'{j}_{strand}_{record.id}'

        if strand == 'bot':
            new_record = new_record.reverse_complement(id=True, description=True)

        outputs.append((ID, new_record))

    return outputs, len(monomers)


def flush(write_buffer, out_dir, suffix):
    """Write idenitifed monomers to respective FASTQ files"""
    total_written = 0

    for ID, records in write_buffer.items():
        file_path = join(out_dir, f'{ID}{suffix}.fastq') if out_dir else f'{ID}{suffix}.fastq'
        with open(file_path, 'a') as out:
            SeqIO.write(records, out, 'fastq')
        total_written += len(records)
    write_buffer.clear()

    return total_written


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq')
    parser.add_argument('barcodes')
    parser.add_argument('--barcode_order', default='1,2,3,4,5')

    parser.add_argument('--suffix', default='')
    parser.add_argument('--out_dir', type=str, default='')

    parser.add_argument('--flush_size', type=int, default=1000)
    parser.add_argument('--workers', type=int, default=1)

    parser.add_argument('--trim', action='store_true')

    parser.add_argument('--min_length', type=int, default=100)
    parser.add_argument('--max_length', type=int, default=500)

    parser.add_argument('--max_edit_dist', type=int, default=0)
    parser.add_argument('--max_deletions', type=int, default=None)
    parser.add_argument('--max_insertions', type=int, default=None)
    parser.add_argument('--max_substitutions', type=int, default=None)

    return parser.parse_args()


def main():
    args = parse_args()

    fastq = args.fastq
    barcodes = args.barcodes
    barcode_order = args.barcode_order
    
    suffix = args.suffix
    out_dir = args.out_dir

    if out_dir:
        makedirs(out_dir, exist_ok=True)

    flush_size = args.flush_size
    workers = args.workers if args.workers > 0 else cpu_count()

    # More worker globals
    _trim = args.trim
    _expected_length = (args.min_length, args.max_length)
    _max_edit_dist = args.max_edit_dist
    _max_deletions = args.max_deletions
    _max_insertions = args.max_insertions
    _max_substitutions = args.max_substitutions

    # Get barcodes from file
    top_barcode_to_ID, bot_barcode_to_ID, top_F_seqs, top_R_seqs, bot_F_seqs, bot_R_seqs = read_barcodes(barcodes, barcode_order, out_dir, suffix)

    # Prepare some for recording some stats
    hit_stats = {}
    write_buffer = {}
    buffer_records = 0

    # Do the thing (in parallel!)
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=init_worker,
        initargs=(
            top_F_seqs,
            top_R_seqs,
            bot_F_seqs,
            bot_R_seqs,
            top_barcode_to_ID,
            bot_barcode_to_ID,
            _trim,
            _expected_length,
            _max_edit_dist,
            _max_deletions,
            _max_insertions,
            _max_substitutions,
        ),
    ) as executor:
        record_iterator = (record_to_tuple(record) for record in SeqIO.parse(fastq, 'fastq'))
        for outputs, monomer_count in executor.map(process_record, record_iterator, chunksize=10):
            if outputs:
                for ID, new_record in outputs:
                    write_buffer.setdefault(ID, []).append(new_record)
                    buffer_records += 1
                    if buffer_records >= flush_size:
                        flushed = flush(write_buffer, out_dir, suffix)
                        buffer_records -= flushed

            hit_stats[monomer_count] = hit_stats.get(monomer_count, 0) + 1

    # Final flush
    if write_buffer:
        flush(write_buffer, out_dir, suffix)

    print(f'file: {fastq}\nsuffix:{suffix}\n')
    print(f'seqs: {sum(hit_stats.values())}')

    hit_stats = dict(sorted(hit_stats.items()))
    total_monomers = sum(i * x for i, x in hit_stats.items())
    print(f'total monomers: {total_monomers}')

    try:
        print(f'junk rate: {hit_stats[0] / sum(hit_stats.values()):.3%}')
    except KeyError:
        print('junk rate: 0.000%')

    for i, x in hit_stats.items():
        print(f'{i}:\t{x}')


if __name__ == '__main__':
    main()
