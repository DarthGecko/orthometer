# Nick Weiner 2017
# pcfg.py
# Get probabilisitic context free grammar from intron files made by intronitator
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
from collections import Counter


def kmers(seq, k):
    return [seq[i:(i + k)] for i in range(0, len(seq) - k)]


def score_kmers(seq, min_k=2, max_k=None):
    data = {}
    if max_k is None:
        max_k = int(round(len(seq)/2, 0))
    if max_k > len(seq) - min_k:
        max_k = len(seq) - min_k
    for k in range(min_k, max_k):
        # for kmer in kmers(seq, k):

        data[k] = Counter(seq[start:start+k] for start in range(len(seq) - k))
    return data

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""Make comparisons between sequences""")
    parser.add_argument("sequence_1",
        help="first sequence")
    parser.add_argument("sequence_2", help="second sequence")
    parser.add_argument("min_k", nargs="?",
                        help="Optional minimum k", default=2)
    parser.add_argument("max_k", nargs="?",
                        help="Optional maximum k", default=None)
    parser.add_argument('--verbose', '-v', action='count', default=0,
                        help='Multiple flags increase verbosity')
    parser.add_argument('-test', '-t', action='store_true',
                        help='Do not store the data')
    # parser.add_argument('-ms', '-multispecies', action='store_true',
    #                    help='Account for multiple species in the FASTA')
    args = parser.parse_args()
    d1 = score_kmers(args.sequence_1, min_k=args.min_k, max_k=args.max_k)
    d2 = score_kmers(args.sequence_2, min_k=args.min_k, max_k=args.max_k)
    print (d1)
    print (d2)
