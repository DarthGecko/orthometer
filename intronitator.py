# Nick Weiner 2017
# introninator.py
# Getting introns from FASTA fiels made by phytozomedler
from Bio import SeqIO
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import IUPAC
don_len = 5
acc_len = 8

# Things to work out:
# How are we storing the data on the introns? In RAM or ROM?
# Make a FASTA file  of introns?
# Set parameters of the record using title2ids
# Reference pages:
#           http://biopython.org/DIST/docs/api/Bio.SeqIO-module.html
#           http://biopython.org/DIST/docs/api/Bio.SeqIO.FastaIO-module.html

# Negative strands??!

def analyze_intron(intron_seq):
    from Bio.SeqUtils import GC
    gc = '{0:.2f}'.format(GC(intron_seq))
    ambiguous = ('W' or 'S' or 'M' or 'K' or 'R' or 'Y' or 'B' or 'D' or 'H' or 'V'
             or 'N' or '-') in intron_seq.upper()
    if ambiguous:
        ambig = 'ambig'
    else:
        ambig = 'unamb'
    len(intron_seq)
    return [gc, ambig, ]

def score_site(seq, model):

    assert isinstance(model, motifs.Motif)
    assert isinstance(seq, Seq)
    pssm = model.counts.normalize(pseudocounts=0.5).log_odds()
    # p = 1
    # for i in range(0, len(seq)):
    #     nt = seq[i]
    #     print(model.counts[nt, i])
    #     p *= model.counts[nt, i]
    # '{0:.2f} compared to {}'.format(log(p / 0.25 ** len(seq)),
    #                              pssm.calculate(seq))
    return '{0:.2f}'.format(pssm.calculate(seq))


def dinucs(seq):
    seq = str(seq)
    if len(seq) < (don_len + acc_len + 2):
        return []
    else:
        return [seq[i:(i + 2)] for i in range(don_len, len(seq) - acc_len - 2)]


def score_dinucleotides(seq, model):
    # sub score_2mer
    # {
    #     my($seq, $model) = @_;
    #     my $score = 0;
    #     for (my $i = 0; $i < length($seq) -1; $i++) {
    #         my $di = substr($seq, $i, 2);
    #         $score += log($model->{$di} / 0.0625);  # ha, because I'm lazy
    #     }
    #     return sprintf "%.1f", 100 * $score;
    # }
    from math import log
    score = 0
    for di in dinucs(seq):
        score += log(model[di] / 0.0625)
    return '{0:.2f}'.format(score * 100)


def pseudo_score_site2(seq, model):
    from random import randrange
    return randrange(2)


def get_exon_id(header):  # Gives each record.name the exon coords septed by |
    return (re.match('.+transcript_name1="([^"]+)"', header).group(1),
            '|'.join(re.findall('exon_chrom_[star|end]+="([\d;]+)"', header)),
            header)


def strip_introns(fasta, verb=None, test=False, min_intron_len=35,
                  max_intron_len=10000, multi_species=False):
    # want the chrom (refers to coordinates)
    intron_file = '{}_introns_1.FASTA'.format(fasta[:-6])
    headline = '# id chr beg end str n/m len gc ambig? seq\n'
    enough_introns = False

    don_motif = {}
    acc_motif = {}
    dinuc_motif = {}
    dinuc_dist = {}
    with open(fasta) as handle:
        o = open(intron_file, 'w')
        o.write(headline)
        example = 0


        don = {}
        acc = {}
        dinuc = {}
        for seq_record in SeqIO.FastaIO.FastaIterator(handle,
                                                      title2ids=get_exon_id):
            if verb:
                print ("Seq Record: " + seq_record.name)
            chrom = re.match('.+chr_name1="([^"]+)"',
                             seq_record.description).group(1)
            if 'scaffold' in chrom:
                if verb:
                    print ('Scaffolding skipped!')
                continue
            exon_positions = {}
            pos = ['beg', 'end']
            r = seq_record.name.split('|')
            for i in range(2):
                exon_positions[pos[i]] = [int(x) for x in r[i].split(';')]
            strand = int(re.match('.+gene_chrom_strand="([^"]+)"',
                                  seq_record.description).group(1))
            species = re.match('.+organism_name="([^"]+)"',
                               seq_record.description).group(1)
            if verb:
                print ('strand: ', strand)
            start = int(re.match('.+transcript_chrom_start="([^"]+)"',
                                 seq_record.description).group(1))

            intron_count = len(exon_positions['beg']) - 1  # Is this right?
            if verb:
                print ('Exons:')
                for i in range(0, intron_count + 1):
                    print (
                    '{} - b: {} e: {}'.format(i + 1, exon_positions['beg'][i],
                                              exon_positions['end'][i]))
                    # print ('There should be {} introns.'.format(intron_count))

            intron_positions = {'beg': [], 'end': []}
            if verb:
                print ('Introns: ')
            for i in range(1, intron_count+1):  # Strand represented by 1 or -1
                # if strand > 0:
                    intron_positions['beg'].append(exon_positions['end'][i-1]+1)
                    intron_positions['end'].append(exon_positions['beg'][i] - 1)
                # else:
                #     intron_positions['beg'].append(exon_positions['end'][i] + 1)
                #     intron_positions['end'].append(exon_positions['beg'][i-1]-1)
            if verb:
                for i in range(0, intron_count):
                    print ('{} - b: {} e: {}'.format(i+1, intron_positions['beg'][i],
                                                     intron_positions['end'][i]))

            # return intron_positions # Is this all I want? Won't work with
            #   per transcript loop

            introns = []

            for i in range(0, intron_count):
                # intron = ''
                if strand > 0:
                    intron = seq_record.seq[intron_positions['beg'][i] -
                                            start:intron_positions['end'][i] -
                                            start]
                else:
                    intron = seq_record.seq[intron_positions['beg'][i] -
                                            start:intron_positions['end'][i] -
                                                  start]
                    # intron = seq_record.seq[intron_positions['end'][i] -
                    #                         start:intron_positions['beg'][i] -
                    #                         start]
                    intron = intron.reverse_complement()
                introns.append(intron)
            if verb:
                print ('The introns of {} are '.format(seq_record.id))
                for x in introns:
                    print (str(x))
            # Gather further info for output

            strand = int(re.match('.+gene_chrom_strand="([^"]+)"',
                                  seq_record.description).group(1))
            if strand > 0:
                strand_sym = '+'
            else:
                strand_sym = '-'

            # Output
            s = 1
            if species not in don:
                don[species] = []
                acc[species] = []
                dinuc[species] = []
                dinuc_motif[species] = []
            for x in introns:
                # If intron is not anomalous...
                if not (len(x) > max_intron_len or len(x) < min_intron_len):
                    #  Setting up donor and acceptor tables
                    # upper is good???
                    don[species].append(x.upper()[:don_len])
                    acc[species].append(x.upper()[-acc_len:])
                    dinuc[species].extend(dinucs(x))

                beg = intron_positions['beg'][s-1]
                end = intron_positions['end'][s-1]
                l = abs(end - beg)
                intron_set = '{}/{}'.format(s, intron_count)
                order = [seq_record.id, species, chrom, str(beg), str(end),
                         strand_sym, intron_set, str(l)]
                order.extend(analyze_intron(x))
                order.append(str(x))
                o.write('\t'.join(order)+'\n')
                s += 1
            example += 1
            if example > 4 and test:
                break
    # delete output file if not enough_introns?
    o.close()

    for species in don:
        don_motif[species] = motifs.create(don[species])
        acc_motif[species] = motifs.create(acc[species])
       # dinuc_motif[species] = motifs.create(dinuc[species])
        dinuc_dist[species] = {}
        for di in dinuc[species]:
            try:
                dinuc_dist[species][di] += 1
            except KeyError:
                dinuc_dist[species][di] = 1

    with open(intron_file) as out1:
        intron_file_2 = '{}_introns_2.FASTA'.format(fasta[:-6])
        out2 = open(intron_file_2, 'w')
        headline = '# id chr beg end str n/m len gc ambig? don acc 2mer seq\n'
        out2.write(headline)
        lines = out1.readlines()
        good_ones = 0
        for line in lines[1:]:
            intron = line.split()[-1]
            if len(intron) > max_intron_len or len(intron) < min_intron_len:
                continue
            species = line.split()[1]
            good_ones += 1
            d = score_site(Seq(intron[:don_len],
                               don_motif[species].alphabet),
                           don_motif[species])
            a = score_site(Seq(intron[-acc_len:],
                               acc_motif[species].alphabet),
                           acc_motif[species])
            di_score = score_dinucleotides(intron, dinuc_dist[species])
            order = ('\t'.join(line.split()[:-1]), d, a, di_score, intron)
            out2.write('\t'.join(order)+'\n')
        out2.close()
        if len(lines) == 0:
            print ('Requires Python 3 for additional processing')
        else:
            print ('Processed {} good introns out of {}'.format(good_ones,
                                                                len(lines) - 1))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""Make intron fasta files""")
    parser.add_argument("file_name",
        help="fasta file input")
    parser.add_argument("outputf", nargs="?",
        help="Optional output file name", default="outputf.FASTA")
    parser.add_argument('--verbose', '-v', action='count', default=0,
                        help='Multiple flags increase verbosity')
    parser.add_argument('-test', '-t', action='store_true',
                        help='Do not store the data')
    # parser.add_argument('-ms', '-multispecies', action='store_true',
    #                    help='Account for multiple species in the FASTA')
    args = parser.parse_args()

    # interp input
    strip_introns(args.file_name, args.verbose, args.test,)


# gene_name1="AT2G32350" transcript_name1="AT2G32350.1" organism_name="Athaliana_Araport11" chr_name1="Chr2" gene_chrom_start="13734945" gene_chrom_end="13735788" gene_chrom_strand="1" transcript_id="37375937" transcript_chrom_start="13734945" transcript_chrom_end="13735788" peptide_name="37375937_peptide" exon_chrom_start="13734945;13735345" exon_chrom_end="13735263;13735788" exon_cds_start="13734979;13735345" exon_cds_end="13735263;13735788" 5_utr_start="13734945" 5_utr_end="13734978" 3_utr_start="" 3_utr_end=""