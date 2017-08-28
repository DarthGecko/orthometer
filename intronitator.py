# Nick Weiner 2017
# introninator.py
# Getting introns from FASTA fiels made by phytozomedler
from Bio import SeqIO
import re
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import IUPAC


# Things to work out:
# How are we storing the data on the introns? In RAM or ROM?
# Make a FASTA file  of introns?
# Set parameters of the record using title2ids
# Reference pages:
#           http://biopython.org/DIST/docs/api/Bio.SeqIO-module.html
#           http://biopython.org/DIST/docs/api/Bio.SeqIO.FastaIO-module.html

# Negative strands??!

'''
However, you can supply a title2ids function to alter this:

    >>> def take_upper(title):
    ...     return title.split(None, 1)[0].upper(), "", title
    >>> with open("Fasta/dups.fasta") as handle:
    ...     for record in FastaIterator(handle, title2ids=take_upper):
    ...         print(record.id)
'''


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


def get_exon_id(header):  # Gives each record.name the exon coords septed by |
    return (re.match('.+transcript_name1="([^"]+)"', header).group(1),
            '|'.join(re.findall('exon_chrom_[star|end]+="([\d;]+)"', header)),
            header)


def strip_introns(fasta, verb=None, test=False):  #want the chrom (refers to coordinates)
    with open(fasta) as handle:
        intron_file = '{}_introns.FASTA'.format(fasta[:-6])
        o = open(intron_file, 'w')
        o.write('# id chr beg end str n/m len gc ambig? seq\n')
        example = 0
        for seq_record in SeqIO.FastaIO.FastaIterator(handle,
                                                      title2ids=get_exon_id):
            if verb:
                print (seq_record.name)
            exon_positions = {}
            pos = ['beg', 'end']
            r = seq_record.name.split('|')
            for i in range(2):
                exon_positions[pos[i]] = [int(x) for x in r[i].split(';')]
            strand = int(re.match('.+gene_chrom_strand="([^"]+)"',
                                  seq_record.description).group(1))
            print ('strand: ', strand)
            start = int(re.match('.+transcript_chrom_start="([^"]+)"',
                                 seq_record.description).group(1))
            # intron fasta file?
            print ('Exons:')

            intron_count = len(exon_positions['beg']) - 1  # Is this right?
            for i in range(0, intron_count+1):
                print ('{} - b: {} e: {}'.format(i+1, exon_positions['beg'][i],
                                                 exon_positions['end'][i]))
            # print ('There should be {} introns.'.format(intron_count))
            intron_positions = {'beg': [], 'end': []}
            print ('Introns: ')
            for i in range(1, intron_count+1):  # Strand represented by 1 or -1
                # if strand > 0:
                    intron_positions['beg'].append(exon_positions['end'][i-1]+1)
                    intron_positions['end'].append(exon_positions['beg'][i] - 1)
                # else:
                #     intron_positions['beg'].append(exon_positions['end'][i] + 1)
                #     intron_positions['end'].append(exon_positions['beg'][i-1]-1)
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
            chrom = re.match('.+chr_name1="([^"]+)"',
                              seq_record.description).group(1)
            data = [chrom, strand_sym]

            # Output
            s = 1

            for x in introns:
                beg = intron_positions['beg'][s-1]
                end = intron_positions['end'][s-1]
                l = abs(end - beg)
                line = '{}\t{}\t{}\t{}\t{}\t{}/{}\t{}\t'.format(seq_record.id,
                                                            data[0], beg, end,
                                                            strand_sym, s,
                                                            intron_count, l) +\
                       '\t'.join(str(d) for d in analyze_intron(x))+'\t'+str(x)
                o.write(line+'\n')
                s += 1
            example += 1
            if example > 4 and test:
                break



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""Make intron fasta files""")
    parser.add_argument("file_name",
        help="fasta file input")
    parser.add_argument("outputf", nargs="?",
        help="Optional output file name", default="outputf.FASTA")
    parser.add_argument('--verbose', '-v', action='count',
                        help='Multiple flags increase verbosity')
    parser.add_argument('-test', '-t', action='store_true',
                        help='Do not store the data')
    args = parser.parse_args()

    # interp input
    strip_introns(args.file_name, args.verbose, args.test)


# gene_name1="AT2G32350" transcript_name1="AT2G32350.1" organism_name="Athaliana_Araport11" chr_name1="Chr2" gene_chrom_start="13734945" gene_chrom_end="13735788" gene_chrom_strand="1" transcript_id="37375937" transcript_chrom_start="13734945" transcript_chrom_end="13735788" peptide_name="37375937_peptide" exon_chrom_start="13734945;13735345" exon_chrom_end="13735263;13735788" exon_cds_start="13734979;13735345" exon_cds_end="13735263;13735788" 5_utr_start="13734945" 5_utr_end="13734978" 3_utr_start="" 3_utr_end=""