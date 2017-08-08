# Nick Weiner 2017
# introninator.py
# Getting introns from FASTA fiels made by phytozomedler
from Bio import SeqIO
import re


# Things to work out:
# How are we storing the data on the introns? In RAM or ROM?
# Make a FASTA file  of introns?
# Set parameters of the record using title2ids
# Reference pages:   http://biopython.org/DIST/docs/api/Bio.SeqIO-module.html
#                    http://biopython.org/DIST/docs/api/Bio.SeqIO.FastaIO-module.html

'''
However, you can supply a title2ids function to alter this:

    >>> def take_upper(title):
    ...     return title.split(None, 1)[0].upper(), "", title
    >>> with open("Fasta/dups.fasta") as handle:
    ...     for record in FastaIterator(handle, title2ids=take_upper):
    ...         print(record.id)
'''

def

def strip_introns(fasta, chrom_or_cds):  #want the chrom (refers to coordinates)
    records = SeqIO.parse(fasta, "fasta")
    for seq_record in records:
        exon_positions = {}

        pos = ['beg', 'end']
        for i in pos:
            r = re.match('exon_{}_{}="(.+)"'.format(chrom_or_cds, i),
                         seq_record.id).split(';')
            exon_positions[i] = r
        # intron fasta file?
        intron_count = len(exon_positions['beg']) - 1  # Is this right?
        print(repr(seq_record.seq))
        print(len(seq_record))
        intron_positions = {}
        intron_positions['beg'] = []
        intron_positions['end'] = []
        for i in xrange(0, intron_count-1):  # Strand is represented by 1 or -1

            intron_positions['beg'].append(exon_positions['end'][i] + 1)
            intron_positions['end'].append(exon_positions['beg'][i+1] - 1)
        introns = []
        for i in xrange(0, intron_count - 1):
            introns.append(seq_record.seq[intron_positions])


