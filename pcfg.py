# Nick Weiner 2017
# pcfg.py
# Get probabilisitic context free grammar from intron files made by intronitator
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs

def recurse1(seq, parts):
    data_set = []
    for i in range(1, len(seq)):
        addition = [seq[:i]]
        for j in range(1,parts):
            addition.append(seq[i:(len(seq)-i)])
        data_set.append(addition)
    return data_set

def recurse2(seq, i):
    retun seq[:i]