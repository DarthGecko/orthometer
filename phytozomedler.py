# Nick Weiner 2017
# phytozomedler.py
# Getting data from the phytozome biomart using edited bioservices' biomart module.
# Use help flag for usage statement.
from bioservices import *
from math import isnan
s = BioMart()
s.host = "phytozome.jgi.doe.gov"

filtSelections = { #protip: DO NOT put spaces into these strings!
                    'ubiquitin': ['pfam_id_list', 'PF00240'],
                    'arabidopsis': ['organism_id','447','phytozome'],  # early release Feb 2017
                    'fiveprime': ['5_utr','excluded="0"'],
                    'ribosome': ['kegg_orth_list','K01977,K01979,K01980,K01981,K01982,K01985,K01986,K02863,K02864,K02865,K02866,K02867,K02868,K02869,K02870,K02871,K02872,K02873,K02874,K02875,K02876,K02877,K02878,K02879,K02880,K02881,K02882,K02883,K02884,K02885,K02886,K02887,K02888,K02889,K02890,K02891,K02892,K02893,K02894,K02895,K02896,K02897,K02898,K02899,K02900,K02901,K02902,K02903,K02904,K02905,K02906,K02907,K02908,K02909,K02910,K02911,K02912,K02913,K02914,K02915,K02916,K02917,K02918,K02919,K02920,K02921,K02922,K02923,K02924,K02925,K02926,K02927,K02928,K02929,K02930,K02931,K02932,K02933,K02934,K02935,K02936,K02937,K02938,K02939,K02940,K02941,K02942,K02943,K02944,K02945,K02946,K02947,K02948,K02949,K02950,K02951,K02952,K02953,K02954,K02955,K02956,K02957,K02958,K02959,K02960,K02961,K02962,K02963,K02964,K02965,K02966,K02967,K02968,K02969,K02970,K02971,K02973,K02974,K02975,K02976,K02977,K02978,K02979,K02980,K02981,K02982,K02983,K02984,K02985,K02986,K02987,K02988,K02989,K02990,K02991,K02992,K02993,K02994,K02995,K02996,K02997,K02998,K07590'], #ribosomal genes
                    'ribosome1': ['kegg_orth_list',
                 'K01977,K01979,K01980,K01981,K01982,K01985,K01986,K02863,K02864,K02865,K02866,K02867,K02868,K02869,K02870,K02871,K02872,K02873,K02874,K02875,K02876,K02877,K02878,K02879,K02880,K02881,K02882,K02883,K02884,K02885,K02886,K02887,K02888,K02889,K02890,K02891,K02892,K02893,K02894,K02895,K02896,K02897,K02898,K02899,K02900,K02901,K02902,K02903,K02904,K02905,K02906,K02907,K02908,K02909,K02910,K02911,K02912,K02913,K02914,K02915,K02916,K02917,K02918,K02919,K02920,K02921,K02922,K02923,K02924,K02925,K02926,K02927,K02928,K02929,K02930,K02931,K02932,K02933,K02934,K02935,K02936,K02937,K02938,K02939,K02940,K02941,K02942,K02943,K02944,K02945,K02946,K02947,K02948,K02949,K02950,K02951,K02952,K02953,K02954,K02955,K02956,K02957,K02958,K02959,K02960,K02961,K02962,K02963,K02964,K02965,K02966,K02967,K02968,K02969,K02970,K02971,K02973,K02974,K02975,K02976,K02977,K02978,K02979,K02980,K02981,K02982,K02983,K02984,K02985,K02986,K02987,K02988,K02989,K02990,K02991,K02992,K02993,K02994,K02995,K02996,K02997,K02998,K07590'],
# ribosomal genes
                }

attrSelections ={
    # ALMOST ALL info
    'unsplicedTranscript': [['gene_name1', 'transcript_name1', 'transcript_exon_intron', 'organism_name', 'chr_name1', 'gene_chrom_start', 'gene_chrom_end', 'gene_chrom_strand', 'transcript_id', 'transcript_chrom_start', 'transcript_chrom_end', 'peptide_name', 'exon_chrom_start', 'exon_chrom_end', 'exon_cds_start', 'exon_cds_end', '5_utr_start', '5_utr_end', '3_utr_start', '3_utr_end'],
            [8, 9, float('nan'), 10, 11, 12, 13, 14, 15, 16, 17, 18, 0, 1, 2, 3, 4, 5, 6, 7]]
    # 8, 9, sequence, 10, 11, 12, 13, 14, 15, 16, 17, 18, 0, 1, 2, 3, 4, 5, 6, 7
    # set2:
    
}


def make_example():
    ret = s.registry()
    s.names
    s.datasets('phytozome_mart')
    s.add_dataset_to_xml('phytozome')
    s.add_filter_to_xml("gene_name_filter","AT2G37550","phytozome")
    s.add_filter_to_xml('organism_id','167','phytozome')
    s.add_attribute_to_xml('gene_name1')
    s.add_attribute_to_xml('transcript_name1')
    s.add_attribute_to_xml('gene_exon_intron')
    return s.get_xml()


def test_selections(filts, attrs):
    filts = filts.split(',')
    for f in filts:
        if f not in filtSelections.keys():
            raise NameError('filter {} not available'.format(f))
    if attrs not in attrSelections.keys():
        raise NameError('attribute set {} not available'.format(attrs))
    return filts, attrs


def make_my_xml(filters, attributes):
    #s.host = "phytozome.jgi.doe.gov"
    #s.datasets('phytozome_mart')
    s.add_dataset_to_xml('phytozome')
    filters, attributes = test_selections(filters, attributes)
    for f in filters:
        print filtSelections[f]
        s.add_filter_to_xml(*filtSelections[f])

    for a in attrSelections[attributes][0]:
        s.add_attribute_to_xml(a)

    return s.get_xml()


# <Dataset name = "phytozome" interface = "default" >
# 		<Filter name = "gene_name_filter" value = "AT2G37550"/>
# 		<Filter name = "organism_id" value = "167"/>
# 		<Attribute name = "gene_name1" />
# 		<Attribute name = "transcript_name1" />
# 		<Attribute name = "transcript_exon_intron" />
def read_my_file(myfile):
    f = open(myfile)
    xml_query = f.read()
    return xml_query


def get_results(qury):
    #if qury.endswith('.xml')
    #make this function also compatible with directly taking in an xml file name
    print "Querying " + s.host + "\n"
    res = s.query(qury)
    return res


def replace_headers(fasta, attributeset):
    import re

    def format_headers(matchobj):
        vs = matchobj.group(1).split('|')
        out_string = '>'
        #put all of the header data into xml style tags
        for i in range(0, len(attrSelections[attributeset][1])):
            if not isnan(attrSelections[attributeset][1][i]):
                out_string = '{}{}="{}" '.format(out_string, attrSelections[attributeset][0][i], vs[attrSelections[attributeset][1][i]])
        return out_string

    new_out = re.sub('^>(.+)', format_headers, fasta, flags=re.M)
    return new_out
    # for seq_record in SeqIO.parse(fasta, "fasta"):
    #     i = i + 1
    #     print(seq_record.id)
    #     print(repr(seq_record.seq))
    #     print(len(seq_record))
    #     if i >3:
    #         break


    # for line in fasta.splitlines():
    #     if line[0] == '>':
    #         values = line.split('|')
    #          in values


def export_results(result, outname, attributeset):
    o = open(outname, 'w')
    o.write(result)
    o.close()


def sample_out(result):
    # The result's data type is usually "unicode"
    print 'Your output looks like:\n'
    print "length: " + str(len(result))
    res = result.split('\n')
    for x in res[0:25]:
        print x

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""Get data from phytozome biomart with select xml elements
    \n example code: python phytozomedler arabidopsis,fiveprime unsplicedTranscript""")
    parser.add_argument("filterSet", help="Select preloaded filters separated by commas:\n{}".format(filtSelections.keys()))
    parser.add_argument("attributeSet", help="Select preloaded attribute set:\n{}".format(attrSelections.keys()))
    parser.add_argument("outputf", nargs="?", help = "Optional output file name", default="outputf.txt")
    parser.add_argument('--verbose', '-v', action='count', help='Multiple flags increase verbosity')
    parser.add_argument('-test', '-t', action='store_true', help='Do not store the data')
    args = parser.parse_args()

    #interp input
    fts, ats = args.filterSet, args.attributeSet
    if args.verbose > 1:
        print "selected filters: {}\n".format(fts)
        print "selected attributes: {}\n".format(ats)
    myQuery = make_my_xml(fts, ats)
    if args.verbose > 2:
        print myQuery
    resulting_data = get_results(myQuery)

    ## output file
    resulting_data = replace_headers(resulting_data, ats)
    if not args.test:
        export_results(resulting_data, args.outputf, ats)
    # see what it looks like
    if args.verbose:
        sample_out(resulting_data)
