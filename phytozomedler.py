# Nick Weiner 2017
# phytozomedler.py
# Getting data from the phytozome biomart using edited bioservices'
#  biomart module.
# Use help flag for usage statement.
from bioservices import *
from math import isnan
# For if we want to to jcvi's biomart
# from jcvi import biomart
# s = Mart(host="phytozome.jgi.doe.gov" name="phytozome")

s = BioMart(secure=True)
ribosomal_kegg_IDs = \
    [
        'K01977', 'K01979', 'K01980', 'K01981', 'K01982', 'K01985', 'K01986',
        'K02863', 'K02864', 'K02865', 'K02866', 'K02867', 'K02868', 'K02869',
        'K02870', 'K02871', 'K02872', 'K02873', 'K02874', 'K02875', 'K02876',
        'K02877', 'K02878', 'K02879', 'K02880', 'K02881', 'K02882', 'K02883',
        'K02884', 'K02885', 'K02886', 'K02887', 'K02888', 'K02889', 'K02890',
        'K02891', 'K02892', 'K02893', 'K02894', 'K02895', 'K02896', 'K02897',
        'K02898', 'K02899', 'K02900', 'K02901', 'K02902', 'K02903', 'K02904',
        'K02905', 'K02906', 'K02907', 'K02908', 'K02909', 'K02910', 'K02911',
        'K02912', 'K02913', 'K02914', 'K02915', 'K02916', 'K02917', 'K02918',
        'K02919', 'K02920', 'K02921', 'K02922', 'K02923', 'K02924', 'K02925',
        'K02926', 'K02927', 'K02928', 'K02929', 'K02930', 'K02931', 'K02932',
        'K02933', 'K02934', 'K02935', 'K02936', 'K02937', 'K02938', 'K02939',
        'K02940', 'K02941', 'K02942', 'K02943', 'K02944', 'K02945', 'K02946',
        'K02947', 'K02948', 'K02949', 'K02950', 'K02951', 'K02952', 'K02953',
        'K02954', 'K02955', 'K02956', 'K02957', 'K02958', 'K02959', 'K02960',
        'K02961', 'K02962', 'K02963', 'K02964', 'K02965', 'K02966', 'K02967',
        'K02968', 'K02969', 'K02970', 'K02971', 'K02973', 'K02974', 'K02975',
        'K02976', 'K02977', 'K02978', 'K02979', 'K02980', 'K02981', 'K02982',
        'K02983', 'K02984', 'K02985', 'K02986', 'K02987', 'K02988', 'K02989',
        'K02990', 'K02991', 'K02992', 'K02993', 'K02994', 'K02995', 'K02996',
        'K02997', 'K02998', 'K07590'
    ]

filtSelections = {  # Protip: DO NOT put spaces into these strings!
                    'none': [],
                    'ubiquitin': ['pfam_id_list', 'PF00240'],
                    'arabidopsis': ['organism_id', '447'],
                    # early release Feb 2017
                    'populus':  ['organism_id', '445,444'],  # poplars
                    'fiveprime': ['5_utr', 'excluded="0"'],
    # random filter for orthologs
                    'b_orth': ["ortholog_organism_name", "Bstricta"],
    # random arabidopsis gene
                    'atgene': ["gene_name_filter", "AT3G11260"],
                    'ribosome': ['kegg_orth_list', ','.join(ribosomal_kegg_IDs)]
                }

attrSelections = {
    # ALMOST ALL info
    'unsplicedTranscript':
        {
            'data_types':
                ['gene_name1', 'transcript_name1', 'transcript_exon_intron',
                 'organism_name', 'chr_name1', 'gene_chrom_start',
                 'gene_chrom_end', 'gene_chrom_strand', 'transcript_id',
                 'transcript_chrom_start', 'transcript_chrom_end',
                 'peptide_name', 'exon_chrom_start', 'exon_chrom_end',
                 'exon_cds_start', 'exon_cds_end', '5_utr_start', '5_utr_end',
                 '3_utr_start', '3_utr_end'],
            'data_positions': [8, 9, float('nan'), 10, 11, 12, 13, 14, 15, 16,
                               17, 18, 0, 1, 2, 3, 4, 5, 6, 7],
            'formats': ['FASTA'],
        },
    'peptideLite':
        {
            'data_types':
                ['peptide_sequence', 'peptide_name', 'transcript_name1',
                 'organism_name', 'gene_name1', 'transcript_id', ],
            'data_positions': [float('nan'), 0, 1, 2, 3, 4],
            'formats': ['FASTA'],
        },
    'allFeatures':
        {
            'data_types':
                [
                    'gene_name1', 'transcript_name1', 'organism_name',
                    'organism_id', 'gene_description', 'chr_name1',
                    'gene_chrom_strand', 'gene_chrom_start', 'gene_chrom_end',
                    'transcript_id', 'peptide_name', 'transcript_chrom_start',
                    'transcript_chrom_end', 'pfam_id', 'pfam_desc', 'smart_id',
                    'smart_desc', 'panther_id', 'panther_desc', 'pathway_id',
                    'pathway_desc', 'kog_id', 'kog_desc', 'kegg_enzyme_id',
                    'kegg_enzyme_desc', 'ko_id', 'keggorth_desc', 'go_id',
                    'go_desc', 'embl_id', 'entrez_gene_id', 'unigene_id',
                    'refseq_id', 'sptrembl_id', 'synonyms',
                ],
            'formats': ['TSV', 'HTML', 'CSV', 'GFF', 'XLS']
        },
    'allStructures':
        {
            'data_types':
                [
                    'organism_name', 'gene_name', 'chr_name',
                    'transcript_name', 'gene_description', 'gene_chrom_start',
                    'gene_chrom_end', 'gene_chrom_strand', 'peptide_name',
                    'transcript_chrom_start', 'transcript_chrom_end',
                    'exon_chrom_start', 'exon_chrom_end', 'exon_cds_start',
                    'exon_cds_end', 'exon_rank', 'exon_phase',
                ],
            'formats': ['TSV', 'HTML', 'CSV', 'GFF', 'XLS']
        },
    'orthologs':
        {
            'data_types':
                [
                    'ortholog_group', 'ortholog__dm_gene_name',
                    'ortholog__dm_organism_name',
                    'ortholog__dm_ortholog_gene_name',
                    'ortholog__dm_ortholog_organism_name',
                    'ortholog__dm_relationship'
                ],
            'formats': ['TSV', 'HTML', 'CSV', 'GFF', 'XLS']
        },
    # Add more sets here
    # set2:
}


def test_selections(filts, attrs, form):
    filts = filts.split(',')
    for f in filts:
        if f not in filtSelections.keys():
            raise NameError('filter {} not available'.format(f))
    if attrs not in attrSelections.keys():
        raise NameError('attribute set {} not available'.format(attrs))
    if form not in attrSelections[attrs]['formats']:
        raise NameError('{} does not have the {} format available'.format(
                attrs, form))
    return filts, attrs


def make_my_xml(filters, attributes, form):
    # s.secure = True
    s.host = "phytozome.jgi.doe.gov"
    # s._set_host(,True)
    s.custom_query(virtualScheme="zome_mart", formatter=form, unique=1)
    s.add_dataset_to_xml('phytozome')
    filters, attributes = test_selections(filters, attributes, form)
    for f in filters:
        print (filtSelections[f])
        if f != 'none':
            s.add_filter_to_xml(*filtSelections[f])
    for a in attrSelections[attributes]['data_types']:
        s.add_attribute_to_xml(a)
    return s.get_xml()


def read_my_xml(my_file):
    f = open(my_file)
    xml_query = f.read()
    f.close()
    return xml_query


def get_results(query):
    # if qury.endswith('.xml')
    # make this function also compatible with directly taking in xml file name?
    print ("Querying " + s.host + "\n")
    res = s.query(query)
    return res


def replace_headers(fasta, attributeset):
    import re

    def format_headers(matchobj):
        vs = matchobj.group(1).split('|')  # values
        out_string = '>'
        # put all of the header data into xml style tags
        for i in range(0, len(attrSelections[attributeset]['data_positions'])):
            if not isnan(attrSelections[attributeset]['data_positions'][i]):
                out_string += '{}="{}" '.format(
                    attrSelections[attributeset]['data_types'][i],
                    vs[attrSelections[attributeset]['data_positions'][i]])
        return out_string

    new_out = re.sub('^>(.+)', format_headers, fasta, flags=re.M)
    return new_out


def export_results(result, outname):
    o = open(outname, 'w')
    o.write(result)
    o.close()


def sample_out(result, form):
    # The result's data type is usually "unicode"
    print ('Your output looks like:\n')
    print ("length: " + str(len(result)) + '\n')
    if form != 'XLS':
        res = result.split('\n')
        for x in res[0:25]:
            print (x)
    else:
        # Maybe add this feature later
        print("Phytozomedler does not yet support XLS previews.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""Get data from phytozome 
    biomart with select xml elements
    \n example code: python phytozomedler arabidopsis,fiveprime
    unsplicedTranscript mfo.FASTA""")
    parser.add_argument("filterSet",
        help="Select preloaded filters separated by commas:\n{}".format(
            ', '.join(filtSelections)))
    parser.add_argument("attributeSet",
        help="Select preloaded attribute set:\n{}".format(
            ', '.join(attrSelections)))
    parser.add_argument("outputf",  # handle both output file name and format
                        help="output file name with format, e.g. output.FASTA")
    parser.add_argument('--verbose', '-v', default=0, action='count',
                        help='Multiple flags increase verbosity')
    parser.add_argument('-test', '-t', action='store_true',
                        help='Do not store the data')
    args = parser.parse_args()

    # interp input
    fts, ats = args.filterSet, args.attributeSet
    fmt = args.outputf.split('.')[1]
    if args.verbose > 1:
        print ("selected filters: {}\n".format(fts))
        print ("selected attributes: {} with {} format\n".format(ats, fmt))
    myQuery = make_my_xml(fts, ats, fmt)
    if args.verbose > 2:
        print (myQuery)
    resulting_data = get_results(myQuery)

    # output file
    if fmt == 'FASTA':
        resulting_data = replace_headers(resulting_data, ats)

    if not args.test:
        export_results(resulting_data, args.outputf)
    # see what it looks like
    if args.verbose:
        sample_out(resulting_data, fmt)
