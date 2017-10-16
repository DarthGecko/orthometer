# Nick Weiner 2017
# cast_the_net.py
# Determine the extent of available homology data

import phytozomedler
import random

my_query = \
    phytozomedler.make_my_xml([['organism_id', '447'], ['has_ko', 'excluded="0"']],
                               ['gene_name1', 'transcript_name1', 'kog_id',
                                'kog_desc', 'kegg_enzyme_id', 'kegg_enzyme_desc',
                                'ko_id'], 'GFF', 'af')

data = phytozomedler.get_results(my_query)
phytozomedler.export_results(data, 'arab_kegg.GFF')

r = random.randint(0, len(data.split('\n')))
print (r)
print(data.split('\n')[r]) #.split('\t')[3]
