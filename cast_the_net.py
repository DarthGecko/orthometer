# Nick Weiner 2017
# cast_the_net.py
# Determine the extent of available homology data

import phytozomedler
import random
org_id = '447' # arabidopsis
path = 'kegg_{}.TSV'.format(org_id)
try:
    with open(path) as f:
        data = f.readlines()
except OSError:
    base_species_query = \
        phytozomedler.make_my_xml([['organism_id', org_id], ['has_ko',
                                   'excluded="0"']], ['ko_id'],
                                  'TSV', 'af')
# 'kog_id', 'kog_desc', 'kegg_enzyme_id', 'kegg_enzyme_desc',
    data_f = phytozomedler.get_results(base_species_query)
    data = data_f.split('\n')
    sorted_data = '\n'.join(sorted(data))
    phytozomedler.export_results(sorted_data, path)

entries = len(data) - 1
print(entries)
r = random.randint(0, entries)
print (r)
my_entry = data[r]
print(my_entry)
kegg_id = my_entry.strip()
net_query = \
    phytozomedler.make_my_xml([['kegg_orth_list', kegg_id]],
                              ['organism_name'], 'TSV', 'af')
kegg_specs = phytozomedler.get_results(net_query)
phytozomedler.export_results(kegg_specs, 'species_with_KEGG_{}.TSV'.format(kegg_id))
species_num = len(kegg_specs.split('\n')) - 1
print ('Kegg ID {} data is available for {} species on Phytozome.'.format(
    kegg_id, species_num))

