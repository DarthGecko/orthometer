# orthometer at korflab
**The project for getting orthologous introns and analyzing them for motifs**

[![Build Status](https://travis-ci.org/DarthGecko/orthometer.svg?branch=master)](https://travis-ci.org/DarthGecko/orthometer)

So far we have:
* A way to download data from Phytozome's biomart with Python.
* Some data files to work with.
* Multiple format options
* Scores 2mers of introns
* Includes peptide chain lengths
* Intron file production with introniator.py which scores acceptor and donor sites.

Working on:
* Dealing with lowercase and unambiguous sequences (see line 195 of intronitator)
* Giving FASTA headers JSON formats
* Motif finding algorithms
* Build up species list through queries
* Make a separate yml file of filters and attributes for each biomart host or database.

