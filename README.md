# orthometer at korflab
**The project for getting orthologous introns and analyzing them for motifs**

[![Build Status](https://travis-ci.org/DarthGecko/orthometer.svg?branch=master)](https://travis-ci.org/DarthGecko/orthometer)

So far we have:
* A way to download data from Phytozome's biomart with Python.
* Some data files to work with.
* Multiple format options
* Intron file production with introniator.py which scores acceptor and donor sites.

Working on:
* Giving FASTA headers JSON formats
* Scoring 2mers of introns
* Motif finding algorithms
* Build up species list through queries
* Make a separate yml file of filters and attributes for each biomart host or database.

