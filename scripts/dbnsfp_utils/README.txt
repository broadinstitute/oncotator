Please follow the steps below in order to create a dbNSFP data source:
1. Download zipped database from https://sites.google.com/site/jpopgen/dbNSFP (example: http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFPv2.4.zip).
2. Unzip using java's jar command (example: jar xvf dbNSFPv2.4.zip)
3. Use concatenatedbNSFPTsvs.py to concatenate the tsv files (example: python concatenatedbNSFPTsvs.py --dir ~/dbSNFP --prefix dbNSFP2.4_variant.chr --out ~/dbNSFP2.4_variant.tsv)
4. Use initializeDatasource to build the data source (example: )