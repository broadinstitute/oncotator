### How to Create the COSMIC Datasources for Oncotator

Three COSMIC datasources make up all of the COSMIC annotations.

Example invocations for v76 are below...

#### Create the basic COSMIC datasource (genomic position and gene-protein position)

```
# This is for v76
python oncotator-index.py cosmic v76/CosmicCompleteTargetedScreensMutantExport.tsv 76``

#  then manually edit the config file, copy files into a new datasource directory (``cosmic/hg19``), and call the ``addDatasourceMd5.py`` utility:
addDatasourceMd5 cosmic/hg19/
```


#### Create the gene datasource (``cosmic_tissue`` gene only):

This one takes a while....

```
 python createCosmicGeneTsv.py v76/CosmicCompleteTargetedScreensMutantExport.tsv CosmicCompleteTargetedScreensMutantExport.v76.geneHistology.tsv
 initializeDatasource --ds_type gene_tsv --ds_file CosmicCompleteTargetedScreensMutantExport.v76.geneHistology.tsv --name COSMIC_Tissue --version v76 --dbDir /home/lichtens/test_dbdir/ --ds_foldername cosmic_tissue --genome_build hg19 --index_columns gene
```

#### Create the fusion gene datasrouce (``cosmic_fusion`` gene only):

```
 python createCosmicFusionGeneTsv.py CosmicFusionExport.tsv CosmicFusionExport.v76.tsv
 initializeDatasource --ds_type gene_tsv --ds_file CosmicFusionExport.v76.tsv --name COSMIC_FusionGenes --version v76 --dbDir /home/lichtens/test_dbdir/ --ds_foldername cosmic_fusion --genome_build hg19 --index_columns gene
```
