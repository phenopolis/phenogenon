Download from `ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/` 20190701
then run

```
zcat clinvar_20190701.vcf.gz | grep -v '#' | grep 'HP:' | grep 'OMIM' | grep 'Pathogenic' > pathogenic_clinvar.vcf && bgzip pathogenic_clinvar.vcf && tabix -p vcf pathogenic_clinvar.vcf.gz
```

The resulting vcf.gz file can be uploaded to https://cadd.gs.washington.edu/ for scoring (with annotation, so that we can get gene names)
