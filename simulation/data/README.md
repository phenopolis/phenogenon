Download from `ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/` 20190701

Annotated with gnomAD and CADD using the CGAT wgs pipeline (slightly modified)

then run

```
zcat clinvar_20190701.anno.vcf.gz | grep '^#' > pathogenic_clinvar.vcf && zcat clinvar_20190701.anno.vcf.gz | grep 'HP:' | grep 'OMIM' | grep -E "reviewed_by_expert_panel|_multiple_submitters" | grep -v '_conflicting_interpretations' | grep 'Pathogenic' >> pathogenic_clinvar.vcf && bgzip pathogenic_clinvar.vcf && tabix -p vcf pathogenic_clinvar.vcf.gz
```

The resulting vcf.gz file can be uploaded to https://cadd.gs.washington.edu/ for scoring (with annotation, so that we can get gene names)
