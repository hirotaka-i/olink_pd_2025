# Introduction

used this [repo](https://github.com/hirotaka-i/new_analysis.git)


# Data required (data folder)
- `OLINK Disaease Blood Atlas.xlsx` - Olink data
- `olink-gene-mapping.tsv` - Olink gene mapping


# Codes
perform_pca.py: input-`OLINK Disaease Blood Atlas.xlsx` output-`olink_pca*.png`, `temp/ds.csv`
```
Rscript code/DEP.r\
 --data "temp/ds.csv"\
 --protein_start 22\
 --explanatory_vars "diagnosis,Sex,Age_at_collection,LEDD,pPC1,pPC2,pPC3"\
 --contrast "diagnosisPD - diagnosisControl"\
 --log2_transform FALSE

Rscript code/DEP.r\
 --data "temp/ds.csv"\
 --protein_start 22\
 --explanatory_vars "LEDD_CAT,Sex,Age_at_collection,LEDD,LEDD_sq,pPC1,pPC2,pPC3"\
 --contrast "LEDD_CATDenovoPD - LEDD_CATControl;LEDD_CATPD - LEDD_CATControl;LEDD_CATPD - LEDD_CATDenovoPD"\
 --log2_transform FALSE
```