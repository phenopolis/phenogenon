# Phenogenon
[![Build Status](https://travis-ci.com/phenopolis/phenogenon.svg?branch=master)](https://travis-ci.com/phenopolis/phenogenon)
[![Coverage Status](https://coveralls.io/repos/github/phenopolis/phenogenon/badge.svg?branch=master)](https://coveralls.io/github/phenopolis/phenogenon?branch=master)

Preprint available on [biorxiv](https://www.biorxiv.org/content/early/2018/07/11/367292).

Phenogenon is a method that combines: the power of Human Phenotype Ontology for describing patient phenotypes, gnomAD for estimating rare variant population frequency, and CADD for variant pathogenicity prediction.

Given:
* a genotype vcf file (such as `tests/data/ABCA4.anonymised.vcf.gz`)
* a list of unrelated individuals (such as `tests/data/test_unrelated.tsv`)
* HPO dictionary (such as `tests/data/new-hpo-hpo.json`)
* a CADD score file (such as `tests/data/CADD_ABCA4_SCN1A.vcf.gz`)
* a patient-hpo file (such as `tests/data/test_patients_hpo_snapshot.tsv`)

Phenogenon is able to uncover true gene to phenotype relations, such as:
* "ABCA4 – Macular dystrophy"
* "SCN1A – Seizures".

Additionally, it accurately infers mode of inheritance (moi), such as:
* recessive mode of inheritance in the case of the "ABCA4 – Macular dystrophy"
* dominant mode of inheritance with the "SCN1A – Seizures" relationship."




## Installation


### Standalone
External programs you need to install:
* tabix
* bcftools

```bash
python setup.py install
```

#### Docker

It is also possible to run Phenogenon with Docker

```bash
# After starting the docker daemon:
docker pull phenopolis/phenogenon:latest
docker run -v $(pwd -P):/phenogenon -it phenogenon:latest bash
```

You can also build a docker image from this repository with the following command:

```bash
docker build -t phenogenon .
```

## Running Phenogenon

Running `phenogenon -h` will display the following help text

```bash
$ phenogenon -h
Usage: phenogenon [options] arg1 arg2

Options:
  -h, --help           show this help message and exit
  --output=OUTPUT      output file name?
  --vcf_file=VCF_FILE  bgzipped and tabix-indexed vcf.gz
  --range=RANGE        genome range to calculate? e.g. 2:4000-6000
```

The genomic range is usually a gene, which can be retrieved from an ensembl gtf file.
**Note** please use human reference build version b37, as this is the version supported by both `gnomAD` (gnomAD now also supports hg38) and `CADD`.

The vcf file should include all the samples that would be tested in Phenogenon. In order to index the file, please ensure it is sorted by both chromosomes and locations.
```
bgzip -c 1.vcf > 1.vcf.gz && tabix -p vcf 1.vcf.gz
```

### Config file
Parameters of Phenogenon are set in `configure.cfg` in the root folder.
It has detailed explanations in comments.

### Output
Run Unit test:
```bash
coverage run --omit=*/site-packages/*,*/tests/* -m unittest discover -s tests
```

##### Example of using Phenogenon:
```bash
phenogenon  --range 1:94458394-94586689 --vcf_file tests/data/ABCA4.anonymised.vcf.gz --output ABCA4.test.json
```
Produces `ABCA4.test.json`.

Explanation of output:
```python
{
  # "r": recessive moi, "d": dominant moi
  # cadd_15_ratio [0,1] look at how effective (relatively) variants with CADD phred score >=15 contribute to hgf
  # 1: very effective; 0: not at all!
  "cadd_15_ratio": {
    "r": {
      "HP:0007754": 0.8851477859029983
    },
    "d": {
      "HP:0007754": 0.9594599287088905
    }
  },
  # MOI_score > 0: recessive; MOI_score < 0: dominant
  "MOI_score": 0.8872573775123538,

  # pop_curse_flags look at per moi per HPO level.
  # it reports if certain POP specific variants are predominantly enriched in a group with moi/HPO
  "pop_curse_flags": {
    "r": {
      # e.g. we found patients with recessive/HP:0004329 tend to carry variants v1,v2,v3,v4
      # and v1,v2,v3,v4 seem to have a Jewish descent (inferred from gnomad), it then raises a 
      # ASJ flag.
      "HP:0004329": [
        "ASJ"
      ],
      "HP:0012145": [
        "AFR"
      ],
      "HP:0001098": [
        "ASJ"
      ],
      "HP:0000478": [
        "ASJ"
      ],
      "HP:0012372": [
        "ASJ"
      ],
      "HP:0000479": [
        "ASJ"
      ],
      "HP:0012374": [
        "ASJ"
      ],
      "HP:0012373": [
        "NFE"
      ],
      "HP:0001871": [
        "AFR"
      ],
      "HP:0000556": [
        "ASJ"
      ],
      "HP:0005561": [
        "AFR"
      ]
    },
    "d": {
      "HP:0000504": [
        "SAS"
      ],
      "HP:0007754": [
        "NFE"
      ],
      "HP:0012103": [
        "OTH"
      ],
      "HP:0001103": [
        "NFE"
      ],
      "HP:0011017": [
        "OTH"
      ],
      "HP:0001574": [
        "EAS"
      ],
      "HP:0001939": [
        "OTH"
      ]
    }
  },
  # this is the goodness_of_fit score. It only reports HPOs that are most relevant.
  "hgf": {
    "r": {
      "HP:0007754": 10.293977466877369
    },
    "d": {
      "HP:0007754": 12.581133080136787
    }
  },
  # number of patients found to carry at least two variants ("r") and at least one variant ("d")
  "NP": {
    "r": 51,
    "d": 187
  }
}
```

##### Another example of using Phenogenon:
```bash
phenogenon  --range 2:166845571-166930215  --vcf_file tests/data/SCN1A.anonymised.vcf.gz --output SCN1A.test.json
```
Explain output:
```python
{
  "cadd_15_ratio": {
    "r": {
      "HP:0012103": 0,
      "HP:0001250": 1
    },
    "d": {
      "HP:0001250": 0.9860913529587728
    }
  },
  # dominant moi, MOI_score < 0
  "MOI_score": -67.76587193983218,
  "pop_curse_flags": {
    "r": {},
    "d": {}
  },
  "hgf": {
    "r": {
      "HP:0012103": 1.2849541828056956,
      "HP:0001250": 2.514407022516984
    },
    "d": {
      # seizures
      "HP:0001250": 74.66686039643339
    }
  },
  "NP": {
    "r": 11,
    "d": 108
  }
}
```
