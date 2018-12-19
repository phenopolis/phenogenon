# Phenogenon
[![Build Status](https://travis-ci.com/phenopolis/phenogenon.svg?branch=master)](https://travis-ci.com/phenopolis/phenogenon)

# Installation


## Standalone
External programs you need to install:
* tabix
* bcftools
* jq

```
pip install -r requirements.txt
```

### Docker

Make sure you have docker running.
The following will build and run docker image with all the dependencies.

Build docker image.
```
bash build-docker.sh
```
Run docker:
```
bash run-docker.sh
```

# Running Phenogenon

Run tests:
```
coverage run --omit=*/site-packages/*,*/tests/* -m unittest discover -s tests
```
Example:
```
python2 lib/goodness_of_fit.py  --range 1:94458394-94586689 --vcf_file tests/data/ABCA4.anonymised.vcf.gz --output ABCA4.test.json
```
Produces `ABCA4.test.json`.

Explanation of output:
```
{
  "genon_hratios": {
    "r": {
      "HP:0007754": 0.8046670602618868
    },
    "d": {
      "HP:0007754": 0.7483526819524269
    }
  },
  "genon_vratios": {
    "r": {
      "HP:0007754": 0.8851477859029983
    },
    "d": {
      "HP:0007754": 0.9594599287088905
    }
  },
  "predicted_moi": 0.8872573775123538,
  "pop_curse_flags": {
    "r": {
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
  "hgf": {
    "r": {
      "HP:0007754": 10.293977466877369
    },
    "d": {
      "HP:0007754": 12.581133080136787
    }
  },
  "NP": {
    "r": 51,
    "d": 187
  },
  "genon_sratios": {
    "r": {
      "HP:0007754": 0.6504054838971388
    },
    "d": {
      "HP:0007754": 0.4616437948047774
    }
  }
}
```

```
python2 lib/goodness_of_fit.py  --range 2:166845571-166930215  --vcf_file tests/data/SCN1A.anonymised.vcf.gz --output SCN1A.test.json
```
