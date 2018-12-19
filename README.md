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
python2 lib/goodness_of_fit.py --range 1:94458394-94586689 --vcf_file chr1.vcf.gz
```
