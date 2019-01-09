from glob import glob

from setuptools import setup, find_packages

setup(name='phenogenon',
      version='0.1.2',
      description='Phenogenon is a method that combines: the power of Human Phenotype Ontology for describing patient phenotypes, gnomAD for estimating rare variant population frequency, and CADD for variant pathogenicity prediction.',
      url='https://github.com/phenopolis/phenogenon',
      author='Jing Yu',
      author_email='info@phenopolis.org',
      license='MIT',
      packages=find_packages(),
      scripts=glob('./scripts/*'),
      install_requires=['pymongo==3.4.0',
                        'biopython==1.68',
                        'pyliftover==0.3',
                        'bs4==0.0.1',
                        'scipy==0.19.0',
                        'pandas==0.20.2',
                        'pysam==0.11.2.2',
                        'fisher==0.1.5',]
)
