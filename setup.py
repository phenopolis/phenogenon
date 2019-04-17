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
      install_requires=['fisher==0.1.5',
                        'numpy==1.16.0',
                        'pandas==0.23.4',
                        'pysam==0.15.2',
                        'scipy==1.1.0',
                        'plotly==3.8.0',
                        'psutil==5.6.1']
      )
