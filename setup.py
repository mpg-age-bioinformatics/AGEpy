from setuptools import setup

setup(name = 'AGEpy',
      version = '0.5.0',
      description = 'Bioinformatics tools for Python developed at the MPI for Biology of Ageing',
      url = 'https://github.com/mpg-age-bioinformatics/AGEpy',
      author = 'Bioinformatics Core Facility of the Max Planck Institute for Biology of Ageing',
      author_email = 'bioinformatics@age.mpg.de',
      license = 'MIT',
      packages = [ 'AGEpy' ],
      install_requires = [ 'Pandas>=0.15.2', 'numpy>=1.9.2','requests==2.10.0', \
      'suds', 'xlrd', 'biomart', 'rpy2', 'matplotlib', \
      'xlsxwriter','pybedtools'],
      zip_safe = False,
      scripts=['bin/david','bin/obo2tsv']
      )

#scripts=['bin/david','bin/bit','bin/obo2tsv']  
