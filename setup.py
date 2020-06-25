import codecs
import os
import re
import sys
from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    # intentionally *not* adding an encoding option to open, See:
    #   https://github.com/pypa/virtualenv/issues/201#issuecomment-3145690
    with codecs.open(os.path.join(here, *parts), 'r') as fp:
        return fp.read()

setup(name = 'AGEpy',
      version = '0.8.1',
      description = 'Bioinformatics tools for Python developed at the MPI for Biology of Ageing',
      long_description = read('README.rst'),
      url = 'https://github.com/mpg-age-bioinformatics/AGEpy',
      author = 'Bioinformatics Core Facility of the Max Planck Institute for Biology of Ageing',
      author_email = 'bioinformatics@age.mpg.de',
      license = 'MIT',
      packages = [ 'AGEpy' ],
      install_requires = [ 'Pandas>=0.15.2', 'numpy>=1.9.2','requests>=2.20.0', \
      'suds-jurko', 'xlrd', 'biomart', 'matplotlib', 'client', \
      'xlsxwriter','pybedtools','wand','paramiko','ipaddress', 'seaborn'],
      zip_safe = False,
      scripts=['bin/obo2tsv','bin/aDiff','bin/abed','bin/david', 'bin/blasto']
      )
