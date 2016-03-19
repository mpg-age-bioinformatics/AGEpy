from setuptools import setup

setup(name = 'AGEpy',
      version = '0.1.5',
      description = 'Bioinformatics tools for Python developed at the MPI for Biology of Ageing',
      url = 'https://github.com/mpg-age-bioinformatics/AGEpy',
      author = 'Jorge Boucas, Sven E. Templer',
      author_email = 'jorge.boucas@age.mpg.de',
      license = 'MIT',
      packages = [ 'AGEpy' ],
      install_requires = [ 'Pandas>=0.15.2', 'numpy>=1.9.2', 'suds', 'xlrd','biomart','rpy2','matplotlib' ],
      zip_safe = False
      )
