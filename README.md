## AGEpy

This python package contains Bioinformatics tools developed at the
Bioinformatics Core Facility of the Max Planck Institute for Biology of Ageing.

> Max Planck Institute for Biology of Ageing  
> Joseph-Stelzmann-Str. 9b  
> D-50931 Cologne  
> Germany

#### Read the Docs

We have started introducing some documentation [here](http://agepy.readthedocs.org).


#### Installation

##### Dependencies

AGEpy requires [R](https://www.r-project.org) and the [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) package for R.

For installing R follow the instructions [here](https://www.r-project.org).

Once you have installed R you are ready to install biomaRt

```R

$ R

> source("http://bioconductor.org/biocLite.R") 

> biocLite()

> biocLite("biomaRt") 

> quit()

```

##### AGEpy

###### github

Get the latest release from github:

```bash
git clone https://github.com/mpg-age-bioinformatics/AGEpy
```

Install:

````bash
cd AGEpy
python setup.py install --user
```

and then update to the latest release whenever required with:

```bash
cd AGEpy
git pull
python setup.py install --user --force

```

Alternatively you can also install the package with a symlink, so that changes
to the source files will be immediately available to users of the package on
your system:

````bash
cd AGEpy
python setup.py develop --user
```

Be aware that his the develop option you won't be able to properly update once new scripts are added.

###### pip

If you have pip installed you can get the current pip version over pip:

```bash
pip install --user AGEpy
```

and upgrade whenever required with:

```bash
pip install --user AGEpy --upgrade
```

#### Help

In bash:

```bash
pydoc AGEpy.AGEpy
```

In python:

```python
help("AGEpy.AGEpy")
```

#### Example usage

```python
import AGEpy.AGEpy as age

gtf=age.readGTF("/path/to/file.gtf")

gtf.head()
```

#### Scripts

* `david` a script to perform enrichment analysis from the DAVID database.
The usage is described in the script's help output called via `david --help`.
More information at: https://github.com/mpg-age-bioinformatics/AGEpy/wiki/david

* `bit` The [b]ermuda [i]nformation [t]riangle is a git-based tool for the management of code and data.
Check out https://github.com/mpg-age-bioinformatics/AGEpy/wiki/bit
