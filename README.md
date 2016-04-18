## AGEpy

This python package contains Bioinformatics tools developed at the
Bioinformatics Core Facility of the Max Planck Institute for Biology of Ageing.

> Max Planck Institute for Biology of Ageing  
> Joseph-Stelzmann-Str. 9b  
> D-50931 Cologne  
> Germany

We are currently on beta and will introduce the needed documention asap.

#### Installation

###### github

Get the latest release from github:

```bash
git clone https://github.com/mpg-age-bioinformatics/AGEpy
```

Install:

````bash
cd AGEpy
python setup.py install --user --install-scripts
```

Alternatively you can also install the package with a symlink, so that changes
to the source files will be immediately available to users of the package on
your system:

````bash
cd AGEpy
python setup.py develop --user --install-scripts
```

and then update to the latest release whenever required with:

```bash
cd AGEpy
git pull
```

and 

```bash
python setup.py develop --user --install-scripts
```

if new scripts have been added.


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


