### Importing

All functions in the AGEpy pakcage can be accessed using:

```python
import AGEpy as age
help(age.readGTF)
```

Alternatively, functions from the different modules can be accessed with for example:

```python
from AGEpy import gtf
help(gtf.readGTF)
```

### Help

In bash:

```bash
pydoc AGEpy.AGEpy
```

In python:

```python
help("AGEpy.AGEpy")
```

### Example usage

```python
import AGEpy as age

gtf=age.readGTF("/path/to/file.gtf")

gtf.head()
```
