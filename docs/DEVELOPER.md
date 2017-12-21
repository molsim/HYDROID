# HYDROID, Developer Docs

## What to do to ship a release:
- increment version file
- tag release
- add example files to release
- push to PiPy
- make new anaconda package


## Deploying HYDROID to PyPi
Increment hydroid/VERSION

```make sdist
make pypi-push```

.pipyrc should look like this
```
index-servers =
  pypi

[pypi]
repository=https://upload.pypi.org/legacy/
username=....
password=....
```

## Building anaconda packages

```conda install conda-build
conda update -n root conda-build
conda install anaconda-client
anaconda login```

```conda skeleton pypi hydroid```

edit meta.yml
comment out `tests: commands:` section
build.sh change to
`$PYTHON setup.py install --single-version-externally-managed --record=record.txt`
bld.bat to
`"%PYTHON%" setup.py install --single-version-externally-managed --record=record.txt`


```conda config --set anaconda_upload yes
conda-build -c conda-forge hydroid
conda convert --platform all /Users/alexsha/miniconda2/conda-bld/osx-64/hydroid-0.0.4.post10-np111py27_0.tar.bz2 -o output/
find output/ -name 'hydroid*' -exec anaconda upload {} \;```

