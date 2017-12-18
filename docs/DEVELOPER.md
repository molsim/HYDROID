# HYDROID, Developer Docs

## Deploying HYDROID to PyPi
Increment hydroid/VERSION

`make sdist`
`make pypi-push`

.pipyrc should look like this
`
index-servers =
  pypi

[pypi]
repository=https://upload.pypi.org/legacy/
username=....
password=....
`

## Building anaconda packages

`conda skeleton pypi hydroid`

edit meta.yml
comment out `tests: commands:` section
build.sh change to
`$PYTHON setup.py install --single-version-externally-managed --record=record.txt`
bld.bat to
`"%PYTHON%" setup.py install --single-version-externally-managed --record=record.txt`

`conda-build -c conda-forge hydroid`
anaconda upload __PATH_to_Module__

