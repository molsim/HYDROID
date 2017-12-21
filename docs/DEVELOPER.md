# HYDROID, Developer Docs

## What to do to ship a release:
- increment version file hydroid/VERSION
- tag release
- add example files to release
- push to PiPy
- make new anaconda package and deploy via this recipe https://github.com/intbio/hydroid-conda


## Deploying HYDROID to PyPi

```
make sdist
make pypi-push
```

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

```
conda install conda-build
conda update -n root conda-build
conda install anaconda-client
anaconda login
```


