# HYDROID, Developer Docs

## What to do to ship a release:
- increment version file hydroid/VERSION
- tag release
- add example files to release
- make new anaconda package and deploy via this recipe https://github.com/intbio/hydroid-conda


### Optionally deploy HYDROID to PyPi (no FreeSASA module)

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


