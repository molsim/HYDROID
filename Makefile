
clean-build:
	rm --force --recursive build/
	rm --force --recursive dist/
	rm --force --recursive *.egg-info

sdist:
	rm -rf dist
	python setup.py sdist

pypi-push:
	twine upload --repository-url https://upload.pypi.org/legacy/ dist/*