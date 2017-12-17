from setuptools import setup, find_packages
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))

version_file = open(os.path.join(here, 'HYDROID', 'VERSION'))
version = version_file.read().strip()

# Get the long description from the relevant file
with codecs.open(os.path.join(here, 'DESCRIPTION.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='HYDROID',

    # Versions should comply with PEP440. For single-sourced versioning, see
    # http://packaging.python.org/en/latest/tutorial.html#version
    version=version,

    description='HYDROID (HYDroxyl-Radical fOotprinting Interpretation for DNA) is a python package for the analysis of experimental data generated by hydroxyl-radical footprinting (HRF) of DNA-protein complexes',
    long_description=long_description,

    # The project URL.
    url='https://github.com/ncbi/HYDROID',

    # Author details
    author='Alexey K. Shaytan, Grigoriy A. Armeev, Anna R. Panchenko',
    author_email='alex@intbio.org',

    # Choose your license
    license='Public Domain',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # https://pypi.python.org/pypi?%3Aaction=list_classifiers
        # Project maturity. 
        'Development Status :: 3 - Alpha',

        # Intended audience
        'Intended Audience :: Science/Research',

        # Topic
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # License should match "license" above
        'License :: Public Domain',

        # Python versions supported
        'Programming Language :: Python :: 2.7',
    ],

    # What does your project relate to?
    keywords='science hydroxyl-radical footprinting image analysis',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['examples', 'docs', 'tests*']),

    # Run-time package dependencies. These will be installed by pip when your
    # project is installed.
    install_requires=[
        'biopython==1.68',
        'cycler==0.10.0',
        'functools32==3.2.3.post2',
        'matplotlib==2.0.0',
        'numpy==1.12.0',
        'olefile==0.44',
        'packaging==16.8',
        'pandas==0.19.2',
        'patsy==0.4.1',
        'PeakUtils==1.0.3',
        'Pillow==4.0.0',
        'pyparsing==2.1.10',
        'python-dateutil==2.6.0',
        'pytz==2016.10',
        'scipy==0.18.1',
        'six==1.10.0',
        'statsmodels==0.8.0',
        'subprocess32==3.2.7',
    ],
    extras_require = {
    'pred': [
        'Cython',
    ],},

    # Data files included in your packages. If using Python 2.6 or less, 
    # then these have to be included in MANIFEST.in as well.
    package_data={
        'HYDROID': ['pkgdata/amber10_rmin.config','pkgdata/charmm36_rmin.config','pkgdata/cnr.otf','pkgdata/cnrb.otf'],
    },

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    #entry_points={
    #    'console_scripts': [
    #        'HYDROID=HYDROID.main:main',
    #    ],
    #},
    python_requires='==2.7.*',
    # MANIFEST.in included entries should be included as package data and
    # installed into site-packages 
    include_package_data=True,

    # Default to False unless you specifically intend the package to be
    # installed as an .egg
    zip_safe=False,
)