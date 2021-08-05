#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='eurocalliopelib',
    version='1.1.2.dev', # additionally defined in __init__.py
    description='Library code of the euro-calliope workflow.',
    maintainer='calliope-project',
    maintainer_email='tim.troendle@usys.ethz.ch',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "xarray",
        "pycountry>=18.12.8"
    ],
    extras_require={
        'geo': [
            "geopandas",
            "rasterio",
            "rasterstats",
        ],
    },
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering'
    ]
)
