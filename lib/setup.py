#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='eurocalliopelib',
    version='1.1.0.dev', # additionally defined in __init__.py
    description='Library code of the Euro-Calliope workflow.',
    maintainer='calliope-project',
    maintainer_email='tim.troendle@usys.ethz.ch',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "xarray",
        "pycountry==18.12.8"
    ],
    extras_require={
        'geo': [
            "geopandas",
            "rasterio",
            "rasterstats",
        ],
        'docs': [
            "pydot",
            "graphviz",
            "mkdocs",
            "jsonschema2md"
        ]
    },
    entry_points={
        'mkdocs.plugins': [
            'dag = eurocalliopelib.docs.dag:DAGPlugin',
            'schema = eurocalliopelib.docs.schema:SchemaPlugin',
            'softlink-file = eurocalliopelib.docs.softlink:SoftLinkPlugin'
        ]
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
