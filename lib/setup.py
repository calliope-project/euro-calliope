#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='eurocalliopelib',
    version='1.2.0.dev', # additionally defined in __init__.py
    description='Library code of the Euro-Calliope workflow.',
    maintainer='calliope-project',
    maintainer_email='tim.troendle@usys.ethz.ch',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        # "numpy", # TODO readd after solving #262
        # "scipy", # TODO readd after solving #262
        # "pandas", # TODO readd after solving #262
        # "xarray", # TODO readd after solving #262
        # "pycountry==18.12.8" # TODO readd after solving #262
    ],
    extras_require={
        'geo': [
            # "geopandas", # TODO readd after solving #262
            # "rasterio", # TODO readd after solving #262
            # "rasterstats", # TODO readd after solving #262
        ],
        'docs': [
            # "pydot", # TODO readd after solving #262
            # "graphviz", # TODO readd after solving #262
            # "mkdocs", # TODO readd after solving #262
            # "jsonschema2md" # TODO readd after solving #262
        ]
    },
    entry_points={
        'mkdocs.plugins': [
            'dag = eurocalliopelib.docs.dag:DAGPlugin',
            'schema = eurocalliopelib.docs.schema:SchemaPlugin',
            'add-file = eurocalliopelib.docs.addfile:AddFilePlugin'
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
