#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='eurocalliopelib',
    version='0.1.0', # additionally defined in __init__.py
    description='Library code of the euro-calliope workflow.',
    maintainer='calliope-project',
    maintainer_email='tim.troendle@usys.ethz.ch',
    packages=find_packages(exclude=['tests*']),
    include_package_data=True,
    install_requires=["pycountry>=18.12.8"],
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering'
    ]
)
