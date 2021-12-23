# eurocalliopelib

Library code of the Euro-Calliope workflow.

The library code contains general-purpose functions and routines that we expect to change at a slow pace -- in contrast to scripts. If you alter library code be aware that this will not trigger reruns of workflow rules. Think of `eurocalliopelib` as any other dependency of this workflow (like NumPy): changing the version of any dependency will not rerun worfklow rules. When you change library code, you will have to rerun rules manually where needed.

## Developer Guide

The following assumes you are in the root folder of this repository, i.e. the parent folder of the folder this file is in.

### Installation

Best install `eurocalliopelib` from the conda environment:

    $ conda env create -f test-requirements.yaml
    $ conda activate test-eurocalliope

### Run the test suite

Run the test suite with py.test:

    $ py.test tests/lib
