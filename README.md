# Euro-Calliope

Models of the European electricity system for _Calliope_.

This repository contains the workflow routines that automatically build models from raw data. Alternatively to building models yourself, you can use [pre-built models](https://doi.org/10.5281/zenodo.3949553) that run out-of-the-box. You can find a more detailed description of the first application in a [scientific article in Joule](https://doi.org/10.1016/j.joule.2020.07.018).

[![article DOI](https://img.shields.io/badge/article-10.1016/j.joule.2020.07.018-blue)](https://doi.org/10.1016/j.joule.2020.07.018)
[![pre-built models DOI](https://img.shields.io/badge/prebuilts-10.5281%2Fzenodo.3949553-blue)](https://doi.org/10.5281/zenodo.3949553)
[![Documentation Status](https://readthedocs.org/projects/euro-calliope/badge/?version=latest)](https://euro-calliope.readthedocs.io/en/latest/?badge=latest)
[![Check Markdown links](https://github.com/calliope-project/euro-calliope/actions/workflows/externallinks.yaml/badge.svg)](https://github.com/calliope-project/euro-calliope/actions/workflows/externallinks.yaml)
[![Tests of YAML configuration and schema](https://github.com/calliope-project/euro-calliope/actions/workflows/schemavalidation.yaml/badge.svg)](https://github.com/calliope-project/euro-calliope/actions/workflows/schemavalidation.yaml)
[![Tests of eurocalliopelib and scripts](https://github.com/calliope-project/euro-calliope/actions/workflows/pythonpackage.yaml/badge.svg)](https://github.com/calliope-project/euro-calliope/actions/workflows/pythonpackage.yaml)

## At a glance

Euro-Calliope models the European electricity system with each location representing an administrative unit. It is built on three spatial resolutions: on the continental level as a single location, on the national level with 34 locations, and on the regional level with 497 locations. At each location, renewable generation capacities (wind, solar, bioenergy) and balancing capacities (battery, hydrogen) can be built. In addition, hydro electricity and pumped hydro storage capacities can be built up to the extent to which they exist today. All capacities are used to satisfy electricity demand on all locations where demand is based on historic data. Locations are connected through transmission lines of either unrestricted capacity or projections. Using [Calliope](https://www.callio.pe), the model is formulated as a linear optimisation problem with total monetary cost of all capacities as the minimisation objective. Due to the flexibility of Calliope and the availability of the routines building the model all components can be adapted to the modeller's needs.

## More information

Here is where you can find more information:

* [Full documentation](https://euro-calliope.readthedocs.io/)
* [Release notes](./CHANGELOG.md)
* [Citation information](./docs/about/citation.md)
* [License](./LICENSE.md)

If you are unable to access the full documentation via ReadTheDocs following the link above, or otherwise want to build the documentation locally, you can run the following from the repository top-level directory (assuming you have [conda](https://conda.io) installed):


```bash
conda env install -f envs/docs.yaml
conda activate docs
mkdocs build --no-directory-urls
```

The documentation can then be accessed by opening `build/docs/index.html`.

## License

Euro-Calliope is developed and maintained within the [Calliope project](https://www.callio.pe). The code in this repository is [MIT licensed](./LICENSE.md).
