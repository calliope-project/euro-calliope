# Customising the workflow

Beyond [customisation options of the models](../model/customisation.md), the workflow gives you full control configuring, adapting, and extending models. There are two options:

1. [Configuration](./customisation.md#configuration)
2. [Adaptation](./customisation.md#adaptation)

## Configuration

You can configure many aspects in the build process by changing configuration parameter values.
You can, for example, change the temporal and spatial scope of the data and model, change data sources, and change the way raw data is preprocessed.
[We list and describe all configuration parameters in this documentation.](./schema.md)

The configuration builds on Snakemake's configuration mechanism and consists of two parts: a default configuration `./config/default.yaml` and a schema declaring all configuration parameters `./config/schema.yaml`.
To override configuration parameters, you can add another configuration file with just your updates or change parameter values on the command line when calling `snakemake`.
For details on how the configuration mechanism works, please read [Snakemake's documention](https://snakemake.readthedocs.io/en/v8.10.7/snakefiles/configuration.html).

## Adaptation

Beyond configuration through parameters, you can adapt and extend the workflow in any possible way.
You can adapt the data pre-processing steps and the way model files are generated, but you can also extend the model by adding your own model files or overrides.
Customising Euro-Calliope in this way requires a solid understanding of the workflow management system [Snakemake](https://snakemake.readthedocs.io/en/v8.10.7/index.html) that we use.

Whenever we applied Euro-Calliope in our research we made use of this option.
Below you will find a list of publications in which we applied Euro-Calliope models.
These may serve as a starting point for you to understand the possibilities of the adaptable and extendable workflow.

* Tröndle, T., Lilliestam, J., Marelli, S., &#38; Pfenninger, S. (2020). Trade-offs between geographic scale, cost, and infrastructure requirements for fully renewable electricity in Europe. <i>Joule</i>, <i>4</i>(9), 1929–1948. [![article DOI](https://img.shields.io/badge/article-10.1016/j.joule.2020.07.018-blue)](https://doi.org/10.1016/j.joule.2020.07.018)[![workflow DOI](https://img.shields.io/badge/workflow-10.5281/zenodo.3950774-blue)](https://doi.org/10.5281/zenodo.3950774)
* Tröndle, T. (2020). Supply-side options to reduce land requirements of fully renewable electricity in Europe. <i>PLOS ONE</i>. [![article DOI](https://img.shields.io/badge/article-10.1371/journal.pone.0236958-blue)](https://doi.org/10.1371/journal.pone.0236958)[![workflow DOI](https://img.shields.io/badge/workflow-10.5281/zenodo.3956530-blue)](https://doi.org/10.5281/zenodo.3956530)
* "Open Source Energiewende" multi-model analysis (2019) [![model DOI](https://img.shields.io/badge/model-10.5281/zenodo.4085047-blue)](https://doi.org/10.5281/zenodo.4085047)[![workflow DOI](https://img.shields.io/badge/workflow-github/ose-blue)](https://github.com/timtroendle/calliope-in-ose-model-comparison)
