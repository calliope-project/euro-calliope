# About

This module is dedicated to getting and processing energy demand for industry.

Basic info of this module:

- Main data sources: JRC IDEES and eurostat
- Spatial resolution: national (because of JRC IDEES)
- Temporal resolution: annual aggregated

Some industries are specified as 'sub-module' in this module, such as iron and steel industry. This allows you to flexibly configure the decarbonisation level of each industry separately. All other industries are grouped into 'other industries'.

For workflow users, i.e. non-developers, `config.yaml` is the main file to look into.
Here, one can specify parameters such as the years to be included in the data processing, assumptions such as the share of recycled steel, and the specific industry sector to be included (such as iron and steel industry).
The structure of the `config.yaml` file is:

- inputs - the needed data sources from other modules in euro-calliope (for now, we assume this module depends on other modules and is not yet a stand-alone module)
- outputs - the output of the whole industry module, usually files, possibly passed on to other modules in euro-calliope.
- params - the parameters that affect the calculation process and result in this module.
By changing the value of the parameters, each user can tailor the workflow to their own needs.
- setup - anything that concerns the general data pipeline of the module.
