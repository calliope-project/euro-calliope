"""Applies config parameters to template files."""
from pathlib import Path

import jinja2

from eurocalliopelib import filters


def parametrise_template(path_to_template, path_to_output_yaml, **kwargs):
    """Applies config parameters to template files."""

    kwargs = _update_kwargs(**kwargs)
    path_to_template = Path(path_to_template)
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(path_to_template.parent),
        lstrip_blocks=True, trim_blocks=True, keep_trailing_newline=True,
        undefined=jinja2.StrictUndefined  # This ensures that missing pandas index elements raise an exception instead of silently returning None
    )
    env.filters['unit'] = filters.unit
    rendered =env.get_template(path_to_template.name).render(**kwargs)

    with open(path_to_output_yaml, "w") as result_file:
        result_file.write(rendered)


def _update_kwargs(**kwargs):
    if "scaling_factors" in kwargs.keys():
        kwargs["scaling_factors"]["specific_costs"] = (
            kwargs["scaling_factors"]["monetary"] / kwargs["scaling_factors"]["power"]
        )
    if "locations" in kwargs.keys():
        kwargs["locations"] = kwargs["locations"].rename(index=lambda x: x.replace(".", "-"))

    return kwargs
