import filecmp

import pandas as pd
import pytest
import jinja2

from eurocalliopelib import parametrise_template

TEMPLATE_NO_PARAMS = """
foo:
    bar: 1
"""

TEMPLATE_LINE_TO_STRIP = """
foo:
    {% for i in [1] %}
    bar: 1
    {% endfor %}
"""

TEMPLATE_PARAMS = """
foo:
    bar: {{ foobar }}
    baz: {{ foobaz }}
"""

TEMPLATE_SCALING_FACTOR_IN_PARAMS = """
foo:
    bar: {{ scaling_factors.power }}
    baz: {{ scaling_factors.specific_costs }}
"""

TEMPLATE_LOCATIONS_IN_PARAMS = """
foo:
    bar: {{ locations.loc['A-B', 'foo'] }}
"""

TEMPLATE_WRONG_LOCATIONS_IN_PARAMS = """
foo:
    bar: {{ locations.loc['A.B', 'foo'] }}
"""

TEMPLATE_LOCATIONS_ITERATE = """
foo:
    {% for idx, data in locations.iterrows() %}
    {{ idx }}: {{ data.foo }}
    {% endfor %}
"""


@pytest.fixture(scope="module")
def template_to_file_obj(tmpdir_factory):
    def _template_to_file_obj(template):
        path_obj = tmpdir_factory.mktemp("model").join("no_params.yaml")
        path_obj.write(template)
        return path_obj
    return _template_to_file_obj


@pytest.fixture(scope="module")
def out_file_obj(tmpdir_factory):
    path_obj = tmpdir_factory.mktemp("model").join("out.yaml")
    return path_obj


@pytest.fixture(scope="module")
def scaling_factors():
    return {
        "power": 0.5,
        "monetary": 2
    }


@pytest.fixture(scope="module")
def locations():
    return pd.DataFrame({
        "foo": {"A.B": 1.0}
    })


def format_template_manually(template, mapping):
    for template_val, replacement_val in mapping.items():
        template = template.replace("{{ " + template_val + " }}", str(replacement_val))
    return template


def test_template_no_params(template_to_file_obj, out_file_obj):
    template_obj = template_to_file_obj(TEMPLATE_NO_PARAMS)
    parametrise_template.parametrise_template(template_obj, str(out_file_obj))

    assert filecmp.cmp(template_obj, out_file_obj)


def test_template_line_strip(template_to_file_obj, out_file_obj):
    """
    Jinja args are set to trim the template to remove empty lines where there are
    jinja2 commands (e.g. a for loop).
    """
    template_obj_extra_line = template_to_file_obj(TEMPLATE_LINE_TO_STRIP)
    template_obj_no_extra_line = template_to_file_obj(TEMPLATE_NO_PARAMS)
    parametrise_template.parametrise_template(template_obj_extra_line, str(out_file_obj))

    assert not filecmp.cmp(template_obj_extra_line, out_file_obj)
    assert filecmp.cmp(template_obj_no_extra_line, out_file_obj)


@pytest.mark.parametrize(("foobar", "foobaz"), [(1, 2), ("a", "b")])
def test_template_params(template_to_file_obj, out_file_obj, foobar, foobaz):
    template_obj = template_to_file_obj(TEMPLATE_PARAMS)

    parametrise_template.parametrise_template(
        template_obj, out_file_obj,
        foobar=foobar, foobaz=foobaz
    )
    expected_result = format_template_manually(
        TEMPLATE_PARAMS, {"foobar": foobar, "foobaz": foobaz}
    )
    assert expected_result == out_file_obj.read()


def test_template_scaling_factors(template_to_file_obj, out_file_obj, scaling_factors):
    template_obj = template_to_file_obj(TEMPLATE_SCALING_FACTOR_IN_PARAMS)

    parametrise_template.parametrise_template(
        template_obj, out_file_obj,
        scaling_factors=scaling_factors
    )
    expected_result = format_template_manually(
        TEMPLATE_SCALING_FACTOR_IN_PARAMS,
        {"scaling_factors.power": 0.5, "scaling_factors.specific_costs": 4.0}
    )
    assert expected_result == out_file_obj.read()


def test_template_locations(template_to_file_obj, out_file_obj, locations):
    template_obj = template_to_file_obj(TEMPLATE_LOCATIONS_IN_PARAMS)

    parametrise_template.parametrise_template(
        template_obj, out_file_obj,
        locations=locations
    )
    expected_result = format_template_manually(
        TEMPLATE_LOCATIONS_IN_PARAMS,
        {"locations.loc['A-B', 'foo']": 1.0}
    )
    assert expected_result == out_file_obj.read()


def test_template_wrong_locations(template_to_file_obj, out_file_obj, locations):
    template_obj = template_to_file_obj(TEMPLATE_WRONG_LOCATIONS_IN_PARAMS)
    with pytest.raises(jinja2.exceptions.UndefinedError, match=r".*('A.B', 'foo')"):
        parametrise_template.parametrise_template(
            template_obj, out_file_obj,
            locations=locations
        )


def test_template_locations_iterate(template_to_file_obj, out_file_obj, locations):
    template_obj = template_to_file_obj(TEMPLATE_LOCATIONS_ITERATE)

    parametrise_template.parametrise_template(
        template_obj, out_file_obj,
        locations=locations
    )
    assert "A-B: 1.0" in out_file_obj.read()
