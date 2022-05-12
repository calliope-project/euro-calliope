import filecmp

import pandas as pd
import pytest
import jinja2

from eurocalliopelib.template import parametrise_template

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
    bar: {{ param1 }}
    baz: {{ param2 }}
"""

TEMPLATE_PARAMS_EXPECTED = """
foo:
    bar: {param1}
    baz: {param2}
"""

TEMPLATE_SCALING_FACTOR_IN_PARAMS = """
foo:
    bar: {{ scaling_factors.power }}
    baz: {{ scaling_factors.specific_costs }}
"""

TEMPLATE_SCALING_FACTOR_IN_PARAMS_EXPECTED = """
foo:
    bar: 0.5
    baz: 4.0
"""

TEMPLATE_LOCATIONS_IN_PARAMS = """
foo:
    bar: {{ locations.loc['A-B', 'foo'] }}
"""

TEMPLATE_LINKS_IN_PARAMS = """
foo:
    bar: {{ links.loc['A-B,C-D', 'foo'] }}
"""

TEMPLATE_LOCATIONS_LINKS_IN_PARAMS_EXPECTED = """
foo:
    bar: 1.0
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
TEMPLATE_LOCATIONS_ITERATE_EXPECTED = """
foo:
    A-B: 1.0
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

@pytest.fixture(scope="module")
def links():
    return pd.DataFrame({
        "foo": {"A.B,C.D": 1.0}
    })


def test_template_no_params(template_to_file_obj, out_file_obj):
    template_obj = template_to_file_obj(TEMPLATE_NO_PARAMS)
    parametrise_template(template_obj, str(out_file_obj))

    filecmp.cmp(template_obj, out_file_obj)


def test_template_line_strip(template_to_file_obj, out_file_obj):
    """
    Jinja args are set to trim the template to remove empty lines where there are
    jinja2 commands (e.g. a for loop).
    """
    not_expected = template_to_file_obj(TEMPLATE_LINE_TO_STRIP)
    expected = template_to_file_obj(TEMPLATE_NO_PARAMS)
    parametrise_template(not_expected, str(out_file_obj))

    not filecmp.cmp(not_expected, out_file_obj)
    filecmp.cmp(expected, out_file_obj)


@pytest.mark.parametrize(("param1", "param2"), [(1, 2), ("a", "b")])
def test_template_params(template_to_file_obj, out_file_obj, param1, param2):
    template_obj = template_to_file_obj(TEMPLATE_PARAMS)
    expected = TEMPLATE_PARAMS_EXPECTED.format(param1=param1, param2=param2)

    parametrise_template(
        template_obj, out_file_obj,
        param1=param1, param2=param2
    )
    assert expected == out_file_obj.read()


def test_template_scaling_factors(template_to_file_obj, out_file_obj, scaling_factors):
    template_obj = template_to_file_obj(TEMPLATE_SCALING_FACTOR_IN_PARAMS)
    expected = TEMPLATE_SCALING_FACTOR_IN_PARAMS_EXPECTED

    parametrise_template(
        template_obj, out_file_obj,
        scaling_factors=scaling_factors
    )

    assert expected == out_file_obj.read()


def test_template_locations(template_to_file_obj, out_file_obj, locations):
    template_obj = template_to_file_obj(TEMPLATE_LOCATIONS_IN_PARAMS)
    expected = TEMPLATE_LOCATIONS_LINKS_IN_PARAMS_EXPECTED

    parametrise_template(
        template_obj, out_file_obj,
        locations=locations
    )

    assert expected == out_file_obj.read()


def test_template_links(template_to_file_obj, out_file_obj, links):
    template_obj = template_to_file_obj(TEMPLATE_LINKS_IN_PARAMS)
    expected = TEMPLATE_LOCATIONS_LINKS_IN_PARAMS_EXPECTED

    parametrise_template(
        template_obj, out_file_obj,
        links=links
    )

    assert expected == out_file_obj.read()


def test_template_wrong_locations(template_to_file_obj, out_file_obj, locations):
    template_obj = template_to_file_obj(TEMPLATE_WRONG_LOCATIONS_IN_PARAMS)
    with pytest.raises(jinja2.exceptions.UndefinedError, match=r".*('A.B', 'foo')"):
        parametrise_template(
            template_obj, out_file_obj,
            locations=locations
        )


def test_template_locations_iterate(template_to_file_obj, out_file_obj, locations):
    template_obj = template_to_file_obj(TEMPLATE_LOCATIONS_ITERATE)
    expected = TEMPLATE_LOCATIONS_ITERATE_EXPECTED

    parametrise_template(
        template_obj, out_file_obj,
        locations=locations
    )

    assert expected == out_file_obj.read()
