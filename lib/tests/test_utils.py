import pytest
import io

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
from eurocalliopelib.utils import convert_country_code, read_eurostat_tsv, remove_digits, add_idx_level, to_numeric


class TestEurostatTSV:
    FILE_NUMERIC = """
nrg_bal,siec,unit,geo\\time	2019	2018	2017
AFC,BIOE,GWH,AL	3211.056	3129	1.0
AFC,BIOE,TW,AL	3211.056	3129	0
    """
    FILE_WITH_STRINGS = """
nrg_bal,siec,unit,geo\\time	2019	2018	2017
AFC,BIOE,GWH,AL	3211.056	3129.617c	-
AFC,BIOE,TW,AL	3211.056	3129	0
    """
    INDEX_LVLS = ["A", "B", "C", "D"]

    @pytest.fixture
    def df(self):
        def _get_df(file):
            filepath = io.StringIO(file)
            return read_eurostat_tsv(filepath, self.INDEX_LVLS)
        return _get_df

    @pytest.mark.parametrize("file", [FILE_NUMERIC, FILE_WITH_STRINGS])
    def test_eurostat_tsv_index_levels(self, df, file):
        assert df(file).index.names == self.INDEX_LVLS

    @pytest.mark.parametrize("file", [FILE_NUMERIC, FILE_WITH_STRINGS])
    def test_eurostat_tsv_index_values(self, df, file):
        assert df(file).index.difference([
            ("AFC", "BIOE", "GWH", "AL"), ("AFC", "BIOE", "TW", "AL")
        ]).empty

    @pytest.mark.parametrize("file", [FILE_NUMERIC, FILE_WITH_STRINGS])
    def test_eurostat_tsv_column_name(self, df, file):
        assert df(file).columns.name == "year"

    @pytest.mark.parametrize("file", [FILE_NUMERIC, FILE_WITH_STRINGS])
    def test_eurostat_tsv_column_values(self, df, file):
        assert set(df(file).columns) == set([2019, 2018, 2017])

    @pytest.mark.parametrize("file", [FILE_NUMERIC, FILE_WITH_STRINGS])
    def test_eurostat_tsv_column_dtype(self, df, file):
        assert df(file).columns.dtype == int

    @pytest.mark.parametrize("file", [FILE_NUMERIC, FILE_WITH_STRINGS])
    def test_eurostat_tsv_dtype(self, df, file):
        assert df(file).stack().dtype.kind in ["i", "f"]

    @pytest.mark.parametrize("file", [FILE_NUMERIC, FILE_WITH_STRINGS])
    def test_eurostat_tsv_slice_index(self, file):
        filepath = io.StringIO(file)
        df = read_eurostat_tsv(
            filepath, self.INDEX_LVLS, slice_lvl="C", slice_idx="TW"
        )
        assert df.index.difference([("AFC", "BIOE", "AL")]).empty
        assert np.allclose(df.values, [[3211.056, 3129, 0]])


class TestAddIndexLevel:
    DF = pd.DataFrame([["foo", "bar"]], columns=["A", "B"], index=pd.Index([1], name="C"))
    SERIES = pd.Series(["foo"], name="A", index=pd.Index([1], name="C"))

    @pytest.mark.parametrize(
        ("level_name", "level_value"), [("foobar", "baz"), ("foobar", 1)]
    )
    @pytest.mark.parametrize("data", [SERIES, DF])
    def test_add_to_series(self, data, level_name, level_value):
        new_data = add_idx_level(data, **{level_name: level_value})
        assert set(new_data.index.names) == set(["C", level_name])
        assert new_data.index.get_level_values(level_name) == [level_value]
        assert new_data.droplevel(level_name).equals(data)

    def test_setting_level_same_as_series_name(self):
        with pytest.raises(KeyError) as excinfo:
            add_idx_level(self.SERIES, A="bar")

@pytest.mark.parametrize(
    "input_country,iso3166,output",
    [
        ("uk", "GBR", "alpha3"),
        ("el", "GRC", "alpha3"),
        ("de", "DEU", "alpha3"),
        ("Germany", "DE", "alpha2"),
        ("United Kingdom", "GB", "alpha2"),
        ("United Kingdom", "UK", "alpha2_eurostat"),
        ("Greece", "GR", "alpha2"),
        ("Greece", "EL", "alpha2_eurostat"),
        ("bh", "BA", "alpha2"),
    ],
)
def test_country_code_conversion(input_country, iso3166, output):
    assert convert_country_code(input_country, output) == iso3166


@pytest.mark.parametrize(
    "input_country,iso3166",
    [
        ("uk", "GBR"),
        ("el", "GRC"),
        ("de", "DEU"),
        ("Germany", "DEU"),
        ("United Kingdom", "GBR"),
    ],
)
def test_country_code_conversion_no_output_specified(input_country, iso3166):
    assert convert_country_code(input_country) == iso3166


@pytest.mark.parametrize("string_with_digits", ["trial01"])
def test_digit_removal_from_string(string_with_digits):
    assert string_with_digits.translate(remove_digits()) == 'trial'

class TestToNumeric:

    def test_dash_to_NaN(self):
        new_series = to_numeric(pd.Series(["-", 1, "-"]))
        assert new_series.isna().sum() == 2
        assert is_numeric_dtype(new_series)

    @pytest.mark.parametrize("numeric_values", ([1, 1, 2], [10, 100, 1000], [1.0, 1.1, 1.2], [1, -1, 2.0]))
    def test_remove_non_numeric_characters(self, numeric_values):
        new_series = to_numeric(pd.Series([
            "-", numeric_values[0], f"{numeric_values[1]}a", f"{numeric_values[2]}-"
        ]))
        assert np.allclose(new_series.iloc[1:], numeric_values)
        assert is_numeric_dtype(new_series)

    def test_no_change_to_numeric_characters(self):
        original_series = pd.Series([0, 1, 100, 2.0, np.nan])
        new_series = to_numeric(original_series.copy())
        assert new_series.equals(original_series)
