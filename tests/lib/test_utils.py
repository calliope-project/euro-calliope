import pytest
import io

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
from eurocalliopelib.utils import convert_country_code, read_eurostat_tsv, remove_digits, add_idx_level, to_numeric, convert_unit, UNIT_CONVERSION_MAPPING


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

    @pytest.mark.parametrize(
        ("raw", "expected"), [
            ([1, "1a", "2-"], [1, 1, 2]),
            ([10, "100a", "1000-"], [10, 100, 1000]),
            ([1.0, "1.0a", "-1-"], [1.0, 1.0, -1])
        ]
    )
    def test_remove_non_numeric_characters(self, raw, expected):
        new_series = to_numeric(pd.Series(raw))
        assert np.allclose(new_series, expected)
        assert is_numeric_dtype(new_series)

    def test_no_change_to_numeric_characters(self):
        original_series = pd.Series([0, 1, 100, 2.0, np.nan])
        new_series = to_numeric(original_series.copy())
        pd.testing.assert_series_equal(new_series, original_series)


class TestConvertUnit:
    @pytest.fixture
    def create_df(self):
        def _create_df(input_units):
            index = pd.MultiIndex.from_product([["foo", "bar"], input_units], names=["baz", "unit"])
            return pd.Series(index=index, data=1)
        return _create_df

    @pytest.fixture
    def create_idx(self):
        def _create_idx(multi_or_single):
            if multi_or_single == "multi":
                return pd.MultiIndex.from_tuples([("foo", "bar")], names=["baz", "foobar"])
            elif multi_or_single == "single":
                return pd.Index(["foo"], name="baz")
        return _create_idx

    def get_conversions():
        return [(val, *from_to) for from_to, val in UNIT_CONVERSION_MAPPING.items()]

    @pytest.mark.parametrize(("conversion_val", "convert_from", "convert_to"), get_conversions())
    @pytest.mark.parametrize("index_type", ("multi", "single"))
    def test_no_unit_in_index_in_or_out(self, conversion_val, convert_from, convert_to, create_idx, index_type):
        df = pd.DataFrame([1], index=create_idx(index_type))
        df_converted = convert_unit(df, convert_to, convert_from, unit_in_output_idx=False)
        pd.testing.assert_frame_equal(df_converted, df * conversion_val)

    @pytest.mark.parametrize(("conversion_val", "convert_from", "convert_to"), get_conversions())
    @pytest.mark.parametrize("index_type", ("multi", "single"))
    def test_no_unit_in_index_in(self, conversion_val, convert_from, convert_to, create_idx, index_type):
        df = pd.DataFrame([1], index=create_idx(index_type))
        df_converted = convert_unit(df, convert_to, convert_from, unit_in_output_idx=True)
        assert "unit" in df_converted.index.names
        assert df_converted.index.get_level_values("unit").difference([convert_to]).empty

        pd.testing.assert_frame_equal(df_converted.droplevel("unit"), df * conversion_val)

    @pytest.mark.parametrize(("unit_in", "unit_out"), [(["TJ"], "TWh"), (["TJ"], "twh"), (["tj"], "TWh"), (["tj", "TJ"], "twh")])
    def test_upper_lower(self, unit_in, unit_out, create_df):
        df = create_df(unit_in)
        df_converted = convert_unit(df, unit_out, unit_in_output_idx=True)
        assert df_converted.index.get_level_values("unit").difference([unit_out.lower()]).empty

    @pytest.mark.parametrize("index_type", ["multi", "single"])
    def test_same_unit_in_and_out_no_unit_in_idx(self, create_idx, index_type):
        df = pd.DataFrame([1], index=create_idx(index_type))
        df_converted = convert_unit(df, "twh", "twh", unit_in_output_idx=False)
        pd.testing.assert_frame_equal(df_converted, df)

    @pytest.mark.parametrize("unit_in", ["TWh", "twh"])
    def test_same_unit_in_and_out_and_unit_in_index_in_and_out(self, unit_in, create_df):
        df = create_df([unit_in])
        df_converted = convert_unit(df, "twh", unit_in, unit_in_output_idx=True)
        assert "unit" in df_converted.index.names
        assert df_converted.index.get_level_values("unit").difference(["twh"]).empty

        pd.testing.assert_series_equal(df_converted.droplevel("unit"), df.droplevel("unit"))

    @pytest.mark.parametrize("unit_in", ["TWh", "twh"])
    def test_same_unit_in_and_out_and_unit_in_index_out(self, unit_in):
        df = pd.DataFrame([1], index=pd.Index(["foo"], name="baz"))
        df_converted = convert_unit(df, "twh", unit_in, unit_in_output_idx=True)
        assert "unit" in df_converted.index.names
        assert df_converted.index.get_level_values("unit").difference(["twh"]).empty

        pd.testing.assert_frame_equal(df_converted.droplevel("unit"), df)

    @pytest.mark.parametrize(("conversion_val", "convert_from", "convert_to"), get_conversions())
    def test_unit_in_index_in_not_out(self, conversion_val, convert_from, convert_to, create_df):
        df = create_df([convert_from])
        df_converted = convert_unit(df, convert_to, convert_from, unit_in_output_idx=False)
        assert "unit" not in df_converted.index.names
        pd.testing.assert_series_equal(df_converted, df.droplevel("unit") * conversion_val)

    @pytest.mark.parametrize(("conversion_val", "convert_from", "convert_to"), get_conversions())
    def test_unit_in_index_in_and_out(self, conversion_val, convert_from, convert_to, create_df):
        df = create_df([convert_from])
        df_converted = convert_unit(df, convert_to, convert_from, unit_in_output_idx=True)
        assert "unit" in df_converted.index.names
        assert df_converted.index.get_level_values("unit").difference([convert_to]).empty

        pd.testing.assert_series_equal(df_converted.droplevel("unit"), df.droplevel("unit") * conversion_val)

    @pytest.mark.parametrize(("conversion_val", "convert_from", "convert_to"), get_conversions())
    def test_infer_unit_in(self, conversion_val, convert_from, convert_to, create_df):
        df = create_df([convert_from])
        df_converted = convert_unit(df, convert_to, unit_in_output_idx=True)
        assert "unit" in df_converted.index.names
        assert df_converted.index.get_level_values("unit").difference([convert_to]).empty

        pd.testing.assert_series_equal(df_converted.droplevel("unit"), df.droplevel("unit") * conversion_val)

    @pytest.mark.parametrize(("conversion_val", "convert_from", "convert_to"), get_conversions())
    def test_infer_unit_in_from_multiple(self, conversion_val, convert_from, convert_to, create_df):
        df = create_df([convert_from, "foo"])
        df_converted = convert_unit(df, convert_to, convert_from, unit_in_output_idx=True)
        assert "unit" in df_converted.index.names
        assert df_converted.index.get_level_values("unit").difference([convert_to, "foo"]).empty

        pd.testing.assert_series_equal(df_converted.xs(convert_to, level="unit"), df.xs(convert_from, level="unit") * conversion_val, check_dtype=False)
        pd.testing.assert_series_equal(df_converted.xs("foo", level="unit"), df.xs("foo", level="unit"), check_dtype=False)

    @pytest.mark.parametrize("index_type", ["multi", "single"])
    def test_cannot_infer_unit_without_idx_level(self, create_idx, index_type):
        df = pd.DataFrame([1], index=create_idx(index_type))
        with pytest.raises(ValueError, match="Cannot infer unit for data"):
            convert_unit(df, "twh")

    def test_cannot_infer_unit_from_multiple_in(self, create_df):
        df = create_df(["tj", "foo"])
        with pytest.raises(AssertionError, match="Cannot infer unit for data with multiple available units"):
            convert_unit(df, "twh")

    def test_cannot_remove_unit_from_multiple_in(self, create_df):
        df = create_df(["tj", "foo"])
        with pytest.raises(AssertionError, match="Cannot drop the index level"):
            convert_unit(df, "twh", input_unit="tj", unit_in_output_idx=False)

