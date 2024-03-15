import numpy as np
import pandas as pd
import pytest

from scripts.transmission.template_transmission_entsoe_tyndp import (
    _split_links_in_index,
)


@pytest.fixture
def sliced_df():
    def _sliced_df(vals, idx):
        return pd.DataFrame(data=vals, index=idx, columns=["Value"])

    return _sliced_df


@pytest.mark.parametrize(
    ("input_idx", "expected_idx", "input_values", "expected_values"),
    (
        [["UK00-IE00", "UK00-FR00"], [("GBR", "IRL"), ("GBR", "FRA")], [1, 2], [1, 2]],
        [["ITN0-CH00", "ITN1-CH00"], [("ITA", "CHE")], [1, 2], [3]],
        [["ITN0-CH00", "ITN0-ITN1", "ITN1-CH00"], [("ITA", "CHE")], [1, 2, 3], [4]],
    ),
)
def test_split_links_in_index(
    sliced_df, input_idx, expected_idx, input_values, expected_values
):
    input_df = sliced_df(input_values, input_idx)
    # Links are split and then summed over all those corresponding to a specific country.
    # All links that are internal to a country are ignored.
    # Expected index/values are based on this expected summation and removal of internal links.
    cleaned_series = _split_links_in_index(input_df)
    assert cleaned_series.index.symmetric_difference(expected_idx).empty
    assert np.allclose(cleaned_series.loc[expected_idx].values, expected_values)
