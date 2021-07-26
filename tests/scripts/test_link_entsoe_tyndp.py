import pytest
import pandas as pd
import numpy as np

from scripts.link_entsoe_tyndp import _split_links_in_index

@pytest.fixture
def sliced_df():
    def _sliced_df(idx):
        data = list(range(10, 10 + len(idx)))
        return pd.DataFrame(data=data, index=idx, columns=["Value"])
    return _sliced_df

@pytest.mark.parametrize(
    ("input_idx", "expected_idx", "expected_values"), (
        [["UK00-IE00", "UK00-FR00"], [("GBR", "IRL"), ("GBR", "FRA")], [10, 11]],
        [["ITN0-CH00", "ITN1-CH00"], [("ITA", "CHE")], [21]],
        [["ITN0-CH00", "ITN1-CH00", "ITN0-ITN1"], [("ITA", "CHE")], [21]],
    )
)
def test_split_links_in_index(sliced_df, input_idx, expected_idx, expected_values):
    input_df = sliced_df(input_idx)
    cleaned_series = _split_links_in_index(input_df)
    assert cleaned_series.index.symmetric_difference(expected_idx).empty
    assert np.allclose(cleaned_series.loc[expected_idx].values, expected_values)
