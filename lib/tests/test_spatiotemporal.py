import pytest
from shapely.geometry import box
import pandas as pd
import geopandas as gpd


from eurocalliopelib import area_weighted_time_series


WGS84 = "EPSG:4326"
EPSG3035 = "EPSG:3035"
EPSG27700 = "EPSG:27700"


@pytest.fixture
def single_shape():
    return gpd.GeoDataFrame(geometry=[box(0, 0, 4, 4)], crs=EPSG3035)


@pytest.fixture
def single_shape_small():
    return gpd.GeoDataFrame(geometry=[box(0, 0, 3, 4)], crs=EPSG3035)


@pytest.fixture
def four_shapes():
    return gpd.GeoDataFrame(
        geometry=[box(0, 0, 2, 2), box(0, 2, 2, 4), box(2, 2, 4, 4), box(2, 0, 4, 2)],
        crs=EPSG3035
    )


@pytest.fixture
def make_timeseries():
    def _make_timeseries(value):
        return pd.Series(
            index=pd.period_range(start="2021-04-18 12:00", periods=2, freq="D", name="timestep"),
            data=value,
            name="value"
        )
    return _make_timeseries


@pytest.fixture
def make_spatiotemporal_data(make_timeseries):
    def _make_spatiotemporal_data(values):
        assert len(values) == 4
        da = (
            pd
            .concat([
                make_timeseries(values[0]).to_frame().assign(x=1, y=1),
                make_timeseries(values[1]).to_frame().assign(x=1, y=3),
                make_timeseries(values[2]).to_frame().assign(x=3, y=3),
                make_timeseries(values[3]).to_frame().assign(x=3, y=1)
            ])
            .reset_index()
            .set_index(["y", "x", "timestep"])
            .to_xarray()
            ["value"]
        )
        da.attrs["crs"] = EPSG3035
        return da
    return _make_spatiotemporal_data


@pytest.fixture
def non_rectangular_spatiotemporal_data(make_timeseries):
    da = (
        pd
        .concat([
            make_timeseries(8).to_frame().assign(y=1, x=1),
            make_timeseries(8).to_frame().assign(y=1, x=3),
            make_timeseries(8).to_frame().assign(y=3, x=1),
            make_timeseries(8).to_frame().assign(y=3, x=3),
            make_timeseries(8).to_frame().assign(y=5, x=3)
        ])
        .reset_index()
        .set_index(["y", "x", "timestep"])
        .to_xarray()
        ["value"]
    )
    da.attrs["crs"] = EPSG3035
    return da


@pytest.fixture
def spatiotemporal_data_with_nodata(make_timeseries):
    da = (
        pd
        .concat([
            make_timeseries(3).to_frame().assign(y=1, x=1),
            make_timeseries(9).to_frame().assign(y=1, x=3),
            make_timeseries(9).to_frame().assign(y=3, x=1)
        ])
        .reset_index()
        .set_index(["y", "x", "timestep"])
        .to_xarray()
        ["value"]
    )
    da.attrs["crs"] = EPSG3035
    return da


@pytest.fixture
def any_spatiotemporal_data(make_spatiotemporal_data):
    return make_spatiotemporal_data([8, 8, 8, 8])


def test_fails_without_crs_in_shapes(single_shape, any_spatiotemporal_data):
    single_shape.crs = None
    with pytest.raises(AssertionError):
        area_weighted_time_series(shapes=single_shape, spatiotemporal=any_spatiotemporal_data)


def test_fails_without_crs_in_spatiotemporal(single_shape, any_spatiotemporal_data):
    del any_spatiotemporal_data.attrs["crs"]
    with pytest.raises(AssertionError):
        area_weighted_time_series(shapes=single_shape, spatiotemporal=any_spatiotemporal_data)


@pytest.mark.parametrize('crs1,crs2', [
    (EPSG3035, EPSG27700),
    (WGS84, EPSG3035),
    (WGS84, EPSG27700)
])
def test_fails_with_diverging_crs(single_shape, any_spatiotemporal_data, crs1, crs2):
    single_shape = single_shape.to_crs(crs1)
    any_spatiotemporal_data.attrs["crs"] = crs2
    with pytest.raises(AssertionError):
        area_weighted_time_series(shapes=single_shape, spatiotemporal=any_spatiotemporal_data)


@pytest.mark.parametrize('crs', [WGS84, EPSG3035, EPSG27700])
def test_allows_all_crs(single_shape, any_spatiotemporal_data, crs):
    single_shape = single_shape.set_crs(crs, allow_override=True)
    any_spatiotemporal_data.attrs["crs"] = crs
    area_weighted_time_series(shapes=single_shape, spatiotemporal=any_spatiotemporal_data)


@pytest.mark.parametrize('correct_name,wrong_name', [
    ("y", "lat"),
    ("x", "lon"),
    ("timestep", "timesteps")
])
def test_fails_with_wrong_dimension_names(single_shape, any_spatiotemporal_data, correct_name, wrong_name):
    any_spatiotemporal_data = any_spatiotemporal_data.rename({correct_name: wrong_name})
    with pytest.raises(AssertionError):
        area_weighted_time_series(shapes=single_shape, spatiotemporal=any_spatiotemporal_data)


def test_equal_values(single_shape, make_spatiotemporal_data):
    spatiotemporal_data = make_spatiotemporal_data([8, 8, 8, 8])
    weighted_ts = area_weighted_time_series(shapes=single_shape, spatiotemporal=spatiotemporal_data)
    assert (weighted_ts.iloc[:, 0].values == 8).all()


def test_not_equal_values(single_shape, make_spatiotemporal_data):
    spatiotemporal_data = make_spatiotemporal_data([8, 8, 8, 0])
    weighted_ts = area_weighted_time_series(shapes=single_shape, spatiotemporal=spatiotemporal_data)
    expected_mean = (8 + 8 + 8 + 0) / 4
    assert (weighted_ts.iloc[:, 0].values == expected_mean).all()


def test_handles_non_rectangular_data(single_shape, non_rectangular_spatiotemporal_data):
    weighted_ts = area_weighted_time_series(
        shapes=single_shape,
        spatiotemporal=non_rectangular_spatiotemporal_data
    )
    assert (weighted_ts.iloc[:, 0].values == 8).all()


def test_four_shapes(four_shapes, make_spatiotemporal_data):
    spatiotemporal_data = make_spatiotemporal_data([0, 1, 2, 3])
    weighted_ts = area_weighted_time_series(shapes=four_shapes, spatiotemporal=spatiotemporal_data)
    assert (weighted_ts.iloc[:, 0].values == 0).all()
    assert (weighted_ts.iloc[:, 1].values == 1).all()
    assert (weighted_ts.iloc[:, 2].values == 2).all()
    assert (weighted_ts.iloc[:, 3].values == 3).all()


def test_handles_nodata_correctly(single_shape, spatiotemporal_data_with_nodata):
    weighted_ts = area_weighted_time_series(shapes=single_shape, spatiotemporal=spatiotemporal_data_with_nodata)
    expected_value = (3 + 9 + 9) / 3
    assert (weighted_ts.iloc[:, 0].values == expected_value).all()


@pytest.mark.xfail(reason="The internal method does not have the required resolution.")
def test_partial_match(single_shape_small, make_spatiotemporal_data):
    spatiotemporal_data = make_spatiotemporal_data([1, 1, 2, 2])
    weighted_ts = area_weighted_time_series(shapes=single_shape_small, spatiotemporal=spatiotemporal_data)
    expected_mean = 1 * (8 / 12) + 2 * (4 / 12)
    assert weighted_ts.iloc[:, 0].values.mean() == expected_mean
