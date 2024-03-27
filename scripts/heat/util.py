import numpy as np
import pandas as pd
import pycountry
import pytz


def get_alpha2(country, eurostat=True):
    if country in ["United Kingdom", "GB", "GBR"] and eurostat is True:
        return "UK"
    elif country in ["Greece", "GR", "GRC"] and eurostat is True:
        return "EL"
    else:
        return pycountry.countries.lookup(country).alpha_2


def get_alpha3(country):
    if country == "UK":
        country = "GB"
    elif country == "EL":
        country = "GR"
    return pycountry.countries.lookup(country).alpha_3


def to_numeric(series):
    series = series.astype(str).str.extract(r"(\-*\d+\.*\d*)")[0]
    return pd.to_numeric(series, errors="coerce")


def pj_to_twh(array):
    return array / 3.6


def tj_to_twh(array):
    return pj_to_twh(array) * 1e-3


def gwh_to_tj(array):
    return array * 3.6


def ktoe_to_twh(array):
    return array * 1.163e-2


def tj_to_ktoe(array):
    return array * 23.88e-3


def update_timeseries_timezone(x, country, model_year):
    """
    Shift a generic profile forward/backward in time based on a country's timezone
    """
    if country == "UK":
        country = "GB"
    elif country == "EL":
        country = "GR"
    tz = pytz.country_timezones[country][0]
    try:
        idx = x.index.tz_localize(tz, nonexistent="shift_forward").tz_convert("UTC")
    except pytz.AmbiguousTimeError as err:
        idx = x.index.tz_localize(
            tz, ambiguous=x.index != err.args[0], nonexistent="shift_forward"
        ).tz_convert("UTC")
    shift = len(idx[idx.year > model_year]) - len(idx[idx.year < model_year])

    x = np.roll(x, shift=shift)

    return x


def read_tdf(filename):
    df = pd.read_csv(filename, header=0)
    tdf = df.set_index([i for i in df.columns[:-1]]).squeeze()
    return tdf


def read_eurostat_tsv(path_to_tsv, index_names, slice_idx=None, slice_lvl=None):
    df = pd.read_csv(path_to_tsv, delimiter="\t", index_col=0)
    df.index = df.index.str.split(",", expand=True).rename(index_names)
    if slice_idx is not None:
        df = df.xs(slice_idx, level=slice_lvl)
    df.columns = df.columns.astype(int).rename("year")
    return df.apply(to_numeric)


def get_timedelta(model_time, model_year):
    model_timedelta = (
        pd.to_datetime(model_time[1].replace("year", model_year))
        - pd.to_datetime(model_time[0].replace("year", model_year))
    ).days + 1
    if pd.to_datetime(model_year).is_leap_year:
        model_timedelta = model_timedelta / 366
    else:
        model_timedelta = model_timedelta / 365
    return model_timedelta


def verify_profiles(profile, key, annual_demand):
    assert (profile <= 0).all().all()
    if isinstance(key, list):
        demand = sum([
            annual_demand.xs(_k, level="end_use").sum(level="id") for _k in key
        ])
    else:
        demand = annual_demand.xs(key, level="end_use").sum(level="id")
    assert np.allclose(profile.sum().abs().reindex(demand.index), demand.abs())


def filter_small_values(data, rel_tol=1e-5):
    """
    Ignore any values in a dataset that are rel_tol smaller than the maximum value
    `data` can be a pandas Series or DataFrame. If DataFrame, each column will be dealt
    with independently.
    """
    data_sum = data.sum()
    if data_sum == 0 and data.min() >= 0:
        return data

    data[abs(data) < abs(data).max() * rel_tol] = 0

    # Resulting sum of deviation should be smaller or equal to the rel_tol, as a crude
    # measure that data has varied negligibly
    try:
        assert all(1 - abs(data.sum() / data_sum) <= rel_tol)
    except TypeError:
        assert 1 - abs(data.sum() / data_sum) <= rel_tol

    return data
