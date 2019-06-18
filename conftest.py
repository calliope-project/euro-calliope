import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--include-regional-resolution",
        action="store_true",
        default=False,
        help="run tests on the regional resolution"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "regional: mark test as running on regional resolution")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--include-regional-resolution"):
        return
    skip_regional = pytest.mark.skip(reason="Use command line flag to run on regional resolution.")
    for item in items:
        if "regional" in item.keywords:
            item.add_marker(skip_regional)
