from textwrap import dedent

import geopandas as gpd
import matplotlib.pyplot as plt

from eurocalliopelib.geo import EPSG3035

LIGHT_BLUE = "#8fa2cf"
EDGE_WIDTH = 0.2
WHITE = "white"

MAP_MIN_X = 2200000
MAP_MIN_Y = 1400000
MAP_MAX_X = 6300000
MAP_MAX_Y = 5500000


def spatial_scope_and_resolutions(path_to_regional_units, path_to_national_units, dpi, path_to_output):
    # read data
    regional = gpd.read_file(path_to_regional_units).to_crs(EPSG3035).simplify(10000)
    national = gpd.read_file(path_to_national_units).to_crs(EPSG3035).simplify(10000)

    # prepare figure
    plt.rcParams.update({'font.size': 6})
    fig = plt.figure(figsize=(4, 1.4))
    axes = fig.subplots(1, 4, gridspec_kw={'width_ratios': [0.15, 0.35, 0.35, 0.15]})
    for ax in axes:
        ax.axis("off")

    # plot
    _plot_text(national, regional, ax=axes[1])
    _plot_map(regional, ax=axes[2], edge_width=0.1)
    _plot_map(national, ax=axes[2], face_color="none", edge_width=0.3)

    # adjust figure
    fig.subplots_adjust(wspace=0, hspace=0, top=1, bottom=0, left=0, right=1)
    fig.savefig(path_to_output, dpi=dpi, transparent=True)


def _plot_text(national, regional, ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(
        x=0.0,
        y=0.65,
        fontdict={'fontweight': 'semibold', 'va': 'top', 'ha': 'left'},
        s="Euro-Calliope's resolutions:"
    )
    ax.text(
        x=0.03,
        y=0.63,
        fontdict={'va': 'top', 'ha': 'left'},
        s=dedent(
            f"""
            · a single continent,
            · {len(national)} countries,
            · or {len(regional)} regions.
            """
        )
    )
    ax.text(
        x=0.0,
        y=0.15,
        fontdict={'va': 'top', 'size': 4},
        s="Source: GADM, NUTS © EuroGeographics."
    )


def _plot_map(data, ax, edge_width, face_color=LIGHT_BLUE, edge_color=WHITE):
    ax.set_aspect('equal')
    data.plot(
        ax=ax,
        facecolor=face_color,
        edgecolor=edge_color,
        linewidth=edge_width
    )
    ax.set_xlim(MAP_MIN_X, MAP_MAX_X)
    ax.set_ylim(MAP_MIN_Y, MAP_MAX_Y)


if __name__ == "__main__":
    spatial_scope_and_resolutions(
        path_to_regional_units=snakemake.input.regional_units,
        path_to_national_units=snakemake.input.national_units,
        dpi=snakemake.params.dpi,
        path_to_output=snakemake.output[0]
    )
