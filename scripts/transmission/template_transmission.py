"""Link all those locations that are neighbours."""
import geopandas as gpd
import shapely
import networkx as nx

from eurocalliopelib.geo import EPSG3035
from eurocalliopelib.parametrise_template import parametrise_template

K_EDGE_CONNECTION = 1 # every component of the graph is connect with at least this amount of edges


def construct_techs_and_links(path_to_units, path_to_output, scaling_factors, sea_connections, path_to_template):
    """Link all those locations that are neighbours."""
    graph = _create_graph_with_land_connections(path_to_units)
    if sea_connections:
        graph.add_edges_from([(loc1, loc2) for loc1, loc2 in sea_connections])

    assert nx.is_connected(graph), "There are electrical islands in the network graph."
    return parametrise_template(path_to_template, path_to_output, graph=graph, scaling_factors=scaling_factors)


def _create_graph_with_land_connections(path_to_units):
    regions = gpd.read_file(path_to_units).to_crs(EPSG3035).set_index("id")
    graph = nx.Graph()
    graph.add_nodes_from([
        (region_id, {"centroid": region_centroid})
        for region_id, region_centroid in regions.centroid.to_dict().items()
    ])
    neighbours = {
        index: _neighbours(region.geometry, index, regions) for index, region in regions.iterrows()
    }
    graph.add_edges_from([
        (region, other)
        for region, all_others in neighbours.items()
        for other in all_others
    ])
    return graph


def _neighbours(region_geometry, region_index, regions):
    region_geometry = shapely.prepared.prep(region_geometry)
    return [other_index for other_index, other_region in regions.iterrows()
            if (other_index is not region_index) and region_geometry.intersects(other_region.geometry)]


if __name__ == "__main__":
    construct_techs_and_links(
        path_to_template=snakemake.input.template,
        path_to_units=snakemake.input.units,
        scaling_factors=snakemake.params.scaling_factors,
        sea_connections=snakemake.params.sea_connections,
        path_to_output=snakemake.output[0]
    )
