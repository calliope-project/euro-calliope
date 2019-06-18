"""Link all those locations that are neighbours."""
import itertools

import geopandas as gpd
import shapely
import jinja2
import networkx as nx

EPSG_3035_PROJ4 = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs "
K_EDGE_CONNECTION = 1 # every component of the graph is connect with at least this amount of edges
TEMPLATE = """
links:
    {% for location1, location2 in graph.edges %}
    {{ location1 | replace(".", "-") }},{{ location2 | replace(".", "-") }}:
        techs:
            ac_transmission:
    {% endfor %}
"""


def link_neighbours(path_to_locations, path_to_result):
    """Link all those locations that are neighbours.

    If that creates unconnected islands, connect the islands with the shortest
    connection possible.
    """
    graph = _create_graph(path_to_locations)
    if len(graph) >= 2:
        graph = _connect_graph(graph)
    links = jinja2.Template(TEMPLATE).render(
        graph=graph
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(links)


def _create_graph(path_to_locations):
    regions = gpd.read_file(path_to_locations).to_crs(EPSG_3035_PROJ4).set_index("id")
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


def _connect_graph(graph):
    # A graph based on direct neighbouring has disconnected components, for example UK.
    # Hence, here the components are connected.
    all_possible_edges = [
        (loc1, loc2, graph.nodes[loc1]["centroid"].distance(graph.nodes[loc2]["centroid"]))
        for loc1, loc2 in itertools.product(graph.nodes, graph.nodes)
    ]
    connecting_edges = nx.k_edge_augmentation(graph, K_EDGE_CONNECTION, avail=all_possible_edges)
    graph.add_edges_from(connecting_edges)
    return graph


if __name__ == "__main__":
    link_neighbours(
        path_to_locations=snakemake.input.units,
        path_to_result=snakemake.output[0]
    )
